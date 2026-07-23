test_that("fit_recursive_mlogit recovers a seeded 3-category system", {
  skip_on_cran()
  set.seed(11)
  n <- 600; Ks <- 2; Kh <- 2
  # truth: mlogit cats 1,2 vs base 0; loading l2 = 0.5 (l1 = 1 fixed)
  g1 <- c(-0.5, 0.8, 0.3); g2 <- c(-1.0, 0.3, 0.6); l2 <- 0.5
  a1 <- -0.5; a2t <- -1.0; b2x <- 0.3; b20 <- -1.0
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig); u <- uv[,1]; v <- uv[,2]
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  e1 <- g1[1] + g1[2]*zs + g1[3]*xs + 1  * v[id1]
  e2 <- g2[1] + g2[2]*zs + g2[3]*xs + l2 * v[id1]
  den <- 1 + exp(e1) + exp(e2)
  r <- runif(n*Ks); p1 <- exp(e1)/den; p2 <- exp(e2)/den
  y1 <- ifelse(r < p1, 1L, ifelse(r < p1 + p2, 2L, 0L))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = factor(y1, 0:2))
  id2 <- rep(1:n, each = Kh)
  xh <- rnorm(n*Kh); dur <- rexp(n*Kh, 1) + 0.5
  y1h <- y1[match(id2, id1) + sample(0:(Ks-1), n*Kh, TRUE)]
  ev <- rpois(n*Kh, dur * exp(b20 + a1*(y1h==1) + a2t*(y1h==2) + b2x*xh + u[id2]))
  dhaz <- data.frame(id = id2, x = xh, y1 = factor(y1h, 0:2), dur, ev)

  fit <- fit_recursive_mlogit(y1 ~ z + x, ev ~ y1 + x,
                              data_sel = dsel, data_haz = dhaz,
                              id = "id", exposure = "dur", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  for (p in c("haz:y11", "haz:y12", "sel1:z", "sel2:z", "load v:2")) {
    truth <- c(`haz:y11` = a1, `haz:y12` = a2t, `sel1:z` = g1[2],
               `sel2:z` = g2[2], `load v:2` = l2)[p]
    expect_lt(abs(co[p] - truth), 1.96 * se[p])
  }
  s <- recover_structural(fit, "`haz:y12` - `haz:y11`")
  expect_equal(s$estimate, unname(co["haz:y12"] - co["haz:y11"]), tolerance = 1e-8)
})

test_that("split_episodes partitions exposure and assigns the event correctly", {
  d <- data.frame(id = 1:3, dur = c(2, 5, 20), ev = c(1, 0, 1), x = 7:9)
  s <- split_episodes(d, event = "ev", exposure = "dur", breaks = c(3, 12))
  # spell 1 (dur 2, ev 1): one segment [0,3), exposure 2, event 1
  s1 <- s[s$id == 1, ]; expect_equal(nrow(s1), 1L)
  expect_equal(s1$dur, 2); expect_equal(s1$ev, 1)
  # spell 2 (dur 5, ev 0): segments [0,3)+[3,12), exposures 3+2, events 0
  s2 <- s[s$id == 2, ]; expect_equal(sort(s2$dur), c(2, 3))
  expect_equal(sum(s2$ev), 0)
  # spell 3 (dur 20, ev 1): three segments, exposures 3+9+8, event only in last
  s3 <- s[s$id == 3, ][order(s[s$id == 3, ]$.t0), ]
  expect_equal(s3$dur, c(3, 9, 8)); expect_equal(s3$ev, c(0, 0, 1))
  expect_equal(sum(s$dur), sum(d$dur))     # exposure conserved
  expect_equal(nlevels(s$interval), 3L)
})
