test_that("three-level fit recovers seeded child-frailty variance and alpha", {
  skip_on_cran()
  set.seed(61)
  n <- 400; Kc <- 2; Ks <- 2; Ksp <- 2      # mothers, children, sel-obs, spells
  a2 <- -0.5; b2x <- 0.3; b20 <- -1.0
  sw2 <- 0.3                                # child-frailty variance (target)
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig); u <- uv[,1]; v <- uv[,2]
  # selection: mother-level RE v
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(-0.3 + 0.7*zs + 0.4*xs + v[id1]))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = y1)
  # hazard: spells in children in mothers; lp = ... + u_mother + w_child
  mo <- rep(1:n, each = Kc * Ksp)
  ch <- rep(rep(1:Kc, each = Ksp), times = n)
  w_child <- rnorm(n * Kc, 0, sqrt(sw2))
  chix <- (mo - 1) * Kc + ch
  xh <- rnorm(n * Kc * Ksp); dur <- rexp(n * Kc * Ksp, 1) + 0.5
  y1h <- y1[match(mo, id1) + sample(0:(Ks-1), n*Kc*Ksp, TRUE)]
  ev <- rpois(n*Kc*Ksp, dur * exp(b20 + a2*y1h + b2x*xh + u[mo] + w_child[chix]))
  dhaz <- data.frame(id = mo, bid = ch, x = xh, y1 = y1h, dur, ev)

  fit <- fit_recursive(y1 ~ z + x, ev ~ y1 + x, dsel, dhaz, id = "id",
                       family = "poisson", exposure = "dur",
                       nested = "bid", Q = 7, Qw = 5)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log s_w(nested)"] - log(sqrt(sw2))),
            1.96 * se["log s_w(nested)"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["log s1"] - log(sqrt(0.4))), 1.96 * se["log s1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
  expect_equal(fit$sw2, unname(exp(co["log s_w(nested)"])^2))
  s <- recover_structural(fit, "`haz:y1` * 2")
  expect_equal(s$estimate, unname(2 * co["haz:y1"]), tolerance = 1e-8)
})
