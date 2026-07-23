test_that("lognormal AFT kernel recovers seeded sigma and alpha", {
  skip_on_cran()
  set.seed(53)
  n <- 600; Ks <- 2; Kh <- 2
  a2 <- -0.5; b2x <- 0.3; b20 <- -0.3; sg <- 0.7
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig)
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(-0.3 + 0.7*zs + 0.4*xs + uv[id1, 2]))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = y1)
  id2 <- rep(1:n, each = Kh)
  xh <- rnorm(n*Kh)
  y1h <- y1[match(id2, id1) + sample(0:(Ks-1), n*Kh, TRUE)]
  lp <- b20 + a2*y1h + b2x*xh + uv[id2, 1]
  Tt <- exp(lp + sg * rnorm(n*Kh))               # log T ~ N(lp, sg^2)
  cens <- 3; d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, t = tt, d = d)

  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x, dsel, dhaz, id = "id",
                       family = "lognormal", time = "t", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log sigma"] - log(sg)), 1.96 * se["log sigma"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
  # entry = 0 must reproduce the no-truncation fit exactly
  dhaz$t0 <- 0
  fit0 <- fit_recursive(y1 ~ z + x, d ~ y1 + x, dsel, dhaz, id = "id",
                        family = "lognormal", time = "t", entry = "t0", Q = 9)
  expect_equal(unname(fit0$coef["haz:y1"]), unname(co["haz:y1"]),
               tolerance = 1e-6)
})
