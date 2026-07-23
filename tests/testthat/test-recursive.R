test_that("fit_recursive recovers a seeded recursive mixed-family system", {
  skip_on_cran()
  set.seed(7)
  n <- 600; Ks <- 2; Kh <- 2
  a2 <- -0.5; b2x <- 0.3; b20 <- -1.0
  g_z <- 0.7; g_x <- 0.4; g0 <- -0.3
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig); u <- uv[,1]; v <- uv[,2]
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(g0 + g_z*zs + g_x*xs + v[id1]))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = y1)
  id2 <- rep(1:n, each = Kh)
  xh <- rnorm(n*Kh); dur <- rexp(n*Kh, 1) + 0.5
  y1h <- y1[match(id2, id1) + sample(0:(Ks-1), n*Kh, TRUE)]
  ev <- rpois(n*Kh, dur * exp(b20 + a2*y1h + b2x*xh + u[id2]))
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, dur = dur, ev = ev)

  fit <- fit_recursive(y1 ~ z + x, ev ~ y1 + x,
                       data_sel = dsel, data_haz = dhaz,
                       id = "id", exposure = "dur", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  # truth inside 95% CI for the structural targets
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"]  - g_z), 1.96 * se["sel:z"])
  expect_lt(abs(co["haz:x"]  - b2x), 1.96 * se["haz:x"])
  # RE variances: truth inside the 95% CI on the log-sd scale (the
  # parameterisation actually estimated), same criterion as the betas
  expect_lt(abs(co["log s1"] - log(sqrt(0.4))), 1.96 * se["log s1"])
  expect_lt(abs(co["log s2"] - log(sqrt(0.8))), 1.96 * se["log s2"])
  # recover_structural works on the fit (list coef/vcov path)
  s <- recover_structural(fit, "`haz:y1` * 2")
  expect_equal(s$estimate, unname(2 * co["haz:y1"]), tolerance = 1e-8)
})
