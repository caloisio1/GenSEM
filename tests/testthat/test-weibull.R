test_that("fit_recursive with family='weibull' recovers seeded shape and alpha", {
  skip_on_cran()
  set.seed(31)
  n <- 600; Ks <- 2; Kh <- 2
  a2 <- -0.5; b2x <- 0.3; b20 <- -1.0; p_shape <- 1.4
  g_z <- 0.7; g_x <- 0.4; g0 <- -0.3
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig); u <- uv[,1]; v <- uv[,2]
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(g0 + g_z*zs + g_x*xs + v[id1]))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = y1)
  id2 <- rep(1:n, each = Kh)
  xh <- rnorm(n*Kh)
  y1h <- y1[match(id2, id1) + sample(0:(Ks-1), n*Kh, TRUE)]
  lam <- exp(b20 + a2*y1h + b2x*xh + u[id2])          # Weibull PH rate
  Tt <- (rexp(n*Kh) / lam)^(1/p_shape)                # S(t)=exp(-lam t^p)
  cens <- 2.5
  d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, t = tt, d = d)

  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x,
                       data_sel = dsel, data_haz = dhaz,
                       id = "id", family = "weibull", time = "t", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log shape"] - log(p_shape)), 1.96 * se["log shape"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - g_z), 1.96 * se["sel:z"])
  expect_lt(abs(co["log s1"] - log(sqrt(0.4))), 1.96 * se["log s1"])
  expect_equal(fit$ancillary, unname(exp(co["log shape"])))
  # left truncation path runs and shifts nothing structurally on t0=0 data
  dhaz$t0 <- 0
  fit0 <- fit_recursive(y1 ~ z + x, d ~ y1 + x, dsel, dhaz, id = "id",
                        family = "weibull", time = "t", entry = "t0", Q = 9)
  expect_equal(unname(fit0$coef["haz:y1"]), unname(co["haz:y1"]),
               tolerance = 1e-6)
})
