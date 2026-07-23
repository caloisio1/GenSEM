sim_sel <- function(n, Ks, Sig) {
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig)
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(-0.3 + 0.7*zs + 0.4*xs + uv[id1, 2]))
  list(u = uv[,1], id1 = id1,
       dsel = data.frame(id = id1, z = zs, x = xs, y1 = y1), y1 = y1, Ks = Ks)
}

test_that("gamma AFT family recovers seeded shape and alpha", {
  skip_on_cran()
  set.seed(41)
  n <- 600; Kh <- 2; Sig <- matrix(c(0.4, 0.45*sqrt(0.32), 0.45*sqrt(0.32), 0.8), 2)
  S <- sim_sel(n, 2, Sig)
  a2 <- -0.5; b2x <- 0.3; b20 <- 0.5; k <- 1.5   # AFT: lp = log location
  id2 <- rep(1:n, each = Kh); xh <- rnorm(n*Kh)
  y1h <- S$y1[match(id2, S$id1) + sample(0:(S$Ks-1), n*Kh, TRUE)]
  lp <- b20 + a2*y1h + b2x*xh + S$u[id2]
  Tt <- rgamma(n*Kh, shape = k, rate = exp(-lp))   # AFT: T = e^lp * Gamma(k,1)
  cens <- quantile(Tt, .8)
  d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, t = tt, d = d)
  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x, S$dsel, dhaz, id = "id",
                       family = "gamma", time = "t", Q = 9)
  expect_equal(fit$convergence, 0); expect_equal(fit$family, "gamma")
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log shape"] - log(k)), 1.96 * se["log shape"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
})

test_that("loglogistic AFT family recovers seeded scale and alpha", {
  skip_on_cran()
  set.seed(43)
  n <- 600; Kh <- 2; Sig <- matrix(c(0.4, 0.45*sqrt(0.32), 0.45*sqrt(0.32), 0.8), 2)
  S <- sim_sel(n, 2, Sig)
  a2 <- -0.5; b2x <- 0.3; b20 <- 0.5; s <- 0.6
  id2 <- rep(1:n, each = Kh); xh <- rnorm(n*Kh)
  y1h <- S$y1[match(id2, S$id1) + sample(0:(S$Ks-1), n*Kh, TRUE)]
  lp <- b20 + a2*y1h + b2x*xh + S$u[id2]
  Tt <- exp(lp + s * rlogis(n*Kh))                 # log T = lp + s*logistic
  cens <- quantile(Tt, .8)
  d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, t = tt, d = d)
  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x, S$dsel, dhaz, id = "id",
                       family = "loglogistic", time = "t", Q = 9)
  expect_equal(fit$convergence, 0); expect_equal(fit$family, "loglogistic")
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log scale"] - log(s)), 1.96 * se["log scale"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
  expect_equal(fit$ancillary, unname(exp(co["log scale"])))
})
