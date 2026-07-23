# Seeded-truth validation for the gamma and loglogistic AFT hazard kernels.
# Simulations match the IMPLEMENTED parameterisations exactly:
#   gamma      : rate = exp(-lp)  =>  T = exp(lp + u) * Gamma(shape k, rate 1)
#   loglogistic: z = (log t - lp)/s  =>  log T = lp + u + s * rlogis()

sim_sel <- function(n, Ks, Sig, g0 = -0.3, g_z = 0.7, g_x = 0.4, seed) {
  set.seed(seed)
  uv <- matrix(rnorm(2 * n), n, 2) %*% chol(Sig)
  id1 <- rep(1:n, each = Ks)
  zs <- rnorm(n * Ks); xs <- rnorm(n * Ks)
  y1 <- rbinom(n * Ks, 1, pnorm(g0 + g_z * zs + g_x * xs + uv[id1, 2]))
  list(u = uv[, 1], dsel = data.frame(id = id1, z = zs, x = xs, y1 = y1),
       id1 = id1, y1 = y1, Ks = Ks)
}

sim_haz_frame <- function(n, Kh, s) {
  id2 <- rep(1:n, each = Kh)
  xh <- rnorm(n * Kh)
  y1h <- s$y1[match(id2, s$id1) + sample(0:(s$Ks - 1), n * Kh, TRUE)]
  list(id2 = id2, xh = xh, y1h = y1h)
}

test_that("gamma AFT kernel recovers seeded shape and alpha", {
  skip_on_cran()
  n <- 600; Kh <- 2; a2 <- -0.5; b2x <- 0.3; b20 <- -0.5; k <- 1.6
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  s <- sim_sel(n, 2, Sig, seed = 41)
  h <- sim_haz_frame(n, Kh, s)
  lp <- b20 + a2 * h$y1h + b2x * h$xh + s$u[h$id2]
  Tt <- exp(lp) * rgamma(n * Kh, shape = k, rate = 1)
  cens <- 2.5; d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = h$id2, x = h$xh, y1 = h$y1h, t = tt, d = d)
  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x, s$dsel, dhaz, id = "id",
                       family = "gamma", time = "t", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log shape"] - log(k)), 1.96 * se["log shape"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
})

test_that("loglogistic AFT kernel recovers seeded scale and alpha", {
  skip_on_cran()
  n <- 600; Kh <- 2; a2 <- -0.5; b2x <- 0.3; b20 <- -0.3; sc <- 0.6
  Sig <- matrix(c(0.4, 0.45*sqrt(0.4*0.8), 0.45*sqrt(0.4*0.8), 0.8), 2)
  s <- sim_sel(n, 2, Sig, seed = 43)
  h <- sim_haz_frame(n, Kh, s)
  lp <- b20 + a2 * h$y1h + b2x * h$xh + s$u[h$id2]
  Tt <- exp(lp + sc * rlogis(n * Kh))
  cens <- 3; d <- as.integer(Tt <= cens); tt <- pmin(Tt, cens)
  dhaz <- data.frame(id = h$id2, x = h$xh, y1 = h$y1h, t = tt, d = d)
  fit <- fit_recursive(y1 ~ z + x, d ~ y1 + x, s$dsel, dhaz, id = "id",
                       family = "loglogistic", time = "t", Q = 9)
  expect_equal(fit$convergence, 0)
  co <- fit$coef; se <- sqrt(diag(fit$vcov))
  expect_lt(abs(co["log scale"] - log(sc)), 1.96 * se["log scale"])
  expect_lt(abs(co["haz:y1"] - a2), 1.96 * se["haz:y1"])
  expect_lt(abs(co["sel:z"] - 0.7), 1.96 * se["sel:z"])
  # AFT truncation sanity: entry column of zeros must not change estimates
  dhaz$t0 <- 0
  fit0 <- fit_recursive(y1 ~ z + x, d ~ y1 + x, s$dsel, dhaz, id = "id",
                        family = "loglogistic", time = "t", entry = "t0", Q = 9)
  expect_equal(unname(fit0$coef["haz:y1"]), unname(co["haz:y1"]),
               tolerance = 1e-6)
})
