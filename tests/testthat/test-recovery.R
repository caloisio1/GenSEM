# Simulation-based validation: seed a nonrecursive two-process system with
# known structural parameters, simulate its reduced form (Bartus 2017, eq. 3)
# on MULTISPELL data (several episodes per person -- the data structure this
# model class is built for, and the regime where the Laplace approximation is
# adequate), fit the stacked model, and recover lambda/beta via eqs. 5-6.

expect_close <- function(actual, expected, tol) {
  expect_lt(max(abs(actual - expected)), tol)
}

test_that("stacked fit + delta method recover seeded structural parameters", {
  skip_if_not_installed("glmmTMB")
  set.seed(20260721)
  n <- 2000; K <- 4                      # persons x spells per process
  # --- seeded structural truth ---------------------------------------------
  lam1 <- 0.5; lam2 <- -0.4              # cross-hazard selection effects
  b1 <- 0.3;  b2 <- -0.2                 # structural effects of common x
  g1 <- 0.6;  g2 <- 0.7                  # excluded instruments z1, z2
  a1 <- -1.0; a2 <- -1.2
  Sig <- matrix(c(0.30, 0.10, 0.10, 0.25), 2)
  d <- 1 / (1 - lam1 * lam2)             # reduced form, Bartus eq. (3)
  p10 <- d * (b1 + lam1 * b2); p11 <- d * g1; p12 <- d * lam1 * g2
  p20 <- d * (b2 + lam2 * b1); p22 <- d * g2; p21 <- d * lam2 * g1
  # --- simulate multispell data --------------------------------------------
  L <- chol(Sig); u <- matrix(rnorm(2 * n), n, 2) %*% L   # cov = t(L)%*%L = Sig
  id <- rep(seq_len(n), each = K)
  x  <- rnorm(n * K); z1 <- rnorm(n * K); z2 <- rnorm(n * K)
  dur <- rexp(n * K, 1) + 0.5
  y1 <- rpois(n * K, dur * exp(a1 + p10*x + p11*z1 + p12*z2 + u[id, 1]))
  y2 <- rpois(n * K, dur * exp(a2 + p20*x + p21*z1 + p22*z2 + u[id, 2]))
  dat <- data.frame(id, x, z1, z2, dur, ev1 = y1, ev2 = y2)

  long <- stack_processes(dat, id = "id",
                          events = c(pA = "ev1", pB = "ev2"),
                          exposure = "dur")
  expect_equal(nrow(long), 2L * n * K)
  expect_setequal(levels(long$process), c("pA", "pB"))

  fit <- fit_multiproc(long, fixed = ~ x + z1 + z2, id = "id")
  bhat <- glmmTMB::fixef(fit)$cond

  # reduced-form recovery (absolute tolerances)
  expect_close(bhat["processpA:x"],  p10, 0.04)
  expect_close(bhat["processpB:x"],  p20, 0.04)
  expect_close(bhat["processpA:z2"], p12, 0.04)
  expect_close(bhat["processpB:z2"], p22, 0.04)

  # random-effect covariance recovery
  Shat <- glmmTMB::VarCorr(fit)$cond$id
  expect_close(diag(Shat), diag(Sig), 0.05)
  expect_close(Shat[1, 2], Sig[1, 2], 0.05)

  # structural recovery, Bartus eqs. (5)-(6)
  lam1_hat <- recover_structural(fit, "`processpA:z2` / `processpB:z2`")
  b1_hat <- recover_structural(fit,
    "`processpA:x` - (`processpA:z2`/`processpB:z2`) * `processpB:x`")
  expect_close(lam1_hat$estimate, lam1, 0.08)
  expect_close(b1_hat$estimate,   b1,   0.08)
  expect_lt(lam1_hat$ci_lo, lam1); expect_gt(lam1_hat$ci_hi, lam1)
  expect_lt(b1_hat$ci_lo,   b1);   expect_gt(b1_hat$ci_hi,   b1)
})

test_that("single-spell sparse data documents the Laplace variance bias", {
  # With ONE spell per person the Laplace approximation overestimates the
  # random-effect variances (adaptive quadrature territory). This test pins
  # the known limitation so it is documented, not silent.
  skip_if_not_installed("glmmTMB")
  set.seed(2)
  n <- 4000
  u <- rnorm(n, 0, sqrt(0.3)); x <- rnorm(n); dur <- rexp(n, 1) + 0.5
  dat <- data.frame(id = seq_len(n), x, dur,
                    ev1 = rpois(n, dur * exp(-1 + 0.3 * x + u)))
  long <- stack_processes(dat, "id", c(pA = "ev1"), "dur")
  fit <- fit_multiproc(long, ~ x, id = "id", shared_re = TRUE)
  v <- attr(glmmTMB::VarCorr(fit)$cond$id, "stddev")^2
  expect_gt(unname(v), 0.3)   # known upward bias in the sparse regime
})

test_that("shared_re handles the recurrent-events case", {
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  n <- 800; K <- 4
  id <- rep(seq_len(n), each = K)
  u <- rnorm(n, 0, 0.5); x <- rnorm(n * K); dur <- rexp(n * K, 1) + 0.5
  dat <- data.frame(id, x, dur,
                    ev1 = rpois(n * K, dur * exp(-1 + 0.3 * x + u[id])),
                    ev2 = rpois(n * K, dur * exp(-1.2 + 0.3 * x + u[id])))
  long <- stack_processes(dat, "id", c(k1 = "ev1", k2 = "ev2"), "dur")
  fit <- fit_multiproc(long, ~ x, id = "id", shared_re = TRUE)
  sd_hat <- attr(glmmTMB::VarCorr(fit)$cond$id, "stddev")
  expect_close(sd_hat, 0.5, 0.06)
})

test_that("GLMMadaptive backend corrects the sparse-regime Laplace bias", {
  skip_if_not_installed("GLMMadaptive")
  set.seed(103)  # a seed from the Monte Carlo band
  n <- 2000
  u <- rnorm(n, 0, sqrt(0.3)); x <- rnorm(n); dur <- rexp(n, 1) + 0.5
  dat <- data.frame(id = seq_len(n), x, dur,
                    ev = rpois(n, dur * exp(-1 + 0.3 * x + u)))
  long <- stack_processes(dat, "id", c(pA = "ev"), "dur")
  f_lap <- fit_multiproc(long, ~ x, id = "id", shared_re = TRUE)
  f_agq <- fit_multiproc(long, ~ x, id = "id", shared_re = TRUE,
                         backend = "GLMMadaptive", nAGQ = 11)
  v_lap <- unname(attr(glmmTMB::VarCorr(f_lap)$cond$id, "stddev")^2)
  v_agq <- unname(f_agq$D[1, 1])
  expect_lt(abs(v_agq - 0.3), abs(v_lap - 0.3))  # AGHQ closer to truth
  # recover_structural must also work on MixMod fits
  s <- recover_structural(f_agq, "`x` * 2")
  expect_equal(s$estimate, 2 * unname(GLMMadaptive::fixef(f_agq)["x"]),
               tolerance = 1e-8)
})
