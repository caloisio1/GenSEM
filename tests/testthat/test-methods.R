test_that("S3 interface works on a small recursive fit", {
  set.seed(99)
  n <- 200; Ks <- 2; Kh <- 2
  Sig <- matrix(c(0.3, 0.1, 0.1, 0.4), 2)
  uv <- matrix(rnorm(2*n), n, 2) %*% chol(Sig)
  id1 <- rep(1:n, each = Ks); zs <- rnorm(n*Ks); xs <- rnorm(n*Ks)
  y1 <- rbinom(n*Ks, 1, pnorm(-0.2 + 0.6*zs + 0.3*xs + uv[id1, 2]))
  dsel <- data.frame(id = id1, z = zs, x = xs, y1 = y1)
  id2 <- rep(1:n, each = Kh); xh <- rnorm(n*Kh)
  dur <- rexp(n*Kh, 1) + 0.5
  y1h <- y1[match(id2, id1) + sample(0:(Ks-1), n*Kh, TRUE)]
  ev <- rpois(n*Kh, dur * exp(-1 - 0.4*y1h + 0.3*xh + uv[id2, 1]))
  dhaz <- data.frame(id = id2, x = xh, y1 = y1h, dur, ev)
  fit <- fit_recursive(y1 ~ z + x, ev ~ y1 + x, dsel, dhaz,
                       id = "id", exposure = "dur", Q = 5)

  expect_identical(coef(fit), fit$coef)
  expect_identical(vcov(fit), fit$vcov)
  ci <- confint(fit)
  expect_equal(rownames(ci), names(coef(fit)))
  expect_true(all(ci[, 1] < coef(fit) & coef(fit) < ci[, 2]))
  ci9 <- confint(fit, parm = "haz:y1", level = 0.9)
  expect_true(ci[ "haz:y1", 1] < ci9[1, 1])   # 90% inside 95%
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -fit$nll)
  expect_equal(attr(ll, "df"), length(coef(fit)))
  expect_output(print(fit), "Recursive multiprocess system")
  s <- summary(fit)
  expect_s3_class(s, "summary.gensem_recursive")
  expect_output(print(s), "Pr\\(>\\|z\\|\\)")
  expect_equal(unname(s$coefficients[, "Estimate"]), unname(coef(fit)))
})
