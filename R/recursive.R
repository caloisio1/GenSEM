# Recursive mixed-family systems: probit or mlogit selection + hazard equation
# (Poisson/piecewise-exponential OR Weibull PH with left truncation) with
# correlated bivariate normal random intercepts (Bartus 2017, secs. 2.3, 5;
# Lillard, Brien & Waite 1995). The joint marginal likelihood is written out
# explicitly and integrated by 2D Gauss-Hermite quadrature (Golub-Welsch
# nodes) -- every step inspectable, no compiled black box.
#
# Hazard kernels (per spell, linear predictor lp = x'b + u):
#   poisson : ev*(lp + log dur) - exp(lp + log dur)          [PWE / counts]
#   weibull : d*(log p + (p-1) log t + lp) - exp(lp)*(t^p - t0^p)
#             (proportional hazards, shape p, left truncation at t0)

gh_nodes <- function(Q) {
  i <- seq_len(Q - 1)
  J <- matrix(0, Q, Q)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

haz_kernel <- function(H, lp2, anc) {
  if (H$fam == "poisson") {
    H$ev * (lp2 + H$ldur) - exp(lp2 + H$ldur)
  } else if (H$fam == "weibull") {          # PH; anc = shape p
    H$d * (log(anc) + (anc - 1) * H$lt + lp2) -
      exp(lp2) * (H$t^anc - H$t0^anc)
  } else if (H$fam == "gamma") {            # AFT; anc = shape k; rate = e^{-lp}
    lp2 <- pmin(pmax(lp2, -30), 30)         # numerical guard: at |lp|>30 the
    rate <- exp(-lp2)                       # contribution is already 0/1
    out <- H$d * stats::dgamma(H$t, shape = anc, rate = rate, log = TRUE) +
      (1 - H$d) * stats::pgamma(H$t, shape = anc, rate = rate,
                                lower.tail = FALSE, log.p = TRUE)
    tr <- H$t0 > 0
    if (any(tr)) out[tr] <- out[tr] -
      stats::pgamma(H$t0[tr], shape = anc, rate = rate[tr],
                    lower.tail = FALSE, log.p = TRUE)
    out
  } else if (H$fam == "loglogistic") {      # AFT; anc = scale s
    lp2 <- pmin(pmax(lp2, -30), 30)
    z <- (H$lt - lp2) / anc
    out <- H$d * (stats::dlogis(z, log = TRUE) - log(anc) - H$lt) +
      (1 - H$d) * stats::plogis(z, lower.tail = FALSE, log.p = TRUE)
    tr <- H$t0 > 0
    if (any(tr)) {
      z0 <- (log(H$t0[tr]) - lp2[tr]) / anc
      out[tr] <- out[tr] - stats::plogis(z0, lower.tail = FALSE, log.p = TRUE)
    }
    out
  } else {                                  # lognormal AFT; anc = sigma
    lp2 <- pmin(pmax(lp2, -30), 30)
    out <- H$d * (stats::dnorm(H$lt, mean = lp2, sd = anc, log = TRUE) - H$lt) +
      (1 - H$d) * stats::pnorm(H$lt, mean = lp2, sd = anc,
                               lower.tail = FALSE, log.p = TRUE)
    tr <- H$t0 > 0
    if (any(tr)) out[tr] <- out[tr] -
      stats::pnorm(log(H$t0[tr]), mean = lp2[tr], sd = anc,
                   lower.tail = FALSE, log.p = TRUE)
    out
  }
}

nll_recursive <- function(theta, X1, y1, id1, X2, H, id2, ids, gh) {
  p1 <- ncol(X1); p2 <- ncol(X2); wb <- H$fam != "poisson"
  b1 <- theta[1:p1]; b2 <- theta[(p1 + 1):(p1 + p2)]
  shape <- if (wb) exp(theta[p1 + p2 + 1]) else 1
  j <- p1 + p2 + as.integer(wb)
  s1 <- exp(theta[j + 1]); s2 <- exp(theta[j + 2]); rho <- tanh(theta[j + 3])
  L <- matrix(c(s1, s2 * rho, 0, s2 * sqrt(1 - rho^2)), 2, 2)
  eta1 <- drop(X1 %*% b1); eta2 <- drop(X2 %*% b2)
  n <- length(ids); Q <- length(gh$z)
  ll <- matrix(-Inf, n, Q * Q); m <- 0L
  for (k in seq_len(Q)) for (l in seq_len(Q)) {
    m <- m + 1L
    uv <- sqrt(2) * (L %*% c(gh$z[k], gh$z[l]))
    lp1 <- eta1 + uv[2]
    c_sel <- y1 * stats::pnorm(lp1, log.p = TRUE) +
             (1 - y1) * stats::pnorm(lp1, lower.tail = FALSE, log.p = TRUE)
    c_haz <- haz_kernel(H, eta2 + uv[1], shape)
    ll[, m] <- rowsum(c_sel, id1)[ids, 1] + rowsum(c_haz, id2)[ids, 1] +
               log(gh$w[k]) + log(gh$w[l]) - log(pi)
  }
  mx <- apply(ll, 1, max)
  v <- -sum(mx + log(rowSums(exp(ll - mx))))
  if (!is.finite(v)) 1e10 else v
}

nll_recursive_mlogit <- function(theta, X1, y1i, C, id1, X2, H, id2, ids, gh) {
  p1 <- ncol(X1); p2 <- ncol(X2); nc <- C - 1L; wb <- H$fam != "poisson"
  B1 <- matrix(theta[1:(nc * p1)], p1, nc)
  b2 <- theta[(nc * p1 + 1):(nc * p1 + p2)]
  lv <- c(1, if (nc > 1) theta[(nc * p1 + p2 + 1):(nc * p1 + p2 + nc - 1)])
  j <- nc * p1 + p2 + nc - 1L
  shape <- if (wb) exp(theta[j + 1]) else 1
  j <- j + as.integer(wb)
  s1 <- exp(theta[j + 1]); s2 <- exp(theta[j + 2]); rho <- tanh(theta[j + 3])
  L <- matrix(c(s1, s2 * rho, 0, s2 * sqrt(1 - rho^2)), 2, 2)
  E1 <- X1 %*% B1; eta2 <- drop(X2 %*% b2)
  n <- length(ids); Q <- length(gh$z)
  pick <- cbind(seq_along(y1i), pmax(y1i, 1L))
  ll <- matrix(-Inf, n, Q * Q); m <- 0L
  for (k in seq_len(Q)) for (l in seq_len(Q)) {
    m <- m + 1L
    uv <- sqrt(2) * (L %*% c(gh$z[k], gh$z[l]))
    A <- sweep(E1, 2, lv * uv[2], "+")
    m0 <- pmax(0, apply(A, 1, max))
    lden <- m0 + log(exp(-m0) + rowSums(exp(A - m0)))
    c_sel <- ifelse(y1i > 0L, A[pick], 0) - lden
    c_haz <- haz_kernel(H, eta2 + uv[1], shape)
    ll[, m] <- rowsum(c_sel, id1)[ids, 1] + rowsum(c_haz, id2)[ids, 1] +
               log(gh$w[k]) + log(gh$w[l]) - log(pi)
  }
  mx <- apply(ll, 1, max)
  v <- -sum(mx + log(rowSums(exp(ll - mx))))
  if (!is.finite(v)) 1e10 else v
}

.build_haz <- function(haz, data_haz, family, exposure, time, entry) {
  X2 <- stats::model.matrix(haz, data_haz)
  yv <- stats::model.response(stats::model.frame(haz, data_haz))
  if (family == "poisson") {
    if (is.null(exposure) || !exposure %in% names(data_haz))
      stop("family='poisson' requires `exposure`.")
    if (any(data_haz[[exposure]] <= 0, na.rm = TRUE))
      stop("All exposures must be strictly positive.")
    H <- list(fam = "poisson", ev = yv, ldur = log(data_haz[[exposure]]))
  } else {
    if (is.null(time) || !time %in% names(data_haz))
      stop(sprintf("family='%s' requires `time` (exit time).", family))
    if (!all(yv %in% c(0, 1)))
      stop("Survival-family response must be a 0/1 event.")
    t <- data_haz[[time]]
    t0 <- if (!is.null(entry)) data_haz[[entry]] else rep(0, length(t))
    if (any(t <= t0)) stop("time must exceed entry for every spell.")
    H <- list(fam = family, d = yv, t = t, t0 = t0, lt = log(t))
  }
  list(X2 = X2, H = H, yv = yv)
}

.start_haz <- function(X2, H) {
  if (H$fam == "poisson")
    return(stats::coef(stats::glm.fit(X2, H$ev, offset = H$ldur,
                                      family = stats::poisson())))
  b <- stats::coef(stats::glm.fit(X2, H$d, offset = log(H$t - H$t0),
                                  family = stats::poisson()))
  # exponential PH approximation; AFT families (gamma, loglogistic) have
  # lp entering as log location, so the exponential AFT start is -b_PH
  if (H$fam == "weibull") b else -b
}

#' Fit a recursive mixed-family system (probit selection + hazard)
#'
#' Jointly estimates a probit equation for an endogenous binary exposure and a
#' hazard equation that includes the observed exposure, with correlated
#' person-level random intercepts -- Bartus's (2017) recursive class, eq. (7).
#' All coefficients are structural. Identification of the random-effect
#' correlation requires at least one selection-only variable (an excluded
#' instrument) in `sel`.
#'
#' Four hazard families. `"poisson"` (piecewise-constant exponential via the
#' Poisson trick; response = event count, plus `exposure`) and `"weibull"`
#' (proportional hazards, estimated shape) are the two specifications Bartus
#' fits in his Example 2. `"gamma"` and `"loglogistic"` are ACCELERATED
#' FAILURE TIME families (they have no proportional-hazards form):
#' coefficients are log TIME ratios -- positive = longer survival = LOWER
#' hazard, the OPPOSITE sign convention of poisson/weibull. The loglogistic
#' allows a non-monotone hazard. Survival families take a 0/1 event response
#' plus `time` and optional `entry` (left truncation).
#'
#' The marginal likelihood is integrated by non-adaptive 2D Gauss-Hermite
#' quadrature with `Q^2` nodes and maximised by BFGS with numerical
#' derivatives; standard errors come from the numerical Hessian. Transparent
#' but not fast: minutes, not seconds, at Q = 15. Start values come from
#' separate single-equation GLMs (exponential approximation for Weibull).
#'
#' @param sel probit formula, e.g. `y1 ~ z + x` (include the instrument z).
#' @param haz hazard formula, e.g. `event ~ y1 + x` (include the endogenous
#'   dummy).
#' @param data_sel data for the selection equation.
#' @param data_haz spell-level data for the hazard equation.
#' @param id person identifier column, present in both data sets.
#' @param family `"poisson"` (default), `"weibull"` (PH), `"gamma"` or
#'   `"loglogistic"` (both AFT -- opposite coefficient sign convention).
#' @param exposure spell-duration column (> 0); poisson family only.
#' @param time exit-time column (> entry); survival families only.
#' @param entry optional entry-time column for left truncation (weibull;
#'   default 0).
#' @param Q quadrature points per dimension (default 15; `Q^2` total).
#' @param nested optional column in `data_haz` identifying a level-2 cluster
#'   within each `id` unit (children within mothers, youths within schools,
#'   jobs within workers). When given, a cluster-level frailty w ~ N(0, sw^2)
#'   is added to the HAZARD linear predictor and integrated out by nested
#'   Gauss-Hermite quadrature; selection and (u, v) stay at the `id` level.
#' @param Qw quadrature points for the nested frailty integral (default 7).
#' @param start optional full start vector for the optimizer (replaces the
#'   GLM-based defaults); length must match the parameter vector.
#' @return list of class `gensem_recursive` with `coef` (named:
#'   `sel:`/`haz:` prefixes, `log shape`/`log scale` for survival families,
#'   `log s1`, `log s2`, `atanh rho`), `vcov`, `Sigma`, `family`,
#'   `ancillary` (shape or scale on the natural scale), `convergence`,
#'   `nll`, `Q`. Compatible with [recover_structural()].
#' @export
fit_recursive <- function(sel, haz, data_sel, data_haz, id,
                          family = c("poisson", "weibull", "gamma",
                                     "loglogistic", "lognormal"),
                          exposure = NULL, time = NULL, entry = NULL, Q = 15,
                          nested = NULL, Qw = 7, start = NULL) {
  family <- match.arg(family)
  stopifnot(inherits(sel, "formula"), inherits(haz, "formula"),
            id %in% names(data_sel), id %in% names(data_haz))
  X1 <- stats::model.matrix(sel, data_sel)
  y1 <- stats::model.response(stats::model.frame(sel, data_sel))
  bh <- .build_haz(haz, data_haz, family, exposure, time, entry)
  ids <- sort(unique(c(data_sel[[id]], data_haz[[id]])))
  id1 <- factor(data_sel[[id]], levels = ids)
  id2 <- factor(data_haz[[id]], levels = ids)
  b1_0 <- stats::coef(stats::glm.fit(X1, y1,
            family = stats::binomial("probit")))
  b2_0 <- .start_haz(bh$X2, bh$H)
  wb <- family != "poisson"
  anc_nm <- c(weibull = "log shape", gamma = "log shape",
              loglogistic = "log scale", lognormal = "log sigma")[family]
  gh <- gh_nodes(Q)
  if (is.null(nested)) {
    theta0 <- c(b1_0, b2_0, if (wb) 0, log(0.5), log(0.5), 0)
    if (!is.null(start)) { stopifnot(length(start) == length(theta0)); theta0 <- start }
    opt <- stats::optim(theta0, nll_recursive, method = "BFGS", hessian = TRUE,
                        X1 = X1, y1 = y1, id1 = id1, X2 = bh$X2, H = bh$H,
                        id2 = id2, ids = seq_along(ids), gh = gh,
                        control = list(maxit = 600))
    nm <- c(paste0("sel:", colnames(X1)), paste0("haz:", colnames(bh$X2)),
            if (wb) unname(anc_nm), "log s1", "log s2", "atanh rho")
  } else {
    if (!nested %in% names(data_haz)) stop("`nested` column not in data_haz.")
    id_child <- factor(interaction(data_haz[[id]], data_haz[[nested]],
                                   drop = TRUE))
    ch2m <- factor(tapply(as.character(id2), id_child, `[`, 1),
                   levels = levels(id2))
    theta0 <- c(b1_0, b2_0, if (wb) 0, log(0.4), log(0.5), log(0.5), 0)
    if (!is.null(start)) { stopifnot(length(start) == length(theta0)); theta0 <- start }
    opt <- stats::optim(theta0, nll_recursive_3l, method = "BFGS",
                        hessian = TRUE, X1 = X1, y1 = y1, id1 = id1,
                        X2 = bh$X2, H = bh$H, id_child = id_child,
                        ch2m = ch2m, ids = seq_along(ids), gh = gh,
                        ghw = gh_nodes(Qw), control = list(maxit = 700))
    nm <- c(paste0("sel:", colnames(X1)), paste0("haz:", colnames(bh$X2)),
            if (wb) unname(anc_nm), "log s_w(nested)", "log s1", "log s2",
            "atanh rho")
  }
  co <- stats::setNames(opt$par, nm)
  V <- solve(opt$hessian); dimnames(V) <- list(nm, nm)
  s1 <- exp(co["log s1"]); s2 <- exp(co["log s2"]); r <- tanh(co["atanh rho"])
  Sigma <- matrix(c(s1^2, r * s1 * s2, r * s1 * s2, s2^2), 2,
                  dimnames = list(c("u(haz)", "v(sel)"), c("u(haz)", "v(sel)")))
  structure(list(coef = co, vcov = V, Sigma = Sigma, family = family,
                 ancillary = if (wb) unname(exp(co[anc_nm])),
                 shape = if (wb) unname(exp(co[anc_nm])),  # alias
                 sw2 = if (!is.null(nested))
                         unname(exp(co["log s_w(nested)"])^2),
                 convergence = opt$convergence, nll = opt$value, Q = Q),
            class = "gensem_recursive")
}

#' Fit a recursive system with multinomial (mlogit) selection
#'
#' As [fit_recursive()], but the endogenous exposure takes C >= 2 unordered
#' categories, modelled by multinomial logit (first level = base outcome).
#' The selection random effect v enters category 1 with loading fixed to 1
#' and categories 2..C-1 with estimated loadings, matching gsem's
#' parameterisation in Bartus (2017, sec. 5.3). Include the exposure in
#' `haz` as a factor. Both hazard families of [fit_recursive()] are
#' available. Identification requires a selection-only instrument.
#'
#' @inheritParams fit_recursive
#' @param start optional full start vector for the optimizer (replaces the
#'   GLM-based defaults); length must match the parameter vector.
#' @param sel mlogit formula; the response must be a factor (or integer
#'   codes 0..C-1), first level = base outcome.
#' @return list of class `gensem_recursive` (see [fit_recursive()]);
#'   coefficient names carry `sel<cat>:` prefixes and `load v:<cat>` for the
#'   estimated loadings.
#' @export
fit_recursive_mlogit <- function(sel, haz, data_sel, data_haz, id,
                                 family = c("poisson", "weibull", "gamma",
                                            "loglogistic", "lognormal"),
                                 exposure = NULL, time = NULL, entry = NULL,
                                 Q = 15, start = NULL) {
  family <- match.arg(family)
  stopifnot(inherits(sel, "formula"), inherits(haz, "formula"),
            id %in% names(data_sel), id %in% names(data_haz))
  X1 <- stats::model.matrix(sel, data_sel)
  yr <- stats::model.response(stats::model.frame(sel, data_sel))
  yf <- if (is.factor(yr)) yr else factor(yr)
  C <- nlevels(yf); if (C < 2) stop("selection response needs >= 2 levels")
  y1i <- as.integer(yf) - 1L
  bh <- .build_haz(haz, data_haz, family, exposure, time, entry)
  ids <- sort(unique(c(data_sel[[id]], data_haz[[id]])))
  id1 <- factor(data_sel[[id]], levels = ids)
  id2 <- factor(data_haz[[id]], levels = ids)
  nc <- C - 1L
  b1_0 <- unlist(lapply(seq_len(nc), function(c) {
    stats::coef(stats::glm.fit(X1[y1i %in% c(0L, c), , drop = FALSE],
      as.integer(y1i[y1i %in% c(0L, c)] == c),
      family = stats::binomial()))
  }))
  b2_0 <- .start_haz(bh$X2, bh$H)
  wb <- family != "poisson"
  anc_nm <- c(weibull = "log shape", gamma = "log shape",
              loglogistic = "log scale", lognormal = "log sigma")[family]
  theta0 <- c(b1_0, b2_0, rep(0.5, nc - 1), if (wb) 0,
              log(0.5), log(0.5), 0)
  if (!is.null(start)) { stopifnot(length(start) == length(theta0)); theta0 <- start }
  gh <- gh_nodes(Q)
  opt <- stats::optim(theta0, nll_recursive_mlogit, method = "BFGS",
                      hessian = TRUE, X1 = X1, y1i = y1i, C = C, id1 = id1,
                      X2 = bh$X2, H = bh$H, id2 = id2,
                      ids = seq_along(ids), gh = gh,
                      control = list(maxit = 700))
  cats <- levels(yf)[-1]
  nm <- c(unlist(lapply(cats, function(cc) paste0("sel", cc, ":", colnames(X1)))),
          paste0("haz:", colnames(bh$X2)),
          if (nc > 1) paste0("load v:", cats[-1]),
          if (wb) unname(anc_nm), "log s1", "log s2", "atanh rho")
  co <- stats::setNames(opt$par, nm)
  V <- solve(opt$hessian); dimnames(V) <- list(nm, nm)
  s1 <- exp(co["log s1"]); s2 <- exp(co["log s2"]); r <- tanh(co["atanh rho"])
  Sigma <- matrix(c(s1^2, r * s1 * s2, r * s1 * s2, s2^2), 2,
                  dimnames = list(c("u(haz)", "v(sel)"), c("u(haz)", "v(sel)")))
  structure(list(coef = co, vcov = V, Sigma = Sigma, family = family,
                 ancillary = if (wb) unname(exp(co[anc_nm])),
                 shape = if (wb) unname(exp(co[anc_nm])),  # alias
                 convergence = opt$convergence, nll = opt$value, Q = Q),
            class = "gensem_recursive")
}

#' Split spells at cutpoints for piecewise-constant duration dependence
#'
#' Episode splitting: each spell is divided at `breaks`, producing one row per
#' occupied segment with adjusted exposure, the (0/1) event assigned to the
#' final occupied segment, and an `interval` factor to include among the
#' hazard covariates. With enough cutpoints, the piecewise-constant
#' exponential approximates arbitrary duration dependence (Bartus 2017,
#' p. 459); for monotone hazards, `family = "weibull"` in the recursive
#' fitters is the parametric alternative.
#'
#' @param data spell-level data.frame with a 0/1 event indicator.
#' @param event name of the 0/1 event column.
#' @param exposure name of the spell-duration column (> 0).
#' @param breaks increasing positive cutpoints, e.g. `c(3, 12)`.
#' @return data.frame with one row per occupied segment: original columns,
#'   exposure and event replaced segment-wise, plus `interval` (factor)
#'   and `.t0` (segment entry time).
#' @export
split_episodes <- function(data, event, exposure, breaks) {
  stopifnot(event %in% names(data), exposure %in% names(data),
            all(breaks > 0), !is.unsorted(breaks, strictly = TRUE),
            all(data[[event]] %in% c(0, 1)))
  lo <- c(0, breaks); hi <- c(breaks, Inf)
  segs <- lapply(seq_along(lo), function(s) {
    dur <- pmin(data[[exposure]], hi[s]) - lo[s]
    keep <- dur > 0
    out <- data[keep, , drop = FALSE]
    out$.t0 <- lo[s]
    out[[exposure]] <- dur[keep]
    last <- data[[exposure]][keep] <= hi[s]
    out[[event]] <- ifelse(last, out[[event]], 0)
    out$interval <- s
    out
  })
  out <- do.call(rbind, segs)
  out$interval <- factor(out$interval, levels = seq_along(lo),
                         labels = paste0("[", lo, ",",
                                         ifelse(is.finite(hi), hi, "Inf"), ")"))
  rownames(out) <- NULL
  out
}

# --- Three-level extension: spells in CHILDREN in MOTHERS -------------------
# Hazard lp = x'b + u_mother + w_child, w ~ N(0, sw^2) independent of (u, v);
# selection stays mother-level. Key structural fact exploited below: with the
# lower Cholesky of Sigma, u depends only on the FIRST quadrature axis, so the
# inner 1D integral over w is evaluated per u-node (Q x Qw kernel passes),
# not per (u, v) pair (Q^2 x Qw).

nll_recursive_3l <- function(theta, X1, y1, id1, X2, H, id_child, ch2m,
                             ids, gh, ghw) {
  p1 <- ncol(X1); p2 <- ncol(X2); wb <- H$fam != "poisson"
  b1 <- theta[1:p1]; b2 <- theta[(p1 + 1):(p1 + p2)]
  anc <- if (wb) exp(theta[p1 + p2 + 1]) else 1
  j <- p1 + p2 + as.integer(wb)
  sw <- exp(theta[j + 1])
  s1 <- exp(theta[j + 2]); s2 <- exp(theta[j + 3]); rho <- tanh(theta[j + 4])
  L <- matrix(c(s1, s2 * rho, 0, s2 * sqrt(1 - rho^2)), 2, 2)
  eta1 <- drop(X1 %*% b1); eta2 <- drop(X2 %*% b2)
  n <- length(ids); Q <- length(gh$z); Qw <- length(ghw$z)
  nch <- nlevels(id_child)
  # hazard block: per u-node, integrate w out child by child
  hazmat <- matrix(0, n, Q)
  for (k in seq_len(Q)) {
    u_k <- sqrt(2) * L[1, 1] * gh$z[k]
    Scj <- matrix(0, nch, Qw)
    for (jw in seq_len(Qw)) {
      w_j <- sqrt(2) * sw * ghw$z[jw]
      Scj[, jw] <- rowsum(haz_kernel(H, eta2 + u_k + w_j, anc),
                          id_child)[, 1] + log(ghw$w[jw]) - 0.5 * log(pi)
    }
    mx <- apply(Scj, 1, max)
    logLc <- mx + log(rowSums(exp(Scj - mx)))            # per-child integral
    hazmat[, k] <- rowsum(logLc, ch2m)[ids, 1]           # sum to mothers
  }
  ll <- matrix(-Inf, n, Q * Q); m <- 0L
  for (k in seq_len(Q)) for (l in seq_len(Q)) {
    m <- m + 1L
    v_kl <- sqrt(2) * (L[2, 1] * gh$z[k] + L[2, 2] * gh$z[l])
    lp1 <- eta1 + v_kl
    c_sel <- y1 * stats::pnorm(lp1, log.p = TRUE) +
             (1 - y1) * stats::pnorm(lp1, lower.tail = FALSE, log.p = TRUE)
    ll[, m] <- hazmat[, k] + rowsum(c_sel, id1)[ids, 1] +
               log(gh$w[k]) + log(gh$w[l]) - log(pi)
  }
  mx <- apply(ll, 1, max)
  -sum(mx + log(rowSums(exp(ll - mx))))
}
