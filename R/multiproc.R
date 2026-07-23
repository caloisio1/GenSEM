# multiproc: multilevel multiprocess hazard models via stacked Poisson estimation
# Core references: Lillard (1993, J. Econometrics); Bartus (2017, Stata Journal 17(2));
# Poisson trick: Holford (1980); Skrondal & Rabe-Hesketh (2004).

#' Stack several hazard processes into long format
#'
#' Each process is a piecewise-constant exponential hazard observed on spell
#' rows. Rows where the process's event indicator is missing are dropped for
#' that process. The result duplicates spell rows once per applicable process
#' and adds a `process` factor, an `.event` outcome, and a `.exposure` column.
#'
#' @param data spell-level data.frame (one row per spell).
#' @param id name of the person identifier column.
#' @param events named character vector: process name -> event indicator column
#'   (0/1, NA = spell not at risk for that process).
#' @param exposure name of the spell-duration column (>0).
#' @return long data.frame with columns `process`, `.event`, `.exposure`,
#'   plus `id` and all covariates.
#' @export
stack_processes <- function(data, id, events, exposure) {
  stopifnot(is.data.frame(data), id %in% names(data),
            exposure %in% names(data), all(events %in% names(data)),
            !is.null(names(events)), all(nzchar(names(events))))
  if (any(data[[exposure]] <= 0, na.rm = TRUE))
    stop("All exposures must be strictly positive.")
  pieces <- lapply(names(events), function(p) {
    ev <- data[[events[[p]]]]
    keep <- !is.na(ev)
    out <- data[keep, setdiff(names(data), unname(events)), drop = FALSE]
    out$process <- p
    out$.event <- as.integer(ev[keep])
    out$.exposure <- out[[exposure]]
    out
  })
  long <- do.call(rbind, pieces)
  long$process <- factor(long$process, levels = names(events))
  rownames(long) <- NULL
  long
}

#' Fit a multiprocess system as a stacked Poisson mixed model
#'
#' Fits `.event ~ 0 + process + process:(covariates) + (0 + process | id)`
#' with `offset(log(.exposure))` and Poisson family in glmmTMB. The
#' `(0 + process | id)` term is the unstructured covariance of the
#' process-specific random intercepts -- the object gsem writes as
#' `cov(U[id], V[id])`. One latent per process: processes sharing an
#' unobserved propensity (e.g., recurrent events of the same kind) should be
#' stacked under ONE process label with rank covariates, per Lillard (1993).
#'
#' Estimation is by Laplace approximation (glmmTMB), not adaptive quadrature
#' (gsem default): small numerical differences against gsem are expected and
#' should be reported, not hidden.
#'
#' @param long output of [stack_processes()].
#' @param fixed one-sided formula of covariates, e.g. `~ x1 + x2`. Every term
#'   is interacted with `process` (process-specific coefficients), matching
#'   the reduced-form specification in Bartus (2017, sec. 3.3).
#' @param id name of the person identifier column.
#' @param shared_re if TRUE, uses a single random intercept `(1 | id)` shared
#'   by all processes (Bartus's recurrent-events case) instead of
#'   process-specific correlated intercepts.
#' @param backend `"glmmTMB"` (Laplace approximation; fast, default) or
#'   `"GLMMadaptive"` (adaptive Gauss-Hermite quadrature, Rabe-Hesketh,
#'   Skrondal & Pickles 2005). Prefer `"GLMMadaptive"` in sparse regimes
#'   (few spells per person), where Laplace overestimates random-effect
#'   variances; quadrature cost grows as `nAGQ^k` with k processes.
#' @param nAGQ number of quadrature points per dimension
#'   (GLMMadaptive backend only).
#' @param ... passed to the backend fitting function.
#' @return a glmmTMB or MixMod fit.
#' @export
fit_multiproc <- function(long, fixed, id, shared_re = FALSE,
                          backend = c("glmmTMB", "GLMMadaptive"),
                          nAGQ = 11, ...) {
  backend <- match.arg(backend)
  stopifnot(inherits(fixed, "formula"), length(fixed) == 2L)
  rhs <- paste(deparse(fixed[[2]]), collapse = " ")
  one_proc <- nlevels(droplevels(long$process)) == 1L
  fix_rhs <- if (one_proc) rhs
             else sprintf("0 + process + process:(%s)", rhs)
  if (backend == "glmmTMB") {
    re <- if (shared_re || one_proc) sprintf("(1 | %s)", id)
          else sprintf("(0 + process | %s)", id)
    f <- as.formula(sprintf(".event ~ %s + %s + offset(log(.exposure))",
                            fix_rhs, re))
    glmmTMB::glmmTMB(f, family = stats::poisson(), data = long, ...)
  } else {
    if (!requireNamespace("GLMMadaptive", quietly = TRUE))
      stop("backend='GLMMadaptive' requires the GLMMadaptive package.")
    ff <- as.formula(sprintf(".event ~ %s + offset(log(.exposure))", fix_rhs))
    rf <- if (shared_re || one_proc) as.formula(sprintf("~ 1 | %s", id))
          else as.formula(sprintf("~ 0 + process | %s", id))
    GLMMadaptive::mixed_model(fixed = ff, random = rf,
                              family = stats::poisson(), data = long,
                              nAGQ = nAGQ, ...)
  }
}

#' Recover structural coefficients by the delta method (nlcom equivalent)
#'
#' Evaluates an arbitrary nonlinear combination of fixed-effect coefficients
#' and computes its standard error by the delta method, exactly as Stata's
#' `nlcom`. For Bartus's nonrecursive two-equation system, the structural
#' selection and covariate effects are
#' `lambda_j = pi_jk / pi_kk` and `beta_j = pi_j0 - lambda_j * pi_k0`
#' (Bartus 2017, eqs. 5-6); pass those expressions here.
#'
#' @param fit a fitted model with `fixef()`-style coefficients (glmmTMB) or a
#'   list with elements `coef` (named vector) and `vcov` (matrix).
#' @param expr a string: an R expression over coefficient names, which may be
#'   backtick-quoted, e.g.
#'   `"`processbirth2:hereduc3` - (`processbirth2:mardur`/`processdivorce:mardur`) * `processdivorce:hereduc3`"`.
#' @return data.frame with estimate, delta-method SE, z, p, and 95% CI.
#' @export
recover_structural <- function(fit, expr) {
  if (inherits(fit, "glmmTMB")) {
    b <- glmmTMB::fixef(fit)$cond
    V <- stats::vcov(fit)$cond
  } else if (inherits(fit, "MixMod")) {
    b <- GLMMadaptive::fixef(fit)
    V <- tryCatch(stats::vcov(fit, parm = "fixed-effects"),
                  error = function(e) stats::vcov(fit))
  } else {
    b <- fit$coef; V <- fit$vcov
  }
  V <- V[names(b), names(b)]
  g <- function(bb) {
    names(bb) <- names(b)
    eval(parse(text = expr), envir = as.list(bb))
  }
  est <- g(b)
  grad <- numDeriv::grad(g, b)
  se <- sqrt(drop(t(grad) %*% V %*% grad))
  z <- est / se
  data.frame(expr = expr, estimate = est, se = se, z = z,
             p = 2 * stats::pnorm(-abs(z)),
             ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
             row.names = NULL)
}
