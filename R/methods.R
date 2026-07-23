# S3 interface for gensem_recursive objects: what a user sees in the first
# five minutes. print() is brief; summary() carries the full inferential
# table; coef/vcov/confint/logLik follow stats conventions.

#' @export
coef.gensem_recursive <- function(object, ...) object$coef

#' @export
vcov.gensem_recursive <- function(object, ...) object$vcov

#' @export
confint.gensem_recursive <- function(object, parm, level = 0.95, ...) {
  co <- object$coef; se <- sqrt(diag(object$vcov))
  if (missing(parm)) parm <- names(co)
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(co[parm] - z * se[parm], co[parm] + z * se[parm])
  colnames(out) <- sprintf("%.1f %%", 100 * c((1 - level) / 2,
                                              1 - (1 - level) / 2))
  out
}

#' @export
#' @note The stored objective drops additive data constants (e.g. log(y!)
#'   in the Poisson kernel), so logLik values are comparable only across
#'   models fit to the SAME data under the SAME design -- valid for LR
#'   tests of nested models, not for absolute comparisons.
logLik.gensem_recursive <- function(object, ...) {
  structure(-object$nll, df = length(object$coef), class = "logLik")
}

#' @export
print.gensem_recursive <- function(x, digits = 4, ...) {
  cat("Recursive multiprocess system (GenSEM)\n")
  cat(sprintf("Hazard family: %s | quadrature Q = %d | convergence: %s\n",
              x$family, x$Q,
              if (x$convergence == 0) "yes" else paste0("NO (code ",
                                                        x$convergence, ")")))
  co <- x$coef; se <- sqrt(diag(x$vcov))
  print(round(cbind(Estimate = co, `Std.Error` = se), digits))
  cat("Random-effect covariance (person level):\n")
  print(round(x$Sigma, digits))
  if (!is.null(x$sw2))
    cat(sprintf("Nested-cluster frailty variance: %.4f\n", x$sw2))
  invisible(x)
}

#' @export
summary.gensem_recursive <- function(object, ...) {
  co <- object$coef; se <- sqrt(diag(object$vcov))
  z <- co / se
  tab <- cbind(Estimate = co, `Std.Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(coefficients = tab, Sigma = object$Sigma,
                 family = object$family, ancillary = object$ancillary,
                 sw2 = object$sw2, nll = object$nll, Q = object$Q,
                 convergence = object$convergence),
            class = "summary.gensem_recursive")
}

#' @export
print.summary.gensem_recursive <- function(x, digits = 4, ...) {
  cat("Recursive multiprocess system (GenSEM)\n")
  cat(sprintf("Hazard family: %s | quadrature Q = %d | convergence: %s\n",
              x$family, x$Q, if (x$convergence == 0) "yes" else "NO"))
  cat(sprintf("Negative marginal log-likelihood (constants dropped): %.3f\n\n",
              x$nll))
  stats::printCoefmat(x$coefficients, digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE)
  cat("\nRandom-effect covariance (person level):\n")
  print(round(x$Sigma, digits))
  if (!is.null(x$ancillary))
    cat(sprintf("Ancillary (shape/scale/sigma): %.4f\n", x$ancillary))
  if (!is.null(x$sw2))
    cat(sprintf("Nested-cluster frailty variance: %.4f\n", x$sw2))
  invisible(x)
}
