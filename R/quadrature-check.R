## -----------------------------------------------------------------------
## pecme_quadrature_check(): sensitivity of the fit to the number of
## Gauss-Hermite quadrature nodes.
##
## The random intercept is integrated out NUMERICALLY (fixed-order
## Gauss-Hermite quadrature approximates the continuous integral; it
## does not evaluate it exactly -- see ?pecme). This diagnostic makes
## that approximation error visible instead of leaving `quadrature`'s
## default as an unexamined assumption.
## -----------------------------------------------------------------------

#' Check sensitivity of a fitted model to the Gauss-Hermite quadrature order
#'
#' Fixed-order Gauss-Hermite quadrature is a numerical approximation to
#' the (otherwise intractable) integral over the random intercept, not
#' an exact evaluation of it. This function reports how much the fit
#' changes as the number of quadrature nodes changes, so `quadrature`
#' is a checked choice rather than an unexamined default. Approximation
#' error is typically largest when `sigma_b` is large relative to
#' `sigma`, or with heavy censoring/many missing responses.
#'
#' @param object A `"pecme"` fit.
#' @param quadrature_values Integer vector of node counts to compare
#'   against `object`'s own `quadrature`, which is always included.
#' @param method `"loglik"` (default; cheap: re-evaluates the marginal
#'   log-likelihood at `object`'s already-fitted `(beta, sigma,
#'   sigma_b)` under each node count, with no re-optimization) or
#'   `"refit"` (expensive: fully refits the model -- including its
#'   `beta`, `sigma`, `sigma_b` -- at each node count, showing whether
#'   the actual parameter *estimates*, not just the likelihood value at
#'   a fixed point, are sensitive to the quadrature order).
#' @return An object of class `"pecme_quadrature_check"`: a data frame
#'   (`quadrature`, `logLik`, and, for `method = "refit"`, `sigma`,
#'   `sigma_b`, `max_abs_beta_change` relative to `object`), with a
#'   `print` method that flags whether the range of `logLik` values
#'   looks small relative to typical evidence thresholds.
#' @examples
#' \donttest{
#' sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 3, seed = 1)
#' fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data)
#' pecme_quadrature_check(fit)
#' }
#' @export
pecme_quadrature_check <- function(object, quadrature_values = c(15L, 20L, 30L, 50L, 80L, 120L),
                                    method = c("loglik", "refit")) {
  if (!inherits(object, "pecme")) stop("`object` must be a \"pecme\" fit.", call. = FALSE)
  method <- match.arg(method)
  qs <- sort(unique(c(as.integer(quadrature_values), as.integer(object$quadrature))))

  if (method == "loglik") {
    ll <- vapply(qs, function(q) {
      pecme_marginal_loglik(object$y, object$X, object$group, object$coefficients,
                             object$sigma, object$sigma_b, object$censor,
                             object$left, object$right, quadrature = q)
    }, numeric(1))
    out <- data.frame(quadrature = qs, logLik = ll)
  } else {
    ## Reuse the exact original analysis data. Reconstructing a data frame
    ## from X can destroy factor levels, transformations, interactions,
    ## splines, and other formula semantics.
    if (is.null(object$analysis_data)) {
      stop("The fitted object does not contain original analysis data required for formula-preserving refit.",
           call. = FALSE)
    }
    d <- object$analysis_data

    rows <- lapply(qs, function(q) {
      f <- tryCatch(
        pecme(object$fixed_formula, random_var = object$random_var, data = d,
              censor = object$censor, left = object$left, right = object$right,
              penalty = object$penalty, lambda = object$lambda, alpha = object$alpha,
              gamma = object$gamma, adaptive_gamma = object$adaptive_gamma,
              missing_y = object$missing_y, missing_x = object$missing_x,
              standardize = object$standardize, quadrature = q,
              variance_step = object$variance_step, refit = FALSE,
              init_method = object$init_method, workers = object$workers),
        error = function(e) NULL
      )
      if (is.null(f)) {
        return(data.frame(quadrature = q, logLik = NA_real_, sigma = NA_real_,
                           sigma_b = NA_real_, max_abs_beta_change = NA_real_))
      }
      data.frame(quadrature = q, logLik = f$logLik, sigma = f$sigma, sigma_b = f$sigma_b,
                 max_abs_beta_change = max(abs(f$coefficients - object$coefficients)))
    })
    out <- do.call(rbind, rows)
  }

  structure(list(table = out, method = method, reference_quadrature = object$quadrature),
            class = "pecme_quadrature_check")
}

#' @export
print.pecme_quadrature_check <- function(x, digits = 4L, ...) {
  cat("Gauss-Hermite quadrature sensitivity check (method = \"", x$method, "\")\n", sep = "")
  cat("Reference (fitted) quadrature order:", x$reference_quadrature, "\n\n")
  tab <- x$table
  tab[-1] <- lapply(tab[-1], function(col) if (is.numeric(col)) round(col, digits) else col)
  print(tab, row.names = FALSE)

  ll_range <- diff(range(x$table$logLik, na.rm = TRUE))
  cat(sprintf("\nRange of logLik across quadrature orders: %.*f\n", digits, ll_range))
  if (is.finite(ll_range) && ll_range > 1) {
    cat("This is larger than a typical 'negligible' threshold (~1 log-likelihood unit): ",
        "consider increasing `quadrature` in the fitted model.\n", sep = "")
  } else if (is.finite(ll_range)) {
    cat("This looks small: the fit does not appear sensitive to the quadrature order ",
        "over the range checked.\n", sep = "")
  }
  invisible(x)
}
