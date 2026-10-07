## -----------------------------------------------------------------------
## Fit-quality metrics
## -----------------------------------------------------------------------

#' Naive parameter count used for AIC/BIC/EBIC/GCV
#'
#' Counts nonzero fixed-effects coefficients plus two variance
#' components (`sigma`, `sigma_b`). This is a convenient, commonly used
#' HEURISTIC, not the effective degrees of freedom in the technical
#' sense the AIC/GCV theory assumes. It is a reasonable approximation
#' for Lasso-type penalties under regularity conditions (Zou, Hastie &
#' Tibshirani, 2007, show `E[#nonzero]` is an unbiased estimate of the
#' Lasso's effective df), but it is NOT generally accurate for Ridge
#' (which keeps essentially every coefficient nonzero, so this count is
#' close to the maximum `p` regardless of how much shrinkage is
#' actually applied -- see [pecme_effective_df()] for a proper,
#' closed-form Ridge alternative) nor for Elastic Net, SCAD, or MCP
#' (whose effective df would need a generalized-df/local-derivative
#' calculation, e.g. Zou et al. 2007's generalized df; not implemented
#' here). AIC/BIC/EBIC/GCV computed from this `df` should therefore be
#' read as approximate, especially for Ridge, Elastic Net, SCAD, and MCP.
#'
#' @keywords internal
pecme_df <- function(beta, has_intercept = TRUE, tol = 1e-8) {
  pen_index <- if (has_intercept) seq_along(beta)[-1L] else seq_along(beta)
  n_selected <- sum(abs(beta[pen_index]) > tol)
  n_selected + as.integer(has_intercept) + 2L  # + sigma + sigma_b
}

#' Working Ridge-type effective degrees-of-freedom approximation
#'
#' For Ridge, the naive nonzero-coefficient count from [pecme_df()] is a
#' poor stand-in for effective df. The classical Ridge smoother trace is
#' used as a CONDITIONAL/WORKING approximation on the standardized fixed-
#' effects design. It is not claimed to be the exact effective df of the
#' complete PECME estimator, which also contains censoring, random effects,
#' estimated variance components, quadrature, and possibly adaptive
#' penalization.
#'
#' No corresponding closed form is implemented here for Elastic Net,
#' SCAD, or MCP (a generalized-df/local-derivative calculation would be
#' needed, e.g. Zou, Hastie & Tibshirani, 2007); for `"none"`/`"lasso"`/
#' `"adaptive"` the naive nonzero count is already a standard, reasonable
#' approximation and is returned unchanged (plus the 2 variance
#' components).
#'
#' @param object A `"pecme"` fit.
#' @return A working/conditional df approximation, including the
#'   intercept if present and the two variance components, for `"ridge"`,
#'   `"none"`, `"lasso"`, or `"adaptive"` fits; `NA_real_` (with an explanatory
#'   `attr(, "note")`) for `"elastic"`, `"scad"`, or `"mcp"`, where no
#'   closed form is implemented.
#' @references
#' Hastie, T. and Tibshirani, R. (1990). *Generalized Additive Models*.
#' Chapman and Hall. (Ridge effective df / "trace of the hat matrix".)
#'
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
#' freedom" of the lasso. *Annals of Statistics*, 35(5), 2173-2192.
#' @export
pecme_effective_df <- function(object) {
  if (!inherits(object, "pecme")) stop("`object` must be a \"pecme\" fit.", call. = FALSE)

  if (object$penalty %in% c("none", "lasso", "adaptive") || object$lambda == 0) {
    return(pecme_df(object$coefficients, object$has_intercept))
  }
  if (object$penalty != "ridge") {
    return(structure(NA_real_, note = paste0(
      "No closed-form effective df is implemented for penalty = \"", object$penalty,
      "\" (would need a generalized-df/local-derivative calculation; see ",
      "Zou, Hastie & Tibshirani, 2007). Falls back to NA rather than the ",
      "misleading nonzero-coefficient count from pecme_df()."
    )))
  }

  std <- pecme_standardize(object$X, object$has_intercept)
  Xs <- std$X
  p <- ncol(Xs)
  D <- diag(p)
  if (object$has_intercept) D[1L, 1L] <- 0

  lambda_star <- object$lambda * object$sigma^2
  XtX <- crossprod(Xs)
  penalized <- XtX + lambda_star * D
  hat_trace <- sum(diag(solve(penalized, XtX)))

  hat_trace + 2  # + sigma, sigma_b (fixed effects df already includes the intercept)
}

#' @keywords internal
pecme_prediction_metrics <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  yt <- y_true[ok]; yp <- y_pred[ok]
  resid <- yt - yp
  n <- length(yt)

  mse  <- mean(resid^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(resid))

  nz <- yt != 0
  mape  <- if (any(nz)) mean(abs(resid[nz] / yt[nz])) * 100 else NA_real_
  denom <- (abs(yt) + abs(yp))
  sd_ok <- denom > .Machine$double.eps
  smape <- if (any(sd_ok)) mean(2 * abs(resid[sd_ok]) / denom[sd_ok]) * 100 else NA_real_

  sst <- sum((yt - mean(yt))^2)
  r2  <- if (sst > 0) 1 - sum(resid^2) / sst else NA_real_

  list(n = n, MSE = mse, RMSE = rmse, MAE = mae, MAPE = mape, SMAPE = smape, R2 = r2)
}

#' @keywords internal
pecme_adjusted_r2 <- function(r2, n, df) {
  if (!is.finite(r2) || n - df - 1 <= 0) return(NA_real_)
  1 - (1 - r2) * (n - 1) / (n - df - 1)
}

#' Information criteria (AIC/BIC/EBIC/GCV)
#'
#' AIC/BIC/GCV use `df`, the naive nonzero-coefficient-based
#' parameter count of [pecme_df()] -- a heuristic, not a general
#' effective-degrees-of-freedom calculation (see [pecme_df()] and
#' [pecme_effective_df()]). This matters most for Ridge (essentially no
#' coefficients are ever exactly zero, so `df` sits near the maximum
#' regardless of shrinkage strength) and, to a lesser extent, Elastic
#' Net/SCAD/MCP; treat AIC/BIC/EBIC/GCV from those fits as approximate,
#' and prefer [pecme_effective_df()]-based comparisons for Ridge, or
#' cross-validated criteria (`criterion = "cv"` in [grid_pecme()],
#' [pecme_loo()]) when the choice of `df` matters for the conclusion.
#'
#' @keywords internal
pecme_ic <- function(logLik, df, n, p_total, gamma_ebic = 0.5,
                      n_selected = NULL) {
  aic <- -2 * logLik + 2 * df
  bic <- -2 * logLik + df * log(n)
  ## Chen & Chen (2008): the combinatorial EBIC term is based on the
  ## number of selected PENALIZED predictors k, not on total df minus
  ## two. The latter can accidentally retain the intercept.
  if (is.null(n_selected)) n_selected <- max(df - 2L, 0L)
  n_selected <- as.integer(n_selected)
  log_choose <- if (p_total > 0L && n_selected >= 0L && n_selected <= p_total) {
    lchoose(p_total, n_selected)
  } else 0
  ebic <- bic + 2 * gamma_ebic * log_choose
  gcv <- if (n > df) (-2 * logLik) / (n * (1 - df / n)^2) else Inf
  list(AIC = aic, BIC = bic, EBIC = ebic, GCV = gcv)
}
