## -----------------------------------------------------------------------
## Data-dependent lambda grid
##
## lambda_max is a data-dependent reference scale based on the
## zero-coefficient/one-step score condition. It is NOT claimed to be a
## formally exact sparsity threshold for the complete censored
## random-effects estimator. Lasso, Adaptive Lasso, SCAD, and MCP all share the
## same zero-boundary condition (their penalty derivative at 0+ equals
## lambda, or lambda*weight_j for Adaptive Lasso, or lambda*alpha for
## Elastic Net), so a single formula covers all of them. Ridge and
## "none" have no exact sparsity threshold; the same value is used
## only as a reference scale for building a comparable grid.
## -----------------------------------------------------------------------

#' @keywords internal
pecme_lambda_max <- function(X, y, group, sigma, sigma_b, censor, left, right,
                              penalty = .pecme_penalties, alpha = 0.5, weights = NULL,
                              has_intercept = TRUE, quadrature = 40L) {
  penalty <- match.arg(penalty)

  p <- ncol(X)

  pen_index <- if (has_intercept) {
    seq_len(p)[-1L]
  } else {
    seq_len(p)
  }

  if (!length(pen_index)) {
    return(0)
  }

  ## For censored mixed models, construct the E-step pseudo-response.
  ## The raw response may contain NA for censored observations.
  beta0 <- rep(0, p)

  moments <- pecme_e_step(
    y = y,
    X = X,
    beta = beta0,
    group = group,
    sigma = sigma,
    sigma_b = sigma_b,
    censor = censor,
    left = left,
    right = right,
    quadrature = quadrature
  )

  y_adj <- as.numeric(moments$Ey - moments$Eb)

  base_line <- if (has_intercept) {
    mean(y_adj)
  } else {
    0
  }

  z <- as.vector(
    crossprod(
      X[, pen_index, drop = FALSE],
      y_adj - base_line
    )
  )

  eff_w <- switch(
    penalty,

    lasso = ,
    scad = ,
    mcp = rep(1, length(z)),

    adaptive = {
      if (is.null(weights) || length(weights) != length(z)) {
        stop(
          "Adaptive LASSO requires `weights` to compute lambda_max.",
          call. = FALSE
        )
      }

      pmax(weights, 1e-8)
    },

    elastic = rep(max(alpha, 1e-8), length(z)),

    ridge = ,
    none = rep(1, length(z))
  )

  lam_max <- max(abs(z) / eff_w) / sigma^2

  max(lam_max, 1e-8)
}


#' Build a log-spaced lambda grid
#'
#' @param lambda_max Upper end of the grid.
#' @param nlambda Number of grid points.
#' @param lambda_min_ratio `lambda_min = lambda_min_ratio * lambda_max`.
#'   Defaults follow the usual `n > p` vs `n <= p` heuristic (`1e-4` /
#'   `1e-2` respectively) unless supplied explicitly.
#' @param n Number of observations.
#' @param p Number of predictors.
#' @return A numeric vector of length `nlambda` containing a decreasing
#'   log-spaced sequence from `lambda_max` to `lambda_min`.
#' @export
pecme_lambda_grid <- function(lambda_max, nlambda = 50L, lambda_min_ratio = NULL,
                               n = NULL, p = NULL) {
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("`lambda_max` must be finite and > 0.", call. = FALSE)
  }
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- if (!is.null(n) && !is.null(p) && p >= n) 1e-2 else 1e-4
  }
  lambda_min <- lambda_min_ratio * lambda_max
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}
