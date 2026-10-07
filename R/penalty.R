## -----------------------------------------------------------------------
## Penalty definitions
##
## All penalties are applied only to the columns in `pen_index` (i.e.
## excluding the intercept, when present). A single source of truth is
## kept here for (a) the penalty VALUE, used in the objective/BIC/etc,
## and (b) the coordinate-wise closed-form UPDATE used by coordinate
## descent, so the two can never drift apart as they did in the
## original draft (which reimplemented MCP inline inside the
## coordinate-descent loop instead of calling the shared helper).
## -----------------------------------------------------------------------

.pecme_penalties <- c("none", "lasso", "adaptive", "scad", "mcp", "elastic", "ridge")

#' Value of the penalty term P_lambda(beta)
#'
#' @param beta Full coefficient vector (including the intercept, if any).
#' @param has_intercept Whether `beta[1]` is an unpenalized intercept.
#' @keywords internal
penalty_value <- function(beta, penalty = .pecme_penalties, lambda,
                           alpha = 0.5, gamma = 3.7, weights = NULL,
                           has_intercept = TRUE) {
  penalty <- match.arg(penalty)
  pen_index <- if (has_intercept) seq_along(beta)[-1L] else seq_along(beta)
  bp <- as.numeric(beta[pen_index])

  if (!is.finite(lambda) || lambda < 0) stop("lambda must be finite and >= 0.", call. = FALSE)
  if (penalty == "none" || lambda == 0) return(0)

  switch(
    penalty,
    lasso = lambda * sum(abs(bp)),
    adaptive = {
      if (is.null(weights) || length(weights) != length(bp)) {
        stop("Adaptive LASSO requires `weights` of length = number of penalized coefficients.", call. = FALSE)
      }
      lambda * sum(weights * abs(bp))
    },
    elastic = lambda * (alpha * sum(abs(bp)) + (1 - alpha) / 2 * sum(bp^2)),
    ridge = lambda / 2 * sum(bp^2),
    mcp = {
      if (!is.finite(gamma) || gamma <= 1) stop("gamma must be > 1 for MCP.", call. = FALSE)
      ab <- abs(bp)
      idx <- ab <= gamma * lambda
      pen <- numeric(length(ab))
      pen[idx]  <- lambda * ab[idx] - ab[idx]^2 / (2 * gamma)
      pen[!idx] <- gamma * lambda^2 / 2
      sum(pen)
    },
    scad = {
      if (!is.finite(gamma) || gamma <= 2) stop("gamma must be > 2 for SCAD.", call. = FALSE)
      ab <- abs(bp)
      pen <- numeric(length(ab))
      i1 <- ab <= lambda
      i2 <- ab > lambda & ab <= gamma * lambda
      i3 <- ab > gamma * lambda
      pen[i1] <- lambda * ab[i1]
      pen[i2] <- (-ab[i2]^2 + 2 * gamma * lambda * ab[i2] - lambda^2) / (2 * (gamma - 1))
      pen[i3] <- (gamma + 1) * lambda^2 / 2
      sum(pen)
    }
  )
}

#' Soft-thresholding operator
#' @keywords internal
soft_threshold <- function(z, a) sign(z) * pmax(abs(z) - a, 0)

#' Closed-form coordinate update for the penalized quadratic subproblem
#'
#' Solves, in the single coordinate `beta_j`,
#' \deqn{\min_{\beta_j} \frac{v_j \beta_j^2 - 2 z_j \beta_j}{2\sigma^2} +
#'   p_\lambda(|\beta_j|)}
#' for each of the supported penalties, where `z_j` is the partial
#' cross-product of column `j` with the partial residual and `v_j` is
#' its sum of squares. `lambda_star = lambda * sigma^2` rescales the
#' penalty onto the same footing as the (unscaled) least-squares term.
#'
#' @keywords internal
pecme_coordinate_update <- function(zj, vj, sigma, lambda, penalty,
                                     alpha = 0.5, gamma = 3.7, wj = 1) {
  lambda_star <- lambda * sigma^2

  if (penalty == "none") return(zj / vj)

  if (penalty == "lasso")    return(soft_threshold(zj, lambda_star) / vj)
  if (penalty == "adaptive") return(soft_threshold(zj, lambda_star * wj) / vj)
  if (penalty == "elastic") {
    return(soft_threshold(zj, lambda_star * alpha) / (vj + lambda_star * (1 - alpha)))
  }
  if (penalty == "ridge") return(zj / (vj + lambda_star))

  if (penalty == "mcp") {
    denom <- vj - sigma^2 / gamma
    if (denom <= 0) {
      stop("MCP coordinate subproblem is locally non-convex (vj - sigma^2/gamma <= 0). ",
           "Increase `gamma`, reduce `lambda`, or standardize predictors.", call. = FALSE)
    }
    if (abs(zj) <= gamma * lambda * vj) {
      return(soft_threshold(zj, lambda_star) / denom)
    }
    return(zj / vj)
  }

  if (penalty == "scad") {
    denom <- (gamma - 1) * vj - sigma^2
    if (denom <= 0) {
      stop("SCAD coordinate subproblem is locally non-convex ((gamma-1)*vj - sigma^2 <= 0). ",
           "Increase `gamma`, reduce `lambda`, or standardize predictors.", call. = FALSE)
    }
    s <- abs(zj)
    b1 <- lambda * (vj + sigma^2)
    b2 <- gamma * lambda * vj
    if (s <= b1) return(soft_threshold(zj, lambda_star) / vj)
    if (s <= b2) return(((gamma - 1) * zj - gamma * lambda * sigma^2 * sign(zj)) / denom)
    return(zj / vj)
  }

  stop("Unknown penalty: ", penalty, call. = FALSE)
}
