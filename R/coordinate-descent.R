## -----------------------------------------------------------------------
## Coordinate descent for the CM-step of beta
##
##   Q_beta(beta) = ||y - X beta||^2 / (2 sigma^2) + P_lambda(beta)
##
## `y` here is the pseudo-response E[Y*] - E[b] from the E-step (see
## objective.R for why this is exactly equivalent, up to a beta-free
## constant, to minimizing the true conditional expected complete-data
## objective). The intercept, if present, is never penalized and is
## updated in closed form as a plain mean of the partial residual.
## -----------------------------------------------------------------------

#' @keywords internal
pecme_coordinate_descent <- function(X, y, beta, sigma, lambda,
                                      penalty = .pecme_penalties,
                                      weights = NULL, alpha = 0.5, gamma = 3.7,
                                      max_iter = 1000L, tol = 1e-7,
                                      monotone_tol = 1e-10,
                                      check_monotonicity = TRUE,
                                      has_intercept = TRUE) {
  penalty <- match.arg(penalty)
  X <- as.matrix(X); y <- as.numeric(y); beta <- as.numeric(beta)
  n <- nrow(X); p <- ncol(X)

  if (length(y) != n) stop("length(y) must equal nrow(X).", call. = FALSE)
  if (length(beta) != p) stop("length(beta) must equal ncol(X).", call. = FALSE)
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be finite and > 0.", call. = FALSE)
  if (!is.finite(lambda) || lambda < 0) stop("lambda must be finite and >= 0.", call. = FALSE)
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("alpha must be in [0, 1].", call. = FALSE)
  if (penalty == "scad" && (!is.finite(gamma) || gamma <= 2)) stop("gamma must be > 2 for SCAD.", call. = FALSE)
  if (penalty == "mcp"  && (!is.finite(gamma) || gamma <= 1)) stop("gamma must be > 1 for MCP.", call. = FALSE)

  if (p == 0L) {
    return(list(beta = beta, iter = 0L, converged = TRUE, reason = "no_columns",
                max_change = 0, objective_before = NA_real_, objective_after = NA_real_,
                objective_history = numeric(0), monotone = TRUE, lambda_star = lambda * sigma^2))
  }

  if (penalty == "adaptive") {
    n_pen <- if (has_intercept) p - 1L else p
    if (is.null(weights) || length(weights) != n_pen) {
      stop("Adaptive LASSO requires `weights` of length ", n_pen, ".", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("Adaptive LASSO weights must be finite and non-negative.", call. = FALSE)
    }
  }

  cond_obj <- function(b) {
    r <- y - as.vector(X %*% b)
    sum(r^2) / (2 * sigma^2) +
      penalty_value(b, penalty, lambda, alpha, gamma, weights, has_intercept)
  }

  obj_before <- cond_obj(beta)
  if (!is.finite(obj_before)) stop("Initial conditional objective is not finite.", call. = FALSE)

  obj_hist <- numeric(max_iter)
  monotone <- TRUE
  converged <- FALSE
  reason <- "max_iterations"
  max_change <- NA_real_
  start_j <- if (has_intercept) 2L else 1L

  for (iter in seq_len(max_iter)) {
    beta_old <- beta

    if (has_intercept) {
      beta[1L] <- mean(y - X[, -1L, drop = FALSE] %*% beta[-1L])
    }

    r <- y - as.vector(X %*% beta)

    if (p >= start_j) {
      for (j in start_j:p) {
        r_j <- r + X[, j] * beta[j]
        zj <- sum(X[, j] * r_j)
        vj <- sum(X[, j]^2)

        if (!is.finite(vj) || vj <= .Machine$double.eps) {
          beta[j] <- 0
          r <- r_j
          next
        }

        wj <- if (penalty == "adaptive") weights[j - start_j + 1L] else 1
        beta_new <- pecme_coordinate_update(zj, vj, sigma, lambda, penalty, alpha, gamma, wj)

        if (!is.finite(beta_new)) {
          stop("Non-finite coordinate update at j = ", j, ".", call. = FALSE)
        }

        r <- r_j - X[, j] * beta_new
        beta[j] <- beta_new
      }
    }

    obj_new <- cond_obj(beta)
    if (!is.finite(obj_new)) {
      reason <- "nonfinite_objective"
      stop("Coordinate descent produced a non-finite objective.", call. = FALSE)
    }
    obj_hist[iter] <- obj_new

    if (iter > 1L) {
      prev <- obj_hist[iter - 1L]
      rel_increase <- (obj_new - prev) / max(abs(prev), 1)
      if (obj_new > prev + monotone_tol && rel_increase > 1e-8) {
        monotone <- FALSE
        if (check_monotonicity) {
          reason <- "objective_increase"
          stop(sprintf("Coordinate descent increased the conditional objective at sweep %d: %.12g -> %.12g.",
                        iter, prev, obj_new), call. = FALSE)
        }
      }
    }

    max_change <- max(abs(beta - beta_old))
    rel_change <- max_change / max(max(abs(beta_old)), 1)
    if (!is.finite(max_change) || !is.finite(rel_change)) {
      reason <- "nonfinite_coefficients"
      stop("Non-finite coefficient update detected.", call. = FALSE)
    }

    if (max_change < tol && rel_change < tol) {
      converged <- TRUE
      reason <- "converged"
      break
    }
  }

  list(
    beta = beta, iter = iter, converged = converged, reason = reason,
    max_change = max_change, objective_before = obj_before,
    objective_after = cond_obj(beta), objective_history = obj_hist[seq_len(iter)],
    monotone = monotone, lambda_star = lambda * sigma^2
  )
}
