## -----------------------------------------------------------------------
## Objective functions
## -----------------------------------------------------------------------

#' Penalized observed-data objective: -logLik(beta, sigma, sigma_b) + P_lambda(beta)
#' @keywords internal
pecme_objective <- function(beta, y, X, sigma, sigma_b, group, censor,
                             left = NULL, right = NULL, penalty, lambda,
                             alpha = 0.5, gamma = 3.7, weights = NULL,
                             has_intercept = TRUE, quadrature = 40L) {
  ll <- pecme_marginal_loglik(y, X, group, beta, sigma, sigma_b,
                               censor, left, right, quadrature)
  pen <- penalty_value(beta, penalty, lambda, alpha, gamma, weights, has_intercept)
  list(objective = -ll + pen, logLik = ll, penalty = pen)
}

#' Conditional expected complete-data objective for beta (the CM-step target)
#'
#' Given E-step moments, the beta-dependent part of the expected
#' complete-data penalized log-likelihood reduces to
#' \deqn{Q(\beta) = \frac{1}{2\sigma^2} \sum_i (E[y^*_i] - E[b_{g(i)}] -
#'   x_i'\beta)^2 + P_\lambda(\beta) + \text{const}}
#' because \eqn{E\|y^* - b - X\beta\|^2 = \mathrm{Var}(y^*-b\mid
#' \text{data}) + (E[y^*-b\mid\text{data}] - X\beta)^2} and the variance
#' term does not depend on beta. This is what makes coordinate descent
#' on the adjusted response `E[y*] - E[b]` exact rather than an
#' approximation.
#' @keywords internal
pecme_q_beta <- function(beta, y_star, X, random_effect = NULL, sigma,
                          penalty = .pecme_penalties, lambda = 0, alpha = 0.5,
                          gamma = 3.7, weights = NULL, has_intercept = TRUE) {
  penalty <- match.arg(penalty)
  beta <- as.numeric(beta); y_star <- as.numeric(y_star); X <- as.matrix(X)
  if (is.null(random_effect)) random_effect <- rep(0, nrow(X))

  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be > 0.", call. = FALSE)
  if (length(beta) != ncol(X)) stop("beta length mismatch.", call. = FALSE)

  resid <- y_star - random_effect - as.vector(X %*% beta)
  loss <- sum(resid^2) / (2 * sigma^2)
  pen <- penalty_value(beta, penalty, lambda, alpha, gamma, weights, has_intercept)
  loss + pen
}
