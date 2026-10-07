## -----------------------------------------------------------------------
## Observed-data log-likelihoods
## -----------------------------------------------------------------------

#' Conditional (given the random intercept) observed-data log-likelihood
#'
#' Sum over observations of `log P(observation | mu, sigma)`, where
#' `mu` already includes the random intercept. Used mainly for testing.
#'
#' @keywords internal
pecme_observed_loglik <- function(y, mu, sigma, censor, left = NULL, right = NULL) {
  n <- length(censor)
  if (is.null(left))  left  <- rep(NA_real_, n)
  if (is.null(right)) right <- rep(NA_real_, n)
  ll <- pecme_conditional_loglik(
    y = y, mu = matrix(mu, ncol = 1L), sigma = sigma,
    censor = censor, left = left, right = right
  )
  sum(ll)
}

#' Marginal observed-data log-likelihood, random intercept integrated out
#'
#' \deqn{\ell(\beta,\sigma,\sigma_b) = \sum_{j=1}^J \log \int
#'   \prod_{i \in j} P(y_{ij} \mid x_{ij}'\beta + b, \sigma)\,
#'   \phi(b; 0, \sigma_b^2)\, db}
#'
#' evaluated by Gauss-Hermite quadrature. This is the objective that the
#' ECME step maximizes (over `sigma`, `sigma_b`) and that the penalized
#' criterion in [pecme()] is built from (over `beta`, given `sigma`,
#' `sigma_b`).
#'
#' @inheritParams pecme_e_step
#' @keywords internal
pecme_marginal_loglik <- function(y, X, group, beta, sigma, sigma_b,
                                   censor, left = NULL, right = NULL,
                                   quadrature = 40L) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(sigma_b) || sigma_b <= 0) {
    stop("sigma and sigma_b must be positive and finite.", call. = FALSE)
  }
  X <- as.matrix(X); beta <- as.numeric(beta)
  y <- as.numeric(y); group <- as.integer(group); censor <- as.character(censor)
  n <- nrow(X)
  if (is.null(left))  left  <- rep(NA_real_, n)
  if (is.null(right)) right <- rep(NA_real_, n)

  gh <- gauss_hermite_normal(quadrature)
  b_nodes <- sigma_b * gh$nodes
  log_w <- log(gh$weights)

  eta <- as.vector(X %*% beta)
  groups <- sort(unique(group))
  total <- 0

  for (g in groups) {
    idx <- which(group == g)
    mu <- outer(eta[idx], b_nodes, "+")
    logL <- pecme_conditional_loglik(
      y = y[idx], mu = mu, sigma = sigma,
      censor = censor[idx], left = left[idx], right = right[idx]
    )
    ## NOTE: must add log_w ONCE per quadrature node (column), not once
    ## per observation (row) -- `colSums(sweep(logL, 2L, log_w, "+"))`
    ## adds log_w[k] to every row before summing, i.e. multiplies it by
    ## the group size, over-weighting the GH prior by a factor equal to
    ## the number of observations in the group (worse for larger
    ## groups). This mirrors a fix already applied to estep.R's
    ## sequential E-step; pecme_marginal_loglik() is a separate,
    ## duplicate implementation of the same integration and had not yet
    ## been corrected to match.
    log_integrand <- colSums(logL) + log_w
    total <- total + logsumexp(log_integrand)
  }

  as.numeric(total)
}
