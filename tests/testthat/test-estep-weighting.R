test_that("pecme_e_step's posterior weights are not group-size-contaminated", {
  ## The E-step computes log_post_group <- colSums(logL) + log_w (added
  ## once per quadrature node). An earlier version used
  ## colSums(sweep(logL, 2L, log_w, "+")), which adds log_w[k] once per
  ## OBSERVATION in the group -- i.e. multiplies it by group size.
  ## Consequence of the bug: the posterior distribution over b (and
  ## hence E[b | data]) would depend spuriously on group size even when
  ## it shouldn't. Check by comparing the E-step's posterior mean of b
  ## for a group against a brute-force numerical-integration reference.
  set.seed(789)
  sigma <- 0.7
  sigma_b <- 0.9
  beta <- c(0.3, -0.6)
  n_i <- 6
  x <- rnorm(n_i)
  X <- cbind(1, x)
  eta <- as.vector(X %*% beta)
  b_true <- 0.4
  y <- eta + b_true + rnorm(n_i, 0, sigma)
  group <- rep(1L, n_i)
  censor <- rep("none", n_i)

  e <- pecme_e_step(y, X, beta, group, sigma, sigma_b, censor,
                     rep(NA, n_i), rep(NA, n_i), quadrature = 80)

  ## Brute-force posterior mean of b via numerical integration.
  log_post <- function(b) {
    sum(dnorm(y, mean = eta + b, sd = sigma, log = TRUE)) + dnorm(b, 0, sigma_b, log = TRUE)
  }
  num <- stats::integrate(function(bb) vapply(bb, function(b) b * exp(log_post(b)), numeric(1)),
                           -10 * sigma_b, 10 * sigma_b, subdivisions = 500, rel.tol = 1e-10)$value
  den <- stats::integrate(function(bb) vapply(bb, function(b) exp(log_post(b)), numeric(1)),
                           -10 * sigma_b, 10 * sigma_b, subdivisions = 500, rel.tol = 1e-10)$value
  eb_ref <- num / den

  expect_equal(unique(e$group_Eb), eb_ref, tolerance = 1e-4)
})
