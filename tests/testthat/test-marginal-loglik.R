test_that("pecme_marginal_loglik matches direct hand computation on a tiny toy example", {
  ## Regression test for a critical bug: pecme_marginal_loglik() (a
  ## separate, duplicate implementation of the same GH integration as
  ## estep.R's E-step) computed colSums(sweep(logL, 2L, log_w, "+")),
  ## which adds log_w[k] to every ROW of logL before summing over rows
  ## -- i.e. multiplies the Gauss-Hermite log-weight by the number of
  ## observations in the group, instead of adding it once per node.
  ## The correct quantity for a group with rows i in g and quadrature
  ## node z_k is log(w_k) + sum_i log P(y_i | b = z_k), added ONCE.
  y <- c(1.2, 0.8)     # one group, two observations
  X <- cbind(1, c(0, 0))
  group <- c(1L, 1L)
  beta <- c(1, 0)
  sigma <- 1
  sigma_b <- 1
  Q <- 5

  gh <- gauss_hermite_normal(Q)
  b_nodes <- sigma_b * gh$nodes
  log_w <- log(gh$weights)

  ## Hand-computed reference: for each node, sum the per-row log
  ## densities, add log(w_k) ONCE, then logsumexp over nodes.
  eta <- as.vector(X %*% beta)
  logL <- sapply(seq_along(b_nodes), function(k) {
    sum(dnorm(y, mean = eta + b_nodes[k], sd = sigma, log = TRUE))
  })
  ref <- logsumexp(logL + log_w)

  got <- pecme_marginal_loglik(y, X, group, beta, sigma, sigma_b,
                                c("none", "none"), c(NA, NA), c(NA, NA), quadrature = Q)

  expect_equal(got, ref, tolerance = 1e-10)

  ## What the OLD (buggy) formula would have produced, for contrast:
  ## adding log_w[k] once PER ROW (here, twice, since the group has 2
  ## observations) rather than once per node.
  buggy <- logsumexp(logL + 2 * log_w)
  expect_false(isTRUE(all.equal(got, buggy)))
})

test_that("pecme_marginal_loglik matches a brute-force numerical integration reference", {
  set.seed(654)
  J <- 8
  n_per <- 5
  beta <- c(0.5, -1)
  sigma <- 0.8
  sigma_b <- 0.6
  group <- rep(seq_len(J), each = n_per)
  x <- rnorm(J * n_per)
  X <- cbind(1, x)
  b_true <- rnorm(J, 0, sigma_b)
  y <- as.vector(X %*% beta) + b_true[group] + rnorm(J * n_per, 0, sigma)
  censor <- rep("none", length(y))
  left <- right <- rep(NA_real_, length(y))

  ll_gh <- pecme_marginal_loglik(y, X, group, beta, sigma, sigma_b, censor, left, right, quadrature = 80)

  ## Brute-force reference: numerically integrate each group's
  ## likelihood over b via stats::integrate().
  eta <- as.vector(X %*% beta)
  ll_ref <- 0
  for (g in seq_len(J)) {
    idx <- which(group == g)
    integrand <- function(b) {
      vapply(b, function(bb) {
        prod(dnorm(y[idx], mean = eta[idx] + bb, sd = sigma)) * dnorm(bb, 0, sigma_b)
      }, numeric(1))
    }
    val <- stats::integrate(integrand, lower = -10 * sigma_b, upper = 10 * sigma_b,
                             rel.tol = 1e-10, subdivisions = 500)$value
    ll_ref <- ll_ref + log(val)
  }

  expect_equal(ll_gh, ll_ref, tolerance = 1e-4)
})

test_that("pecme_marginal_loglik is insensitive to (does not spuriously scale with) group size", {
  ## A second, less precise but more intuitive check alongside the two
  ## above: replacing g near-identical low-noise replicate rows within
  ## a group with 2g such replicates should scale the total
  ## log-likelihood by close to a factor of 2 for a well-separated
  ## case, not by something contaminated by the (buggy) per-row GH
  ## log-weight multiplication.
  set.seed(321)
  Jg <- 15
  beta <- c(1, 2)
  sigma <- 0.01
  sigma_b <- 0.5
  b <- rnorm(Jg, 0, sigma_b)

  build <- function(reps_per_group) {
    group <- rep(seq_len(Jg), each = reps_per_group)
    x <- rep(rnorm(Jg), each = reps_per_group)
    X <- cbind(1, x)
    eta <- as.vector(X %*% beta) + b[group]
    y <- eta + rnorm(length(group), 0, sigma)
    list(y = y, X = X, group = group)
  }

  d1 <- build(4)
  d2 <- build(8)

  ll1 <- pecme_marginal_loglik(d1$y, d1$X, d1$group, beta, sigma, sigma_b,
                                rep("none", length(d1$y)), rep(NA, length(d1$y)), rep(NA, length(d1$y)),
                                quadrature = 60)
  ll2 <- pecme_marginal_loglik(d2$y, d2$X, d2$group, beta, sigma, sigma_b,
                                rep("none", length(d2$y)), rep(NA, length(d2$y)), rep(NA, length(d2$y)),
                                quadrature = 60)

  expect_lt(ll2 / ll1, 2.5)
  expect_gt(ll2 / ll1, 1.5)
})
