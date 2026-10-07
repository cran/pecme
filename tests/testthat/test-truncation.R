test_that("uncensored ('none') moments just return y and y^2", {
  mu <- matrix(c(0, 1, -1), ncol = 1)
  y <- c(2, -0.5, 3)
  m <- pecme_truncated_moments(mu, sigma = 1, censor = c("none", "none", "none"),
                                left = rep(NA, 3), right = rep(NA, 3), y = y)
  expect_equal(as.vector(m$mean), y)
  expect_equal(as.vector(m$second), y^2)
})

test_that("'missing' moments reduce to the unconstrained normal mean/variance", {
  mu <- matrix(c(0.5, -1.2), ncol = 1)
  sigma <- 2
  m <- pecme_truncated_moments(mu, sigma = sigma, censor = c("missing", "missing"),
                                left = rep(NA, 2), right = rep(NA, 2), y = rep(NA, 2))
  expect_equal(as.vector(m$mean), as.vector(mu))
  expect_equal(as.vector(m$second), as.vector(mu)^2 + sigma^2)
})

test_that("'missing' rows contribute exactly 0 to the conditional log-likelihood", {
  mu <- matrix(c(0.5, -1.2, 10), ncol = 1)
  ll <- pecme_conditional_loglik(y = rep(NA, 3), mu = mu, sigma = 1.5,
                                  censor = rep("missing", 3), left = rep(NA, 3), right = rep(NA, 3))
  expect_equal(as.vector(ll), c(0, 0, 0))
})

test_that("left/right censoring log-likelihood matches pnorm at the boundary", {
  mu <- matrix(0, 1, 1)
  ll_left <- pecme_conditional_loglik(y = NA, mu = mu, sigma = 1, censor = "left",
                                       left = 1, right = NA)
  expect_equal(as.vector(ll_left), stats::pnorm(1, log.p = TRUE))

  ll_right <- pecme_conditional_loglik(y = NA, mu = mu, sigma = 1, censor = "right",
                                        left = NA, right = -1)
  expect_equal(as.vector(ll_right), stats::pnorm(-1, lower.tail = FALSE, log.p = TRUE))
})

test_that("interval censoring collapses to left/right in the limit", {
  ## As right -> Inf, interval(left, right) should match left-censoring's
  ## complement region... check log P(left < Y) matches survival at left.
  mu <- matrix(0, 1, 1)
  ll_int <- pecme_conditional_loglik(y = NA, mu = mu, sigma = 1, censor = "interval",
                                      left = 0.3, right = 1e6)
  ll_right_tail <- stats::pnorm(0.3, lower.tail = FALSE, log.p = TRUE)
  expect_equal(as.vector(ll_int), ll_right_tail, tolerance = 1e-6)
})

test_that("truncated first/second moments are consistent (variance >= 0)", {
  set.seed(2)
  mu <- matrix(stats::rnorm(20), ncol = 1)
  left <- mu[, 1] - abs(stats::rnorm(20))
  right <- left + abs(stats::rnorm(20)) + 0.1
  m <- pecme_truncated_moments(mu, sigma = 1, censor = rep("interval", 20),
                                left = left, right = right, y = rep(NA, 20))
  var_est <- as.vector(m$second) - as.vector(m$mean)^2
  expect_true(all(var_est >= -1e-8))
})

## -----------------------------------------------------------------------
## Extreme-tail stability (reviewer point 3.1): z = +-8, +-12, and very
## narrow tail intervals. Left/right censoring were already computed on
## the log scale throughout (exp() of a single BOUNDED log-difference,
## analogous to a standard inverse-Mills-ratio calculation) and were not
## actually at risk; the interval-censoring moment calculation, in
## contrast, previously computed Zinv = exp(-logZ) as a standalone
## number (overflow risk for a narrow/deep-tail interval, where logZ is
## very negative) and multiplied it by ordinary-scale dnorm(a)/dnorm(b)
## (underflow risk for extreme a/b) -- an Inf * 0 = NaN failure mode.
## These tests exercise exactly that scenario.
## -----------------------------------------------------------------------

test_that("left/right censoring moments stay finite and sensible at z = +-8, +-12", {
  for (z in c(-12, -8, 8, 12)) {
    mu <- matrix(0, 1, 1)
    sigma <- 1

    m_left <- pecme_truncated_moments(mu, sigma, "left", left = z, right = NA, y = NA)
    expect_true(all(is.finite(unlist(m_left))), info = paste("left, z =", z))
    var_left <- as.vector(m_left$second) - as.vector(m_left$mean)^2
    expect_true(var_left >= -1e-8 && var_left <= sigma^2 + 1e-6, info = paste("left var, z =", z))

    m_right <- pecme_truncated_moments(mu, sigma, "right", left = NA, right = z, y = NA)
    expect_true(all(is.finite(unlist(m_right))), info = paste("right, z =", z))
    var_right <- as.vector(m_right$second) - as.vector(m_right$mean)^2
    expect_true(var_right >= -1e-8 && var_right <= sigma^2 + 1e-6, info = paste("right var, z =", z))

    ## Deep left truncation: E[Y | Y <= z] should lie below z (mean is
    ## pulled further into the tail than the truncation point itself)
    ## for z very negative, and deep right truncation should push the
    ## mean above z for z very positive.
    if (z < 0) expect_lt(as.vector(m_left$mean), z)
    if (z > 0) expect_gt(as.vector(m_right$mean), z)
  }
})

test_that("left/right censoring moments at extreme z match a brute-force integration reference", {
  for (z in c(-12, -8, 8, 12)) {
    mu <- 0; sigma <- 1
    if (z < 0) {
      dens <- function(y) dnorm(y, mu, sigma)
      denom <- stats::integrate(dens, -Inf, z, rel.tol = 1e-12)$value
      mean_ref <- stats::integrate(function(y) y * dens(y), -Inf, z, rel.tol = 1e-12)$value / denom
      m <- pecme_truncated_moments(matrix(mu), sigma, "left", left = z, right = NA, y = NA)
      expect_equal(as.vector(m$mean), mean_ref, tolerance = 1e-6)
    } else {
      dens <- function(y) dnorm(y, mu, sigma)
      denom <- stats::integrate(dens, z, Inf, rel.tol = 1e-12)$value
      mean_ref <- stats::integrate(function(y) y * dens(y), z, Inf, rel.tol = 1e-12)$value / denom
      m <- pecme_truncated_moments(matrix(mu), sigma, "right", left = NA, right = z, y = NA)
      expect_equal(as.vector(m$mean), mean_ref, tolerance = 1e-6)
    }
  }
})

test_that("interval-censoring moments stay finite for a very narrow, deep-tail interval", {
  ## This is exactly the scenario that overflowed under the previous
  ## implementation: Zinv = exp(-logZ) is huge (the interval probability
  ## is astronomically small), multiplied by dnorm(a)/dnorm(b), which
  ## individually underflow to (near) 0 for a, b this extreme.
  cases <- list(
    c(left = -12.01, right = -12.00),
    c(left = -20.001, right = -20.000),
    c(left = 12.00, right = 12.01),
    c(left = 30.000, right = 30.001)
  )
  for (cc in cases) {
    m <- pecme_truncated_moments(matrix(0), sigma = 1, censor = "interval",
                                  left = cc[["left"]], right = cc[["right"]], y = NA)
    expect_true(all(is.finite(unlist(m))), info = paste(cc, collapse = ","))
    var_est <- as.vector(m$second) - as.vector(m$mean)^2
    expect_true(var_est >= -1e-6, info = paste(cc, collapse = ","))
    ## The truncated mean must lie within [left, right].
    expect_gte(as.vector(m$mean), cc[["left"]] - 1e-6)
    expect_lte(as.vector(m$mean), cc[["right"]] + 1e-6)
  }
})

test_that("interval-censoring log-likelihood also stays finite at the same extreme, narrow intervals", {
  cases <- list(c(-12.01, -12.00), c(-20.001, -20.000), c(12.00, 12.01), c(30.000, 30.001))
  for (cc in cases) {
    ll <- pecme_conditional_loglik(y = NA, mu = matrix(0), sigma = 1, censor = "interval",
                                    left = cc[1], right = cc[2])
    expect_true(is.finite(as.vector(ll)), info = paste(cc, collapse = ","))
  }
})

test_that("interval-censoring moments at a narrow deep-tail interval match brute-force integration", {
  left <- -12.05; right <- -12.00; mu <- 0; sigma <- 1
  dens <- function(y) dnorm(y, mu, sigma)
  denom <- stats::integrate(dens, left, right, rel.tol = 1e-12, subdivisions = 500)$value
  mean_ref <- stats::integrate(function(y) y * dens(y), left, right,
                                rel.tol = 1e-12, subdivisions = 500)$value / denom

  m <- pecme_truncated_moments(matrix(mu), sigma, "interval", left = left, right = right, y = NA)
  expect_equal(as.vector(m$mean), mean_ref, tolerance = 1e-4)
})
