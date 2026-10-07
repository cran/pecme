test_that("coordinate descent reports converged = TRUE when tolerance is met before max_iter", {
  set.seed(3)
  n <- 200; p <- 5
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  beta_true <- c(1, 2, 0, -1, 0)
  y <- as.vector(X %*% beta_true) + rnorm(n, sd = 0.1)

  fit <- pecme_coordinate_descent(X, y, beta = rep(0, p), sigma = 1, lambda = 0.05,
                                   penalty = "lasso", max_iter = 500, tol = 1e-8,
                                   has_intercept = TRUE)
  expect_true(fit$converged)
  expect_identical(fit$reason, "converged")
  expect_lt(fit$iter, 500)
})

test_that("coordinate descent conditional objective is non-increasing across sweeps", {
  set.seed(4)
  n <- 150; p <- 6
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  beta_true <- c(0.5, 1, -1, 0, 0, 2)
  y <- as.vector(X %*% beta_true) + rnorm(n, sd = 0.3)

  for (pen in c("lasso", "ridge", "elastic", "scad", "mcp", "adaptive", "none")) {
    w <- if (pen == "adaptive") rep(1, p - 1) else NULL
    fit <- pecme_coordinate_descent(X, y, beta = rep(0, p), sigma = 1, lambda = 0.1,
                                     penalty = pen, weights = w, max_iter = 200, tol = 1e-9,
                                     has_intercept = TRUE)
    h <- fit$objective_history
    if (length(h) > 1) {
      expect_true(all(diff(h) <= 1e-6), info = paste("penalty =", pen))
    }
  }
})

test_that("lasso coordinate descent recovers sparsity for zero true coefficients", {
  set.seed(5)

  n <- 300
  p <- 8

  X <- cbind(
    1,
    matrix(rnorm(n * (p - 1)), n, p - 1)
  )

  beta_true <- c(0, 3, -3, 0, 0, 0, 0, 0)

  y <- as.vector(X %*% beta_true) + rnorm(n, sd = 0.2)

  fit <- pecme_coordinate_descent(
    X = X,
    y = y,
    beta = rep(0, p),
    sigma = 1,
    lambda = 2,
    penalty = "lasso",
    max_iter = 500,
    tol = 1e-9,
    has_intercept = TRUE
  )

  # The two truly active predictors should remain clearly non-zero.
  expect_gt(abs(fit$beta[2]), 0.5)
  expect_gt(abs(fit$beta[3]), 0.5)

  # LASSO should set at least one penalized coefficient exactly to zero.
  expect_true(any(abs(fit$beta[-1]) < 1e-6))

  # The coordinate-descent algorithm should converge before max_iter.
  expect_true(fit$converged)
  expect_lt(fit$iter, 500)
})
