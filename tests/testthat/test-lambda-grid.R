test_that("at lambda_max, lasso coordinate descent yields all-zero penalized coefficients", {
  set.seed(6)

  n <- 100
  p <- 5

  X <- cbind(
    1,
    matrix(rnorm(n * (p - 1)), n, p - 1)
  )

  beta_true <- c(0.5, 1, -0.7, 0.3, 0)

  y <- as.vector(X %*% beta_true) + rnorm(n, sd = 1)

  group <- rep(1:20, each = 5)
  censor <- rep("none", n)

  sigma <- 1
  sigma_b <- 0.5

  ## Construct the same E-step pseudo-response used by
  ## pecme_lambda_max().
  beta0 <- rep(0, p)

  moments <- pecme:::pecme_e_step(
    y = y,
    X = X,
    beta = beta0,
    group = group,
    sigma = sigma,
    sigma_b = sigma_b,
    censor = censor,
    left = rep(NA_real_, n),
    right = rep(NA_real_, n),
    quadrature = 15
  )

  y_adj <- as.numeric(moments$Ey - moments$Eb)

  base_line <- mean(y_adj)

  pen_index <- 2:p

  z <- as.vector(
    crossprod(
      X[, pen_index, drop = FALSE],
      y_adj - base_line
    )
  )

  lam_max <- max(abs(z)) / sigma^2

  expect_true(is.finite(lam_max))
  expect_gt(lam_max, 0)

  ## At lambda_max, the KKT condition implies that all
  ## penalized coefficients are zero.
  fit <- pecme_coordinate_descent(
    X = X,
    y = y_adj,
    beta = rep(0, p),
    sigma = sigma,
    lambda = lam_max * 1.001,
    penalty = "lasso",
    max_iter = 300,
    tol = 1e-9,
    has_intercept = TRUE
  )

  expect_true(
    all(abs(fit$beta[-1]) < 1e-6)
  )
})



test_that("pecme_lambda_grid produces a decreasing, correctly-sized sequence", {
  g <- pecme_lambda_grid(lambda_max = 5, nlambda = 20, lambda_min_ratio = 1e-3)
  expect_length(g, 20)
  expect_equal(g[1], 5)
  expect_true(all(diff(g) < 0))
  expect_equal(g[20], 5 * 1e-3, tolerance = 1e-8)
})
