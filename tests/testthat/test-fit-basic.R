test_that("pecme() recovers approximately correct coefficients with no censoring, no penalty", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 40, n_per_group = 8, p = 4, n_active = 4,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.4, seed = 42)

  fit <- pecme(y ~ x1 + x2 + x3 + x4, random_var = "group",
               data = sim$data, penalty = "none", lambda = 0,
               quadrature = 20, max_iter = 100, verbose = FALSE)

  expect_s3_class(fit, "pecme")
  expect_true(fit$converged)
  err <- abs(fit$coefficients - sim$beta_true[names(fit$coefficients)])
  expect_true(all(err < 0.5))
  expect_true(fit$sigma > 0 && fit$sigma_b > 0)
})

test_that("pecme() with the default variance_step = 'ecme' converges in more than one iteration", {
  ## Regression test for a critical bug: pecme_variance_ecme() returned
  ## a field called `convergence` (optim's integer code), but the ECME
  ## loop checked `vout$converged` (a logical) -- since that field did
  ## not exist, isTRUE(NULL) was always FALSE, so every fit using the
  ## default variance_step = "ecme" broke after exactly one outer
  ## iteration with reason = "variance_optimizer_failure".
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 8, p = 3, n_active = 2,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.3, seed = 4242)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 15, max_iter = 100)

  expect_identical(fit$variance_step, "ecme")
  expect_gt(fit$iterations, 1L)
  expect_false(identical(fit$reason, "variance_optimizer_failure"))
})

test_that("pecme() handles left censoring without error and stays in a plausible range", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 8, p = 3, n_active = 2,
                              beta_active = 1.2, sigma = 0.6, sigma_b = 0.3,
                              censor_type = "left", censor_prob = 0.25, seed = 7)

  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               censor = "censor", left = "left", right = "right",
               penalty = "none", lambda = 0, quadrature = 20, max_iter = 100)

  expect_s3_class(fit, "pecme")
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(fit$sigma > 0 && fit$sigma_b > 0)
})

test_that("penalized fits (lasso) strongly shrink noise coefficients", {
  skip_if_not_installed("lme4")

  sim <- simulate_pecme_data(
    n_groups = 40,
    n_per_group = 10,
    p = 10,
    n_active = 3,
    beta_active = 2,
    sigma = 0.5,
    sigma_b = 0.4,
    seed = 11
  )

  fixed <- stats::reformulate(
    paste0("x", 1:10),
    response = "y"
  )

  fit <- pecme(
    fixed,
    random_var = "group",
    data = sim$data,
    penalty = "lasso",
    lambda = 0.5,
    quadrature = 20,
    max_iter = 100,
    refit = FALSE
  )

  beta <- fit$coefficients[-1]

  ## The three simulated active predictors should remain clearly nonzero.
  expect_gt(abs(beta[1]), 1)
  expect_gt(abs(beta[2]), 1)
  expect_gt(abs(beta[3]), 1)

  ## The inactive predictors should be strongly shrunk.
  expect_true(all(abs(beta[4:10]) < 0.1))

  ## The penalized fit should converge.
  expect_true(fit$converged)
})

test_that("the automatic post-selection relaxed refit ($refit) actually runs (not silently NULL)", {
  ## Regression test for a critical bug: pecme_run()'s internal refit
  ## call still passed response = ".refit_y" to pecme(), which no
  ## longer has a `response` parameter -- pecme() errored with "unused
  ## argument", silently caught by tryCatch(), so $refit was always
  ## NULL for every penalized fit regardless of `refit = TRUE`.
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 40, n_per_group = 10, p = 8, n_active = 3,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.4, seed = 4243)
  fixed <- stats::reformulate(paste0("x", 1:8), response = "y")
  fit <- pecme(fixed, random_var = "group", data = sim$data,
               penalty = "lasso", lambda = 0.4, quadrature = 15, max_iter = 60, refit = TRUE)

  if (length(fit$selected) > 0) {
    expect_false(is.null(fit$refit))
    expect_s3_class(fit$refit, "pecme")
    expect_true(fit$refit$se_valid)
  }
})

test_that("adaptive lasso, SCAD, MCP, elastic net, and ridge all run without error", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 6, p = 5, n_active = 2,
                              sigma = 0.6, sigma_b = 0.3, seed = 21)
  fixed <- stats::reformulate(paste0("x", 1:5), response = "y")

  for (pen in c("adaptive", "scad", "mcp", "elastic", "ridge")) {
    fit <- pecme(fixed, random_var = "group", data = sim$data,
                 penalty = pen, lambda = 0.3, quadrature = 15, max_iter = 60, refit = FALSE)
    expect_s3_class(fit, "pecme")
    expect_true(all(is.finite(fit$coefficients)))
  }
})
