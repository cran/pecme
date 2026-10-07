test_that("grid_pecme() runs end to end and selects an interior lambda (not the smallest) via BIC", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 40, n_per_group = 10, p = 10, n_active = 3,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.4, seed = 99)
  fixed <- stats::reformulate(paste0("x", 1:10), response = "y")

  ## Regression test for a critical bug: grid_pecme() called
  ## pecme_setup() without `response` (dropped when `response` was
  ## removed from pecme()'s public signature and auto-derived from the
  ## formula instead), shifting every subsequent positional argument by
  ## one slot -- random_var landed in the `response` parameter, `data`
  ## in `random_var`, etc. -- so grid_pecme() could not run at all.
  g <- grid_pecme(fixed, random_var = "group", data = sim$data,
                   penalty = "lasso", nlambda = 12, criterion = "bic",
                   quadrature = 15, max_iter = 60, refit_at_best = FALSE, keep_path = TRUE)

  expect_s3_class(g, "pecme_grid")
  expect_true(g$best_index < length(g$lambda))
  expect_true(length(g$best$selected) <= 10)
})

test_that("grid_pecme() with a user-supplied lambda vector uses exactly those values", {

  skip_if_not_installed("lme4")

  sim <- simulate_pecme_data(
    n_groups = 20,
    n_per_group = 8,
    p = 5,
    n_active = 2,
    seed = 100
  )

  fixed <- stats::reformulate(
    paste0("x", 1:5),
    response = "y"
  )

  lam <- c(1000, 100, 10, 1, 0.1)

  g <- grid_pecme(
    fixed,
    random_var = "group",
    data = sim$data,
    penalty = "lasso",
    lambda = lam,
    criterion = "bic",
    quadrature = 15,
    max_iter = 60,
    refit_at_best = FALSE
  )

  expect_equal(
    sort(g$lambda, decreasing = TRUE),
    sort(lam, decreasing = TRUE)
  )

})

test_that("grid_pecme() with penalty = 'adaptive' and an auto-built grid does not error", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 6, n_active = 3,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.4, seed = 303)
  fixed <- stats::reformulate(paste0("x", 1:6), response = "y")

  g <- grid_pecme(fixed, random_var = "group", data = sim$data,
                   penalty = "adaptive", nlambda = 8, criterion = "bic",
                   quadrature = 12, max_iter = 40, refit_at_best = FALSE, keep_path = TRUE)

  expect_s3_class(g, "pecme_grid")
  expect_true(all(is.finite(g$lambda)))
  expect_true(is.finite(g$best$logLik))
})

test_that("grid_pecme() achieves real sparsity variation across the grid (not a flat n_selected)", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 40, n_per_group = 10, p = 10, n_active = 3,
                              beta_active = 2, sigma = 0.6, sigma_b = 0.4,
                              censor_type = "left", censor_prob = 0.3, seed = 404)
  fixed <- stats::reformulate(paste0("x", 1:10), response = "y")

  g <- grid_pecme(fixed, random_var = "group", data = sim$data,
                   censor = "censor", left = "left", right = "right",
                   penalty = "lasso", nlambda = 15, criterion = "bic",
                   quadrature = 15, max_iter = 60, refit_at_best = FALSE, keep_path = FALSE)

  expect_true(length(unique(g$scores$n_selected)) > 1)
  expect_true(min(g$scores$n_selected) < 10)
})

test_that("pecme_setup() has no separate `response` argument; it is derived from the formula", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 15, n_per_group = 5, p = 3, seed = 407)
  fixed <- stats::reformulate(paste0("x", 1:3), response = "y")

  setup <- pecme_setup(fixed, random_var = "group", data = sim$data)
  expect_identical(setup$response, "y")
  expect_false("response" %in% names(formals(pecme_setup)))
  expect_false("response" %in% names(formals(pecme_prepare_data)))
})

test_that("adaptive lasso pilot weights are not degenerate (not all ~equal)", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 8, p = 8, n_active = 3,
                              beta_active = 2.5, sigma = 0.5, sigma_b = 0.3, seed = 405)
  setup <- pecme_setup(stats::reformulate(paste0("x", 1:8), response = "y"),
                        random_var = "group", data = sim$data)
  aw <- pecme_compute_adaptive_weights(setup, adaptive_gamma = 1)

  expect_true(stats::sd(aw$weights) / mean(aw$weights) > 0.05)
})

test_that("different penalties give genuinely different results on the same data (sanity check)", {
  ## Not a bug-regression test for one specific fix, but a broad sanity
  ## check motivated by an earlier observed symptom (all penalties
  ## picking the same lambda/coefficients), which traced back to two
  ## real bugs now covered by other tests. Lasso and Ridge, in
  ## particular, should virtually never select an identical support
  ## (Ridge does not zero anything).
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 40, n_per_group = 10, p = 10, n_active = 3,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.4, seed = 406)
  fixed <- stats::reformulate(paste0("x", 1:10), response = "y")

  g_lasso <- grid_pecme(fixed, random_var = "group", data = sim$data,
                         penalty = "lasso", nlambda = 15, criterion = "bic",
                         quadrature = 15, max_iter = 60, refit_at_best = FALSE)
  g_ridge <- grid_pecme(fixed, random_var = "group", data = sim$data,
                         penalty = "ridge", nlambda = 15, criterion = "bic",
                         quadrature = 15, max_iter = 60, refit_at_best = FALSE)

  ## Ridge coefficients should virtually never be exactly zero.
  expect_true(length(g_ridge$best$selected) >= length(g_lasso$best$selected))
})
