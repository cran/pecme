test_that("pecme_metrics() returns all the thesis-standard quantities and a properly NA/labelled WAIC", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 4, n_active = 2,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.4, seed = 200)
  fit <- pecme(y ~ x1 + x2 + x3 + x4, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 15, max_iter = 60)

  m <- pecme_metrics(fit)
  expect_s3_class(m, "pecme_metrics")

  for (nm in c("MSE", "RMSE", "MAE", "MAPE", "SMAPE",
               "R2_predictive", "R2_predictive_adj", "R2_marginal", "R2_conditional",
               "logLik", "AIC", "BIC", "EBIC", "GCV")) {
    expect_true(is.finite(m[[nm]]), info = nm)
  }

  expect_true(is.na(m$WAIC))
  expect_true(!is.null(attr(m$WAIC, "note")))

  ## marginal R2 <= conditional R2 (random effect can only add explained variance)
  expect_lte(m$R2_marginal, m$R2_conditional + 1e-8)
  expect_true(m$R2_marginal >= 0 && m$R2_conditional <= 1)

  ## RMSE^2 should equal MSE
  expect_equal(m$RMSE^2, m$MSE, tolerance = 1e-8)
})

test_that("pecme_metrics(loo = TRUE, loo_method = 'kfold') runs and returns finite LOO/LOOIC", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, n_active = 2,
                              beta_active = 1.2, sigma = 0.6, sigma_b = 0.3, seed = 201)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 12, max_iter = 40)

  m <- pecme_metrics(fit, loo = TRUE, loo_method = "kfold", loo_folds = 4)
  expect_true(is.finite(m$LOO))
  expect_equal(m$LOOIC, -2 * m$LOO, tolerance = 1e-8)
})

test_that("pecme_metrics() auto-extracts the response from newdata when the column is present", {
  ## Regression test: newdata almost always already contains the
  ## response column (e.g. a held-out data frame built the same way as
  ## the training data), so requiring a redundant, separately-supplied
  ## `newresponse` in that case was unnecessary friction inconsistent
  ## with the rest of the package (response has no independent role;
  ## see ?pecme). This should now just work without `newresponse`.
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, seed = 202)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 12, max_iter = 40)

  m_auto <- pecme_metrics(fit, newdata = sim$data)
  m_explicit <- pecme_metrics(fit, newdata = sim$data)

  expect_true(m_auto$out_of_sample)
  expect_true(is.finite(m_auto$MSE))
  expect_equal(m_auto$MSE, m_explicit$MSE)
})

test_that("pecme_metrics() errors clearly when newdata lacks the response column and newresponse is not supplied", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, seed = 202)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 12, max_iter = 40)

  newdata_no_y <- sim$data[, setdiff(names(sim$data), "y")]
  expect_error(pecme_metrics(fit, newdata = newdata_no_y), "newresponse")

  m <- pecme_metrics(fit, newdata = newdata_no_y, newresponse = sim$data$y)
  expect_true(m$out_of_sample)
  expect_true(is.finite(m$MSE))
})

test_that("pecme_effective_df(): Ridge uses the closed-form trace formula, others use the naive count or NA", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 8, p = 6, n_active = 3,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.3, seed = 203)
  fixed <- stats::reformulate(paste0("x", 1:6), response = "y")

  fit_ridge <- pecme(fixed, random_var = "group", data = sim$data,
                      penalty = "ridge", lambda = 0.5, quadrature = 12, max_iter = 40, refit = FALSE)
  df_eff <- pecme_effective_df(fit_ridge)
  expect_true(is.finite(df_eff))
  ## Effective df for a genuinely shrinking ridge fit should be less
  ## than the naive (near-maximal, since ridge rarely zeroes anything)
  ## nonzero-coefficient count.
  expect_lt(df_eff, pecme_df(fit_ridge$coefficients, fit_ridge$has_intercept))

  fit_scad <- pecme(fixed, random_var = "group", data = sim$data,
                     penalty = "scad", lambda = 0.5, quadrature = 12, max_iter = 40, refit = FALSE)
  df_scad <- pecme_effective_df(fit_scad)
  expect_true(is.na(df_scad))
  expect_true(!is.null(attr(df_scad, "note")))

  fit_none <- pecme(fixed, random_var = "group", data = sim$data,
                     penalty = "none", lambda = 0, quadrature = 12, max_iter = 40)
  expect_equal(pecme_effective_df(fit_none), pecme_df(fit_none$coefficients, fit_none$has_intercept))
})
