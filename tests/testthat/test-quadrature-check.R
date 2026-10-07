test_that("pecme_quadrature_check() (method = 'loglik') runs and includes the reference quadrature", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 3, n_active = 2,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.4, seed = 501)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 30, max_iter = 60)

  qc <- pecme_quadrature_check(fit, quadrature_values = c(10, 20, 40))
  expect_s3_class(qc, "pecme_quadrature_check")
  expect_true(30 %in% qc$table$quadrature)
  expect_true(all(is.finite(qc$table$logLik)))

  ## At the model's own converged parameters, logLik should already be
  ## close to its large-quadrature limit -- differences across a
  ## reasonable range of node counts should be small for a
  ## well-behaved, moderate sigma_b/sigma ratio like this one.
  expect_lt(diff(range(qc$table$logLik)), 1)
})

test_that("pecme_quadrature_check() (method = 'refit') runs and reports parameter changes", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 2, n_active = 1,
                              beta_active = 1.2, sigma = 0.5, sigma_b = 0.3, seed = 502)
  fit <- pecme(y ~ x1 + x2, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 20, max_iter = 40)

  qc <- pecme_quadrature_check(fit, quadrature_values = c(10, 30), method = "refit")
  expect_s3_class(qc, "pecme_quadrature_check")
  expect_true(all(c("sigma", "sigma_b", "max_abs_beta_change") %in% names(qc$table)))
})

test_that("print.pecme_quadrature_check() does not error", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 2, seed = 503)
  fit <- pecme(y ~ x1 + x2, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 20, max_iter = 40)
  qc <- pecme_quadrature_check(fit, quadrature_values = c(15, 25))
  expect_output(print(qc), "quadrature")
})
