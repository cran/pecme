test_that("missing_y = TRUE estimates a plausible model without dropping rows via listwise deletion", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 10, p = 3, n_active = 2,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.3,
                              missing_prob = 0.15, seed = 55)
  expect_true(anyNA(sim$data$y))

  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               missing_y = TRUE, penalty = "none", lambda = 0, quadrature = 20, max_iter = 100)

  expect_s3_class(fit, "pecme")
  expect_true(fit$n_missing_y > 0)
  expect_equal(fit$nobs, nrow(sim$data))  # rows kept, not dropped
})

test_that("a completely missing response contributes exactly 0 to the observed log-likelihood", {
  ## Documentation/behavioral check for reviewer point 3.5: a missing
  ## response is imputed (predicted) but must NOT itself inform (beta,
  ## sigma, sigma_b). Concretely: adding more missing-response rows
  ## (with valid covariates/group) to an otherwise-identical data set
  ## should leave the fitted logLik of the ORIGINAL rows' contribution
  ## unchanged, since each added "missing" row contributes log(1) = 0.
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 3, n_active = 2,
                              beta_active = 1.5, sigma = 0.5, sigma_b = 0.3, seed = 5501)

  d_extra <- sim$data
  extra_rows <- d_extra[1:5, ]
  extra_rows$y <- NA_real_
  extra_rows$censor <- "none"
  d_with_extra <- rbind(d_extra, extra_rows)

  fit_base <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
                     penalty = "none", lambda = 0, quadrature = 15, max_iter = 60)
  fit_extra <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = d_with_extra,
                      missing_y = TRUE, penalty = "none", lambda = 0,
                      quadrature = 15, max_iter = 60)

  ## logLik should match closely: the 5 added missing rows contribute 0.
  expect_equal(fit_extra$logLik, fit_base$logLik, tolerance = 0.05)
})

test_that("missing_y = FALSE errors informatively when the response has NA", {
  sim <- simulate_pecme_data(n_groups = 10, n_per_group = 5, p = 2, missing_prob = 0.2, seed = 56)
  expect_error(
    pecme(y ~ x1 + x2, random_var = "group", data = sim$data,
          missing_y = FALSE, penalty = "none", lambda = 0),
    "NA"
  )
})
test_that("missing_x handles NA covariates via cca/mean/reg without error", {

  skip_if_not_installed("lme4")

  sim <- simulate_pecme_data(
    n_groups = 20,
    n_per_group = 8,
    p = 3,
    seed = 57
  )

  d <- sim$data

  na_idx <- sample(nrow(d), 10)
  d$x1[na_idx] <- NA

  # missing_x = "error" should fail explicitly
  expect_error(
    pecme(
      y ~ x1 + x2 + x3,
      random_var = "group",
      data = d,
      missing_x = "error",
      penalty = "none",
      lambda = 0
    ),
    "missing"
  )

  # Complete-case analysis
  expect_warning(
    {
      fit_cca <- pecme(
        y ~ x1 + x2 + x3,
        random_var = "group",
        data = d,
        missing_x = "cca",
        penalty = "none",
        lambda = 0,
        quadrature = 15,
        max_iter = 60
      )
    },
    regexp = "row\\(s\\) dropped by complete-case analysis"
  )

  expect_equal(
    fit_cca$nobs,
    nrow(d) - length(na_idx)
  )

  # Mean imputation
  expect_warning(
    {
      fit_mean <- pecme(
        y ~ x1 + x2 + x3,
        random_var = "group",
        data = d,
        missing_x = "mean",
        penalty = "none",
        lambda = 0,
        quadrature = 15,
        max_iter = 60
      )
    },
    regexp = "imputes the numeric model matrix"
  )

  expect_equal(
    fit_mean$nobs,
    nrow(d)
  )

  # Regression imputation
  expect_warning(
    {
      fit_reg <- pecme(
        y ~ x1 + x2 + x3,
        random_var = "group",
        data = d,
        missing_x = "reg",
        penalty = "none",
        lambda = 0,
        quadrature = 15,
        max_iter = 60
      )
    },
    regexp = "imputes the numeric model matrix"
  )

  expect_equal(
    fit_reg$nobs,
    nrow(d)
  )
})
