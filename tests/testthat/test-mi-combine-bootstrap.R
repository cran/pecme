test_that("pecme_mi_combine() works for unpenalized fits (the always-valid case)", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, n_active = 2, seed = 601)
  d1 <- sim$data; d2 <- sim$data  # trivial "imputations" (identical data) for a mechanical check

  f1 <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = d1,
              penalty = "none", lambda = 0, quadrature = 15, max_iter = 40)
  f2 <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = d2,
              penalty = "none", lambda = 0, quadrature = 15, max_iter = 40)

  pooled <- pecme_mi_combine(list(f1, f2))
  expect_true(all(is.finite(pooled$se)))
  expect_true(all(is.finite(pooled$estimate)))
  ## Identical inputs: pooled estimate should equal the (identical) per-fit estimate.
  expect_equal(pooled$estimate, as.numeric(f1$coefficients[pooled$term]), tolerance = 1e-6)
})

test_that("pecme_mi_combine() errors clearly (not silently NA) for penalized fits by default", {
  ## Reviewer point 1.4: this used to silently return NA within-
  ## imputation variance and NA pooled SE for penalized fits.
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 6, n_active = 2,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.3, seed = 602)
  fixed <- stats::reformulate(paste0("x", 1:6), response = "y")

  f1 <- pecme(fixed, random_var = "group", data = sim$data,
              penalty = "lasso", lambda = 0.3, quadrature = 15, max_iter = 40, refit = FALSE)
  f2 <- pecme(fixed, random_var = "group", data = sim$data,
              penalty = "lasso", lambda = 0.3, quadrature = 15, max_iter = 40, refit = FALSE)

  expect_false(f1$se_valid)
  expect_error(pecme_mi_combine(list(f1, f2)), "valid standard errors")
})

test_that("pecme_mi_combine(use_refit = TRUE) pools relaxed refits, including across differing selected sets", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 30, n_per_group = 8, p = 6, n_active = 2,
                              beta_active = 2.5, sigma = 0.4, sigma_b = 0.3, seed = 603)
  fixed <- stats::reformulate(paste0("x", 1:6), response = "y")

  f1 <- pecme(fixed, random_var = "group", data = sim$data,
              penalty = "lasso", lambda = 0.4, quadrature = 15, max_iter = 40, refit = TRUE)
  f2 <- pecme(fixed, random_var = "group", data = sim$data,
              penalty = "lasso", lambda = 0.35, quadrature = 15, max_iter = 40, refit = TRUE)

  skip_if(is.null(f1$refit) || is.null(f2$refit), "no predictors selected in this configuration")

  pooled <- pecme_mi_combine(list(f1, f2), use_refit = TRUE)
  expect_true(all(is.finite(pooled$se)))
  ## The term set should be the union of both fits' selected predictors
  ## plus the intercept.
  all_terms <- union(names(f1$refit$coefficients), names(f2$refit$coefficients))
  expect_setequal(pooled$term, all_terms)
})

test_that("pecme_mi_combine(use_refit = TRUE) errors clearly if any fit lacks a $refit", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, seed = 604)
  f_none <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
                   penalty = "none", lambda = 0, quadrature = 12, max_iter = 30)
  expect_error(pecme_mi_combine(list(f_none, f_none), use_refit = TRUE), "\\$refit")
})

test_that("pecme_bootstrap() runs and produces sensible SE/CI", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 25, n_per_group = 8, p = 4, n_active = 2,
                              beta_active = 2, sigma = 0.5, sigma_b = 0.3, seed = 605)
  fit <- pecme(y ~ x1 + x2 + x3 + x4, random_var = "group", data = sim$data,
               penalty = "lasso", lambda = 0.3, quadrature = 12, max_iter = 30, refit = FALSE)

  bs <- pecme_bootstrap(fit, R = 15, seed = 1)
  expect_s3_class(bs, "pecme_bootstrap")
  expect_equal(nrow(bs$table), length(fit$coefficients))
  expect_true(all(bs$table$ci_lower <= bs$table$ci_upper))
  expect_true(all(is.finite(bs$table$boot_se)))
})

test_that("pecme_bootstrap(seed = ...) is reproducible and preserves the caller's RNG state", {
  skip_if_not_installed("lme4")
  sim <- simulate_pecme_data(n_groups = 20, n_per_group = 6, p = 3, seed = 606)
  fit <- pecme(y ~ x1 + x2 + x3, random_var = "group", data = sim$data,
               penalty = "none", lambda = 0, quadrature = 12, max_iter = 30)

  set.seed(5)
  before <- .Random.seed
  bs1 <- pecme_bootstrap(fit, R = 10, seed = 77)
  after <- .Random.seed
  expect_identical(before, after)

  bs2 <- pecme_bootstrap(fit, R = 10, seed = 77)
  expect_equal(bs1$table$boot_mean, bs2$table$boot_mean)
})
