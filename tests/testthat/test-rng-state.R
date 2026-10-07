test_that("simulate_pecme_data(seed = ...) does not alter the caller's global RNG state", {
  ## Reviewer point 2.8: a bare set.seed(seed) with no restoration
  ## permanently alters .Random.seed for the rest of the session.
  set.seed(1)
  before <- .Random.seed
  invisible(simulate_pecme_data(n_groups = 5, n_per_group = 3, p = 2, seed = 999))
  after <- .Random.seed
  expect_identical(before, after)
})

test_that("simulate_pecme_data(seed = ...) does not change subsequent unrelated random draws", {
  set.seed(1)
  draws_without_call <- runif(5)

  set.seed(1)
  invisible(simulate_pecme_data(n_groups = 5, n_per_group = 3, p = 2, seed = 123))
  draws_with_call <- runif(5)

  expect_identical(draws_without_call, draws_with_call)
})

test_that("simulate_pecme_data(seed = X) is itself reproducible across repeated calls", {
  sim1 <- simulate_pecme_data(n_groups = 10, n_per_group = 4, p = 3, seed = 42)
  sim2 <- simulate_pecme_data(n_groups = 10, n_per_group = 4, p = 3, seed = 42)
  expect_identical(sim1$data, sim2$data)
})

test_that("simulate_pecme_data() without a seed still works and consumes the RNG normally", {
  set.seed(1)
  s1 <- simulate_pecme_data(n_groups = 5, n_per_group = 3, p = 2)
  s2 <- simulate_pecme_data(n_groups = 5, n_per_group = 3, p = 2)
  ## Without an explicit seed, two successive calls should (almost
  ## certainly) differ, and should NOT be seed-preserved/reset between
  ## calls (that would defeat the purpose of "no seed" = keep drawing).
  expect_false(isTRUE(all.equal(s1$data$y, s2$data$y)))
})
