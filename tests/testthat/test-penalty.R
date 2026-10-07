test_that("soft_threshold behaves as expected", {
  expect_equal(soft_threshold(5, 2), 3)
  expect_equal(soft_threshold(-5, 2), -3)
  expect_equal(soft_threshold(1, 2), 0)
  expect_equal(soft_threshold(0, 2), 0)
})

test_that("penalty_value: none/lasso/ridge/elastic on a simple vector", {
  beta <- c(intercept = 0.5, x1 = 2, x2 = -3, x3 = 0)
  expect_equal(penalty_value(beta, "none", lambda = 5), 0)
  expect_equal(penalty_value(beta, "lasso", lambda = 2, has_intercept = TRUE), 2 * (2 + 3 + 0))
  expect_equal(penalty_value(beta, "ridge", lambda = 2, has_intercept = TRUE), 1 * (4 + 9 + 0))
  ev <- penalty_value(beta, "elastic", lambda = 2, alpha = 0.5, has_intercept = TRUE)
  expect_equal(ev, 2 * (0.5 * 5 + 0.5 / 2 * 13))
})

test_that("penalty_value: SCAD and MCP reduce to lasso-like behavior near zero", {
  lambda <- 1
  beta <- c(0, 0.3, -0.3)  # no intercept
  ev_scad <- penalty_value(beta, "scad", lambda = lambda, gamma = 3.7, has_intercept = FALSE)
  ev_lasso <- penalty_value(beta, "lasso", lambda = lambda, has_intercept = FALSE)
  expect_equal(ev_scad, ev_lasso)  # |beta| <= lambda region: SCAD == lasso

  ev_mcp <- penalty_value(beta, "mcp", lambda = lambda, gamma = 3, has_intercept = FALSE)
  expect_true(ev_mcp <= ev_lasso)  # MCP <= lasso everywhere
})

test_that("SCAD/MCP penalty is flat (constant) beyond gamma*lambda", {
  lambda <- 1; gamma_scad <- 3.7; gamma_mcp <- 3
  big <- 100
  p1 <- penalty_value(c(big), "scad", lambda = lambda, gamma = gamma_scad, has_intercept = FALSE)
  p2 <- penalty_value(c(big * 2), "scad", lambda = lambda, gamma = gamma_scad, has_intercept = FALSE)
  expect_equal(p1, p2)

  p1 <- penalty_value(c(big), "mcp", lambda = lambda, gamma = gamma_mcp, has_intercept = FALSE)
  p2 <- penalty_value(c(big * 2), "mcp", lambda = lambda, gamma = gamma_mcp, has_intercept = FALSE)
  expect_equal(p1, p2)
})

test_that("pecme_coordinate_update matches direct optimization for lasso/ridge/elastic/mcp/scad", {
  set.seed(1)
  vj <- 4.2; sigma <- 1.3; lambda <- 0.7; zj <- 3.1

  obj <- function(b, penalty, alpha = 0.5, gamma = 3.7, wj = 1) {
    (vj * b^2 - 2 * zj * b) / (2 * sigma^2) +
      penalty_value(c(b), penalty, lambda, alpha, gamma, weights = wj, has_intercept = FALSE)
  }

  for (pen in c("lasso", "ridge", "elastic", "mcp", "scad", "none", "adaptive")) {
    wj <- if (pen == "adaptive") 1.5 else 1
    est <- pecme_coordinate_update(zj, vj, sigma, lambda, pen, alpha = 0.4, gamma = 3.7, wj = wj)
    grid <- seq(est - 1, est + 1, length.out = 4001)
    obj_grid <- vapply(grid, obj, numeric(1), penalty = pen, alpha = 0.4, gamma = 3.7, wj = wj)
    expect_lt(obj(est, pen, 0.4, 3.7, wj), min(obj_grid) + 1e-6)
  }
})
