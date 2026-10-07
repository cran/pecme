#' pecme: Penalized ECME Estimation for Censored Linear Mixed Models
#'
#' Fits Gaussian random-intercept linear mixed models with a left-,
#' right-, and/or interval-censored (and optionally missing) response
#' via the ECME algorithm, with simultaneous penalized variable
#' selection (Lasso, Adaptive Lasso, SCAD, MCP, Elastic Net, Ridge, or
#' no penalty). See [pecme()] for a single fit and [grid_pecme()] for
#' penalty-parameter selection.
#'
#' @section Model:
#' \deqn{Y^*_{ij} = x_{ij}'\beta + b_j + e_{ij}, \quad e_{ij}\sim
#'   N(0,\sigma^2), \quad b_j \sim N(0,\sigma_b^2)}
#' with \eqn{Y^*} observed exactly, left-censored, right-censored,
#' interval-censored, or completely missing, depending on the row.
#'
#' @section Algorithm:
#' Liu and Rubin's (1994) ECME algorithm: an E-step (Gauss-Hermite
#' quadrature for the random intercept -- a numerical approximation to
#' the continuous integral, not an exact evaluation of it; see
#' [pecme_quadrature_check()] -- with closed-form truncated-normal
#' moments for censored/missing rows) alternates with two
#' conditional-maximization steps -- a penalized coordinate-descent
#' CM-step for \eqn{\beta} that maximizes the expected complete-data
#' objective, and an ECME CM-step for \eqn{(\sigma,\sigma_b)} that
#' maximizes the actual observed-data marginal likelihood directly.
#'
#' @importFrom stats AIC BIC coef fitted logLik nobs predict residuals
#' @importFrom graphics plot
#' @keywords internal
"_PACKAGE"
