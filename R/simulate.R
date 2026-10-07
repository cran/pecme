## -----------------------------------------------------------------------
## simulate_pecme_data(): synthetic censored/missing random-intercept
## LMM data, for examples, tests, and Monte Carlo checks of the
## estimator (bias, sparsistency, etc). Ships as a function rather than
## a bundled .rda so every example is self-contained and reproducible
## from a single `seed` argument.
## -----------------------------------------------------------------------

#' Simulate data from a censored/missing random-intercept linear mixed model
#'
#' Generates \eqn{Y^*_{ij} = \beta_0 + x_{ij}'\beta + b_j + e_{ij}},
#' \eqn{e_{ij}\sim N(0,\sigma^2)}, \eqn{b_j\sim N(0,\sigma_b^2)}, with
#' `p` standard-normal (optionally correlated) predictors, of which
#' only `n_active` have a non-zero coefficient -- convenient for
#' checking variable-selection performance. Censoring and/or missingness
#' can be injected on top of the latent `y_star`.
#'
#' @param n_groups Number of clusters/groups.
#' @param n_per_group Observations per group (recycled if a vector).
#' @param p Number of predictors (in addition to the intercept).
#' @param n_active Number of predictors with a non-zero true
#'   coefficient (the rest are exactly zero -- useful for checking
#'   variable selection).
#' @param beta_active Non-zero coefficient magnitude used for the
#'   `n_active` active predictors (sign alternates).
#' @param rho Pairwise predictor correlation (compound symmetry), `0`
#'   for independent predictors.
#' @param sigma,sigma_b Residual and random-intercept standard deviations.
#' @param censor_type `"none"`, `"left"`, `"right"`, or `"interval"`
#'   applied to a `censor_prob` fraction of observations (chosen
#'   independently at random); the rest are `"none"`.
#' @param censor_prob Fraction of observations to censor.
#' @param missing_prob Fraction of (non-censored) observations to set
#'   completely missing (`NA`) in the response.
#' @param seed Optional integer seed for reproducibility. The caller's
#'   global RNG state is saved and restored on exit, so passing `seed`
#'   has no side effect on random draws made elsewhere in the calling
#'   session (unlike a bare `set.seed()` call).
#'
#' @return A list with `data` (a data frame with columns `y`, `x1..xp`,
#'   `group`, `censor`, `left`, `right`), `beta_true` (length `p + 1`,
#'   intercept first), `sigma_true`, `sigma_b_true`, and `y_star` (the
#'   uncensored, non-missing latent response, for validation).
#' @examples
#' sim <- simulate_pecme_data(n_groups = 20, n_per_group = 8, p = 6,
#'                             n_active = 3, censor_type = "left",
#'                             censor_prob = 0.2, seed = 1)
#' head(sim$data)
#' @export
simulate_pecme_data <- function(n_groups = 30L, n_per_group = 10L, p = 10L,
                                 n_active = 4L, beta_active = 1.5, rho = 0,
                                 sigma = 1, sigma_b = 0.7,
                                 censor_type = c("none", "left", "right", "interval"),
                                 censor_prob = 0, missing_prob = 0, seed = NULL) {
  censor_type <- match.arg(censor_type)

  if (!is.null(seed)) {
    return(withr::with_seed(
      seed,
      simulate_pecme_data(
        n_groups = n_groups,
        n_per_group = n_per_group,
        p = p,
        n_active = n_active,
        beta_active = beta_active,
        rho = rho,
        sigma = sigma,
        sigma_b = sigma_b,
        censor_type = censor_type,
        censor_prob = censor_prob,
        missing_prob = missing_prob,
        seed = NULL
      )
    ))
  }

  n_per_group <- rep_len(as.integer(n_per_group), n_groups)
  group <- rep(seq_len(n_groups), times = n_per_group)
  n <- length(group)

  if (rho == 0) {
    X <- matrix(stats::rnorm(n * p), n, p)
  } else {
    Sigma <- matrix(rho, p, p); diag(Sigma) <- 1
    R <- chol(Sigma)
    X <- matrix(stats::rnorm(n * p), n, p) %*% R
  }
  colnames(X) <- paste0("x", seq_len(p))

  beta <- rep(0, p)
  n_active <- min(n_active, p)
  if (n_active > 0) {
    active_idx <- seq_len(n_active)
    beta[active_idx] <- beta_active * (-1)^(seq_len(n_active) + 1)
  }
  beta_full <- c(intercept = 1, beta)
  names(beta_full) <- c("(Intercept)", colnames(X))

  b <- stats::rnorm(n_groups, 0, sigma_b)
  e <- stats::rnorm(n, 0, sigma)
  y_star <- as.vector(beta_full[1] + X %*% beta + b[group] + e)

  y <- y_star
  left <- rep(NA_real_, n); right <- rep(NA_real_, n)
  censor <- rep("none", n)

  if (!is.finite(censor_prob) || censor_prob < 0 || censor_prob > 1) {
    stop("`censor_prob` must be finite and in [0, 1].", call. = FALSE)
  }

  if (censor_type != "none" && censor_prob > 0) {
    if (censor_type == "left") {
      thr <- as.numeric(stats::quantile(y_star, probs = censor_prob,
                                        names = FALSE, type = 7))
      eligible <- which(y_star <= thr)
      censor[eligible] <- "left"
      left[eligible] <- thr
      y[eligible] <- NA_real_
    } else if (censor_type == "right") {
      thr <- as.numeric(stats::quantile(y_star, probs = 1 - censor_prob,
                                        names = FALSE, type = 7))
      eligible <- which(y_star >= thr)
      censor[eligible] <- "right"
      right[eligible] <- thr
      y[eligible] <- NA_real_
    } else {
      n_cens <- min(n, max(1L, round(censor_prob * n)))
      idx <- sample.int(n, size = n_cens, replace = FALSE)
      width <- stats::sd(y_star, na.rm = TRUE) * 0.5
      delta <- stats::runif(length(idx), -width/2, width/2)
      centre <- y_star[idx] + delta
      left[idx] <- centre - width/2
      right[idx] <- centre + width/2
      censor[idx] <- "interval"
      y[idx] <- NA_real_
    }
  }

  if (missing_prob > 0) {
    eligible <- which(censor == "none")
    idx <- sample(eligible, size = round(missing_prob * length(eligible)))
    y[idx] <- NA
  }

  data <- data.frame(X)
  data$y <- y
  data$group <- factor(group)
  data$censor <- censor
  data$left <- left
  data$right <- right

  list(data = data, beta_true = beta_full, sigma_true = sigma, sigma_b_true = sigma_b,
       y_star = y_star)
}
