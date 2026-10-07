## -----------------------------------------------------------------------
## pecme(): user-facing single-model fit
##
## Internally split into:
##   pecme_setup()  -- formula parsing, missing-data handling,
##                      standardization, lme4-based initial values.
##                      Done ONCE per data set.
##   pecme_run()    -- engine call + adaptive weights + refit + object
##                      assembly for one (penalty, lambda, alpha, gamma).
## grid_pecme() (grid-search.R) calls pecme_setup() once and pecme_run()
## once per grid point, so a 50-point lambda grid does not repeat
## formula parsing / lmer initialization 50 times.
## -----------------------------------------------------------------------

#' @keywords internal
pecme_setup <- function(fixed_formula, random_var, data,
                         censor = "none", left = NULL, right = NULL,
                         missing_y = FALSE,
                         missing_x = c("error", "cca", "mean", "reg"),
                         standardize = TRUE, quadrature = 40L,
                         init_method = c("lmer", "zero", "random")) {
  missing_x <- match.arg(missing_x)
  init_method <- match.arg(init_method)
  ## `response` has no independent role: it is always the left-hand
  ## side of `fixed_formula`, derived fresh here rather than accepted
  ## as a separate argument that callers would have to keep in sync
  ## with the formula themselves.
  response_expression <- paste(deparse(fixed_formula[[2L]]), collapse = "")
  response <- response_expression

  prep <- pecme_prepare_data(data, fixed_formula, random_var, censor, left, right)
  if (any(prep$censor %in% c("left", "right", "interval")) && !is.symbol(fixed_formula[[2L]])) {
    stop("Transformed responses are not supported with censoring in `pecme()`. ",
         "Censoring bounds are not automatically transformed; use the raw response, ",
         "or create a transformed response variable yourself and supply censoring bounds ",
         "already on that transformed scale.", call. = FALSE)
  }
  y <- prep$y; group <- prep$group; censor_vec <- prep$censor
  left_vec <- prep$left; right_vec <- prep$right; has_intercept <- prep$has_intercept

  n_missing_y <- 0L
  if (missing_y) {
    flagged <- pecme_flag_missing_response(y, left_vec, right_vec, censor_vec)
    censor_vec <- flagged$censor
    n_missing_y <- flagged$n_missing
  }

  mx <- pecme_handle_missing_x(prep$X, method = missing_x, has_intercept = has_intercept)
  X_raw <- mx$X
  keep <- mx$keep
  if (!all(keep)) {
    y <- y[keep]; group <- group[keep]; censor_vec <- censor_vec[keep]
    left_vec <- left_vec[keep]; right_vec <- right_vec[keep]
  }

  ## Complete-case filtering can remove every row belonging to a group.
  ## Rebuild integer group IDs and levels from the retained observations so
  ## empty original levels are never counted as fitted groups.
  group_chr <- prep$group_levels[group]
  group_levels <- sort(unique(as.character(group_chr)))
  group <- match(as.character(group_chr), group_levels)
  if (!length(group_levels)) stop("No groups remain after missing-X handling.", call. = FALSE)

  ## NOTE: LAPACK=TRUE would silently ignore `tol` here (base R: with
  ## LAPACK=TRUE, qr()$rank is always min(nrow, ncol), i.e. rank
  ## deficiency is never detected) -- use the default LINPACK-based QR,
  ## which does respect `tol`.
  qrX <- qr(X_raw, tol = 1e-7)
  if (qrX$rank < ncol(X_raw)) {
    keep_idx <- sort(qrX$pivot[seq_len(qrX$rank)])
    dropped <- setdiff(seq_len(ncol(X_raw)), keep_idx)
    warning("Rank-deficient fixed-effects design matrix: dropping ", length(dropped),
            " aliased column(s): ", paste(colnames(X_raw)[dropped], collapse = ", "),
            call. = FALSE)
    X_raw <- X_raw[, keep_idx, drop = FALSE]
  }
  coef_names <- colnames(X_raw)
  model_matrix_cols <- colnames(X_raw)
  p <- ncol(X_raw); n <- nrow(X_raw)
  n_informative <- sum(censor_vec != "missing" & !is.na(censor_vec))
  if (n_informative < 1L) stop("No informative response rows remain for likelihood-based fitting.", call. = FALSE)

  ## Preserve the exact analysis data used to construct the model frame.
  ## Refit utilities must not try to reverse-engineer factor levels,
  ## transformations, interactions, or splines from X alone.
  analysis_data <- data[keep, , drop = FALSE]

  if (standardize) {
    std <- pecme_standardize(X_raw, has_intercept)
  } else {
    std <- list(X = X_raw, center = rep(0, p), scale = rep(1, p))
  }
  Xs <- std$X

  init_data <- data.frame(Xs)
  colnames(init_data) <- paste0("v", seq_len(p))
  init_data$.y <- y
  init_data$.g <- group
  init_formula <- stats::reformulate(
    if (has_intercept) paste0("v", seq_len(p)[-1L]) else paste0("v", seq_len(p)),
    response = ".y", intercept = has_intercept
  )
  init <- tryCatch(
    pecme_init_random_effects(init_data, ".y", init_formula, ".g",
                               has_intercept = has_intercept, init_method = init_method),
    error = function(e) {
      warning("Initialization via `lmer` failed (", conditionMessage(e),
              "); falling back to `init_method = \"zero\"`.", call. = FALSE)
      pecme_init_random_effects(init_data, ".y", init_formula, ".g",
                                 has_intercept = has_intercept, init_method = "zero")
    }
  )
  beta0 <- init$beta
  if (is.null(beta0) || length(beta0) != p) beta0 <- rep(0, p)

  ## A random start must perturb quantities actually passed to the ECME
  ## optimizer. Randomizing an unused vector of group effects has no effect.
  if (init_method == "random") {
    beta0 <- beta0 + stats::rnorm(p, mean = 0, sd = 0.05 * max(stats::sd(y, na.rm = TRUE), 1))
    sigma0 <- max(init$sigma * exp(stats::rnorm(1L, 0, 0.05)), 1e-8)
    sigma_b0 <- max(init$sigma_b * exp(stats::rnorm(1L, 0, 0.05)), 1e-8)
  } else {
    sigma0 <- init$sigma
    sigma_b0 <- init$sigma_b
  }

  list(
    y = y, X_raw = X_raw, Xs = Xs, group = group, censor = censor_vec,
    left = left_vec, right = right_vec, has_intercept = has_intercept,
    coef_names = coef_names, model_matrix_cols = model_matrix_cols,
    p = p, n = n, center = std$center, scale = std$scale,
    beta0 = beta0, sigma0 = sigma0, sigma_b0 = sigma_b0,
    terms = prep$terms, xlevels = prep$xlevels, contrasts = prep$contrasts,
    group_levels = group_levels, fixed_formula = fixed_formula,
    n_informative = n_informative,
    response = response, response_expression = response_expression,
    random_var = random_var, analysis_data = analysis_data,
    missing_y = missing_y,
    n_missing_y = n_missing_y, missing_x = missing_x, n_missing_x = mx$n_missing,
    standardize = standardize, quadrature = quadrature, init_method = init_method
  )
}

#' @keywords internal
pecme_compute_pilot <- function(setup, ridge_lambda_ratio = 1e-3,
                                 variance_step = c("ecme", "em"),
                                 max_iter = 200L, tol = 1e-6,
                                 cd_max_iter = 1000L, cd_tol = 1e-8,
                                 workers = 1L, cluster = NULL) {
  variance_step <- match.arg(variance_step)
  y <- setup$y; Xs <- setup$Xs; group <- setup$group
  censor_vec <- setup$censor; left_vec <- setup$left; right_vec <- setup$right
  has_intercept <- setup$has_intercept; quadrature <- setup$quadrature
  beta0 <- setup$beta0; sigma0 <- setup$sigma0; sigma_b0 <- setup$sigma_b0

  ## A light ridge penalty, scaled to the natural size of vj = sum(x_j^2)
  ## for STANDARDIZED columns (~ n), not to any other penalty's
  ## zero-threshold (a different, unrelated quantity/unit). Large enough
  ## to keep the pilot numerically stable under collinearity/p close to
  ## n, small enough that it barely shrinks anything -- the point of
  ## this fit is a well-calibrated (sigma, sigma_b, beta) under the
  ## ACTUAL censored/missing-data likelihood (lmer()'s own starting
  ## values ignore censoring entirely and can be a poor guide to where
  ## the model actually converges), not variable selection.
  ridge_lambda <- ridge_lambda_ratio * setup$n

  pilot_fit <- pecme_fit_engine(
    Xs, y, group, censor_vec, left_vec, right_vec, beta0, sigma0, sigma_b0,
    penalty = "ridge", lambda = ridge_lambda, has_intercept = has_intercept,
    quadrature = quadrature, max_iter = max_iter, tol = max(tol, 1e-5),
    cd_max_iter = cd_max_iter, cd_tol = cd_tol, variance_step = variance_step,
    workers = workers, cluster = cluster
  )
  if (!isTRUE(pilot_fit$converged)) {
    stop(
      "Pilot ECME fit did not converge; lambda grid calibration was not performed. ",
      "Reason: ", pilot_fit$reason, ".",
      call. = FALSE
    )
  }
  if (!all(is.finite(c(pilot_fit$beta, pilot_fit$sigma, pilot_fit$sigma_b)))) {
    stop("Pilot ECME fit returned non-finite starting values.", call. = FALSE)
  }
  list(beta0 = pilot_fit$beta, sigma0 = pilot_fit$sigma,
       sigma_b0 = pilot_fit$sigma_b, converged = TRUE,
       reason = pilot_fit$reason,
       diagnostics = pilot_fit)
}

#' @keywords internal
pecme_compute_adaptive_weights <- function(setup, adaptive_gamma = 1,
                                            variance_step = c("ecme", "em"),
                                            max_iter = 200L, tol = 1e-6,
                                            cd_max_iter = 1000L, cd_tol = 1e-8,
                                            workers = 1L, cluster = NULL,
                                            pilot = NULL) {
  variance_step <- match.arg(variance_step)
  has_intercept <- setup$has_intercept; p <- setup$p

  if (is.null(pilot)) {
    pilot <- pecme_compute_pilot(setup, ridge_lambda_ratio = 1e-3, variance_step = variance_step,
                                  max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter,
                                  cd_tol = cd_tol, workers = workers, cluster = cluster)
  }
  pen_index <- if (has_intercept) seq_len(p)[-1L] else seq_len(p)
  eps <- 1e-4 * stats::sd(setup$y, na.rm = TRUE)
  weights <- 1 / (abs(pilot$beta0[pen_index]) + eps)^adaptive_gamma

  list(weights = weights, beta0 = pilot$beta0, sigma0 = pilot$sigma0, sigma_b0 = pilot$sigma_b0)
}

#' @keywords internal
pecme_run <- function(setup, penalty = .pecme_penalties, lambda = 0, alpha = 0.5,
                       gamma = NULL, adaptive_gamma = 1, variance_step = c("ecme", "em"),
                       refit = TRUE, max_iter = 200L, tol = 1e-6,
                       cd_max_iter = 1000L, cd_tol = 1e-8, workers = 1L, cluster = NULL,
                       verbose = FALSE, call = NULL, precomputed_adaptive = NULL) {
  penalty <- match.arg(penalty)
  variance_step <- match.arg(variance_step)
  if (is.null(gamma)) gamma <- if (penalty == "scad") 3.7 else if (penalty == "mcp") 3 else 3.7

  y <- setup$y; Xs <- setup$Xs; X_raw <- setup$X_raw; group <- setup$group
  censor_vec <- setup$censor; left_vec <- setup$left; right_vec <- setup$right
  has_intercept <- setup$has_intercept; p <- setup$p; n <- setup$n
  n_informative <- setup$n_informative
  coef_names <- setup$coef_names; quadrature <- setup$quadrature

  pecme_validate_inputs(y, X_raw, group, censor_vec, left_vec, right_vec,
                         penalty, lambda, alpha, gamma, quadrature, workers,
                         allow_missing_response = setup$missing_y)

  own_cluster <- FALSE
  if (is.null(cluster) && workers > 1L) {
    cluster <- pecme_make_cluster(workers)
    own_cluster <- TRUE
  }
  if (own_cluster) on.exit(pecme_stop_cluster(cluster), add = TRUE)

  beta0 <- setup$beta0; sigma0 <- setup$sigma0; sigma_b0 <- setup$sigma_b0

  weights <- NULL
  if (penalty == "adaptive") {
    ## Reuse a precomputed ridge-pilot fit (e.g. from grid_pecme(), which
    ## computes this ONCE for the whole lambda grid) when supplied,
    ## instead of refitting an identical ridge pilot at every lambda.
    aw <- if (!is.null(precomputed_adaptive)) {
      precomputed_adaptive
    } else {
      pecme_compute_adaptive_weights(setup, adaptive_gamma, variance_step, max_iter, tol,
                                      cd_max_iter, cd_tol, workers, cluster)
    }
    weights <- aw$weights
    beta0 <- aw$beta0; sigma0 <- aw$sigma0; sigma_b0 <- aw$sigma_b0
  }

  fit <- pecme_fit_engine(
    Xs, y, group, censor_vec, left_vec, right_vec, beta0, sigma0, sigma_b0,
    penalty = penalty, lambda = lambda, alpha = alpha, gamma = gamma, weights = weights,
    has_intercept = has_intercept, quadrature = quadrature,
    max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter, cd_tol = cd_tol,
    variance_step = variance_step, workers = workers, cluster = cluster, verbose = verbose
  )

  beta_orig <- fit$beta / setup$scale
  if (has_intercept) {
    beta_orig[1L] <- fit$beta[1L] - sum(fit$beta[-1L] * setup$center[-1L] / setup$scale[-1L])
  }
  names(beta_orig) <- coef_names

  pen_index <- if (has_intercept) seq_len(p)[-1L] else seq_len(p)
  selected <- if (length(pen_index)) pen_index[abs(beta_orig[pen_index]) > 1e-8] else integer(0)
  selected_names <- coef_names[selected]

  final_moments <- pecme_e_step_dispatch(
    y, Xs, fit$beta, group, fit$sigma, fit$sigma_b, censor_vec, left_vec, right_vec,
    quadrature = quadrature, workers = workers, cluster = cluster
  )
  random_coef <- final_moments$group_Eb
  names(random_coef) <- setup$group_levels[final_moments$group_ids]

  se <- rep(NA_real_, p); names(se) <- coef_names
  covariance <- NULL
  hessian_ok <- FALSE
  if (penalty == "none" || lambda == 0) {
    se_out <- tryCatch(
      pecme_hessian_se(y, X_raw, group, beta_orig, fit$sigma, fit$sigma_b,
                        censor_vec, left_vec, right_vec, quadrature),
      error = function(e) NULL
    )
    if (!is.null(se_out)) {
      se <- se_out$se
      covariance <- se_out$covariance
      hessian_ok <- isTRUE(se_out$valid)
    }
  }

  refit_out <- NULL
  if (isTRUE(refit) && penalty != "none" && lambda > 0 && length(selected)) {
    refit_cols <- if (has_intercept) c(1L, selected) else selected
    Xr <- X_raw[, refit_cols, drop = FALSE]
    predictor_idx <- if (has_intercept) seq_len(ncol(Xr))[-1L] else seq_len(ncol(Xr))
    predictor_names <- colnames(Xr)[predictor_idx]
    safe_names <- paste0("v", seq_along(predictor_idx))
    refit_name_map <- stats::setNames(predictor_names, safe_names)

    d_refit <- data.frame(Xr[, predictor_idx, drop = FALSE], check.names = FALSE)
    colnames(d_refit) <- safe_names
    d_refit$.refit_y <- y
    d_refit$.refit_g <- setup$group_levels[group]

    refit_out <- tryCatch(
      pecme(
        fixed_formula = stats::reformulate(safe_names, response = ".refit_y", intercept = has_intercept),
        random_var = ".refit_g", data = d_refit,
        censor = censor_vec, left = left_vec, right = right_vec,
        penalty = "none", lambda = 0, missing_y = FALSE, missing_x = "error",
        standardize = setup$standardize, quadrature = quadrature, workers = workers,
        variance_step = variance_step, refit = FALSE, init_method = setup$init_method,
        max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter, cd_tol = cd_tol
      ),
      error = function(e) { warning("Post-selection refit failed: ", conditionMessage(e), call. = FALSE); NULL }
    )
    if (!is.null(refit_out)) {
      names(refit_out$coefficients) <-
        c(if (has_intercept) "(Intercept)", predictor_names)
      
      names(refit_out$se) <- names(refit_out$coefficients)
      if (!is.null(refit_out$covariance)) {
        dimnames(refit_out$covariance) <- list(names(refit_out$coefficients),
                                                 names(refit_out$coefficients))
      }
      
      refit_out$refit_name_map <- refit_name_map
      ## Metadata used by predict.pecme() so that users can call
      ## predict(fit$refit, newdata = original_data) without needing
      ## the synthetic v1/v2/... columns used internally by the refit.
      refit_out$refit_is_relaxed <- TRUE
      refit_out$refit_parent_formula <- setup$fixed_formula
      refit_out$refit_parent_terms <- setup$terms
      refit_out$refit_parent_xlevels <- setup$xlevels
      refit_out$refit_parent_contrasts <- setup$contrasts
      refit_out$refit_parent_model_matrix_cols <- setup$coef_names
      refit_out$refit_parent_random_var <- setup$random_var
    }
  }

  df_current <- pecme_df(beta_orig, has_intercept)
  n_selected_current <- length(selected)
  ic <- pecme_ic(fit$logLik, df_current, n_informative, length(pen_index),
                 gamma_ebic = 0.5, n_selected = n_selected_current)

  structure(list(
    call = call, coefficients = beta_orig, se = se, covariance = covariance, se_valid = hessian_ok,
    sigma = fit$sigma, sigma_b = fit$sigma_b, random_coef = random_coef,
    logLik = fit$logLik, penalty_value = fit$objective - (-fit$logLik),
    objective = fit$objective, AIC = ic$AIC, BIC = ic$BIC, EBIC = ic$EBIC, GCV = ic$GCV,
    df = pecme_df(beta_orig, has_intercept), nobs = n, n_informative = n_informative,
    n_groups = length(setup$group_levels),
    iterations = fit$iterations, converged = fit$converged, reason = fit$reason,
    outer_converged = fit$outer_converged,
    inner_converged = fit$inner_converged,
    inner_failed = fit$inner_failed,
    inner_failure_count = fit$inner_failure_count,
    inner_converged_history = fit$inner_converged_history,
    inner_iterations = fit$inner_iterations,
    variance_converged = fit$variance_converged,
    variance_failed = fit$variance_failed,
    variance_failure_count = fit$variance_failure_count,
    variance_converged_history = fit$variance_converged_history,
    variance_iterations = fit$variance_iterations,
    monotonicity_violations = fit$monotonicity_violations,
    monotonicity_violation_history = fit$monotonicity_violation_history,
    monotonicity_checked = fit$monotonicity_checked,
    monotonicity_tol = fit$monotonicity_tol,
    logLik_history = fit$logLik_history,
    objective_history = fit$objective_history,
    penalty = penalty, lambda = lambda, alpha = alpha, gamma = gamma,
    adaptive_weights = weights, adaptive_gamma = adaptive_gamma,
    selected = selected, selected_names = selected_names,
    has_intercept = has_intercept, standardize = setup$standardize,
    center = setup$center, scale = setup$scale, coef_names = coef_names,
    terms = setup$terms, xlevels = setup$xlevels, contrasts = setup$contrasts,
    fixed_formula = setup$fixed_formula, response = setup$response, random_var = setup$random_var,
    group_levels = setup$group_levels, variance_step = variance_step,
    quadrature = quadrature, workers = workers, missing_y = setup$missing_y,
    n_missing_y = setup$n_missing_y, missing_x = setup$missing_x, n_missing_x = setup$n_missing_x,
    fitted_group = final_moments, refit = refit_out,
    y = y, X = X_raw, group = group, censor = censor_vec, left = left_vec, right = right_vec,
    analysis_data = setup$analysis_data,
    response_expression = setup$response_expression
  ), class = "pecme")
}

#' Fit a penalized ECME censored/missing-response linear mixed model
#'
#' Fits \eqn{Y^*_{ij} = x_{ij}'\beta + b_j + e_{ij}}, \eqn{e_{ij} \sim
#' N(0,\sigma^2)}, \eqn{b_j \sim N(0,\sigma_b^2)}, where the response
#' \eqn{Y^*} may be left-, right-, or interval-censored (and/or
#' missing), by the ECME algorithm of Liu and Rubin (1994). The random
#' intercept is integrated out numerically by fixed-order Gauss-Hermite
#' quadrature (a numerical approximation to the continuous integral,
#' not an exact evaluation of it; see [pecme_quadrature_check()] for a
#' sensitivity diagnostic and `quadrature` below).
#' \eqn{\beta} is estimated by a penalized CM-step (coordinate descent);
#' \eqn{(\sigma,\sigma_b)} are estimated, by default, with a true ECME
#' CM-step that directly maximizes the observed-data marginal
#' likelihood (not merely its EM lower bound), following the strategy
#' of Vaida and Liu (2009) for censored linear mixed models.
#'
#' @param fixed_formula A formula for the fixed effects, e.g. `y ~ x1 + x2`.
#' @param random_var Character: name of the grouping (cluster) column in
#'   `data`. Currently a random INTERCEPT only (no random slopes).
#' @param data A data frame.
#' @param censor Either `"none"`/`"left"`/`"right"`/`"interval"` applied
#'   to every row, a column name in `data` holding per-row labels from
#'   that set, or a vector of the same length as `nrow(data)`. Rows do
#'   not need to be pre-labelled `"missing"`: set `missing_y = TRUE`
#'   and leave the response `NA` for those rows.
#' @param left,right Column names or numeric vectors of censoring bounds
#'   (used only for the corresponding `censor` labels). With censoring,
#'   the response must be supplied on its analysis scale (transformed
#'   responses such as `log(y)` are rejected); censoring bounds are never
#'   transformed automatically.
#' @param penalty One of `"none"`, `"lasso"`, `"adaptive"`, `"scad"`,
#'   `"mcp"`, `"elastic"`, `"ridge"`.
#' @param lambda Non-negative penalty strength.
#' @param alpha Elastic Net mixing parameter in `[0,1]` (`1` = Lasso,
#'   `0` = Ridge); ignored unless `penalty = "elastic"`.
#' @param gamma Extra shape parameter for SCAD (`> 2`, default `3.7`) or
#'   MCP (`> 1`, default `3`); ignored otherwise.
#' @param adaptive_gamma Exponent used to build Adaptive Lasso weights
#'   from a Ridge pilot fit: `weights = 1 / |beta_pilot|^adaptive_gamma`.
#' @param missing_y If `TRUE`, rows with a missing response (or a
#'   missing bound for a labelled-censored row) are treated as
#'   completely missing at random and IMPUTED via the E-step's
#'   posterior mean at every iteration, instead of raising an error.
#'   Important: such a row contributes exactly `log(1) = 0` to the
#'   observed-data log-likelihood -- it carries NO direct information
#'   about `(beta, sigma, sigma_b)`; the imputed value is a model-
#'   implied *prediction* for that row given the parameters, not
#'   evidence that helped estimate them. The parameters themselves are
#'   identified entirely by the observed and censored rows.
#' @param missing_x How to handle `NA` in the fixed-effects predictors:
#'   `"error"` (default), `"cca"`, `"mean"`, or `"reg"`. The `"mean"`
#'   and `"reg"` options impute the numeric model matrix and are intended
#'   only for suitable numeric predictors. For factors, interactions,
#'   splines, or nonlinear terms, impute raw variables before constructing
#'   the model matrix (e.g. with `mice`). See [pecme_handle_missing_x()]
#'   and [pecme_mi_combine()].
#' @param standardize Standardize predictors before penalization
#'   (recommended; default `TRUE`).
#' @param quadrature Number of Gauss-Hermite nodes for the random
#'   intercept (default `40`). This is a fixed-order numerical
#'   approximation to the random-intercept integral, not an exact
#'   evaluation of it; use [pecme_quadrature_check()] on a fitted model
#'   to check sensitivity to this choice, and increase it if `sigma_b`
#'   is large relative to `sigma`, or convergence looks unstable.
#' @param workers Number of parallel worker processes for the E-step
#'   (default `1`). Set to e.g. `32` on a 32-core machine.
#' @param variance_step `"ecme"` (default) or `"em"`; see Details.
#' @param refit If `TRUE` (default) and `penalty != "none"`, an
#'   additional unpenalized refit restricted to the selected
#'   coefficients ("relaxed refit") is stored in `$refit`. This reduces
#'   the shrinkage bias of the penalized estimates for reporting, but
#'   it is NOT generally valid unconditional post-selection inference:
#'   the selection step itself is not accounted for, so treat its `se`
#'   and p-values as approximate/exploratory, not confirmatory.
#' @param init_method `"lmer"` (default), `"zero"`, or `"random"`. With
#'   `"random"`, starting values are drawn from the caller's ongoing
#'   global random-number stream (unlike [simulate_pecme_data()]'s
#'   `seed`, this is NOT locally seed-preserved -- by design, since the
#'   point of `"random"` is to give a genuinely different starting
#'   point on repeated calls, e.g. for a multi-start convergence check).
#'   Call `set.seed()` yourself beforehand if you need the specific
#'   starting values to be reproducible.
#' @param max_iter,tol Outer ECME iteration limit/tolerance.
#' @param cd_max_iter,cd_tol Inner coordinate-descent iteration limit/tolerance.
#' @param verbose Print iteration progress.
#'
#' @return An object of class `"pecme"`.
#' @seealso [grid_pecme()], [pecme_mi_combine()]
#' @references
#' Liu, C. and Rubin, D. B. (1994). The ECME algorithm. *Biometrika*, 81(4), 633-648.
#'
#' Vaida, F. and Liu, L. (2009). Fast implementation for normal mixed
#' effects models with censored response. *JCGS*, 18(1), 151-162.
#'
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave
#' penalized likelihood and its oracle properties. *JASA*, 96(456), 1348-1360.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under
#' minimax concave penalty. *Annals of Statistics*, 38(2), 894-942.
#' @export
pecme <- function(fixed_formula, random_var, data,
                   censor = "none", left = NULL, right = NULL,
                   penalty = c("none", "lasso", "adaptive", "scad", "mcp", "elastic", "ridge"),
                   lambda = 0, alpha = 0.5, gamma = NULL, adaptive_gamma = 1,
                   missing_y = FALSE,
                   missing_x = c("error", "cca", "mean", "reg"),
                   standardize = TRUE, quadrature = 40L, workers = 1L,
                   variance_step = c("ecme", "em"), refit = TRUE,
                   init_method = c("lmer", "zero", "random"),
                   max_iter = 200L, tol = 1e-6, cd_max_iter = 1000L, cd_tol = 1e-8,
                   verbose = FALSE) {
  penalty <- match.arg(penalty)
  call <- match.call()

  setup <- pecme_setup(fixed_formula, random_var, data, censor, left, right,
                        missing_y, missing_x, standardize, quadrature, init_method)

  pecme_run(setup, penalty, lambda, alpha, gamma, adaptive_gamma, variance_step, refit,
            max_iter, tol, cd_max_iter, cd_tol, workers, cluster = NULL,
            verbose = verbose, call = call)
}

#' @keywords internal
pecme_hessian_se <- function(y, X, group, beta, sigma, sigma_b, censor, left, right, quadrature) {
  fn <- function(b) pecme_marginal_loglik(y, X, group, b, sigma, sigma_b, censor, left, right, quadrature)
  H <- stats::optimHess(beta, fn)
  V <- tryCatch(solve(-H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V)) || any(diag(V) < 0)) return(NULL)
  V <- (V + t(V)) / 2
  eig <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                  error = function(e) NA_real_)
  if (any(!is.finite(eig)) || min(eig) <= 0) return(NULL)
  se <- sqrt(pmax(diag(V), 0))
  names(se) <- names(beta)
  dimnames(V) <- list(names(beta), names(beta))
  list(se = se, covariance = V, valid = all(is.finite(se)))
}
