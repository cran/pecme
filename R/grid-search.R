## -----------------------------------------------------------------------
## grid_pecme(): lambda grid search
##
## Replaces the original grid_ECME.R / grid_ECME_parallel.R pair (two
## ~1100-1400 line, near-duplicate files that had already started to
## diverge -- e.g. a broken lambda_max in the parallel copy) with a
## single implementation. Data preparation and lme4 initialization
## (pecme_setup()) run ONCE regardless of grid size; only the ECME fit
## itself (pecme_run()) is repeated per lambda.
##
## Model selection is by information criterion (BIC/EBIC/AIC/GCV) or,
## optionally, leave-groups-out K-fold cross-validated predictive
## log-likelihood -- NEVER by raw observed log-likelihood, which is
## monotone in model complexity and would always select the smallest
## lambda (i.e. no regularization at all).
## -----------------------------------------------------------------------

#' Grid search over the penalty parameter lambda
#'
#' @inheritParams pecme
#' @param lambda Either `NULL` (default: build a data-dependent grid
#'   via [pecme_lambda_max()]/[pecme_lambda_grid()]) or a user-supplied
#'   numeric vector of candidate lambda values.
#' @param nlambda Number of grid points when `lambda = NULL`.
#' @param lambda_min_ratio `lambda_min = lambda_min_ratio * lambda_max`
#'   when `lambda = NULL`; defaults to the usual `n > p` / `n <= p`
#'   heuristic (see [pecme_lambda_grid()]).
#' @param criterion One of `"bic"` (default), `"ebic"`, `"aic"`,
#'   `"gcv"`, or `"cv"`. `"cv"` performs leave-groups-out K-fold
#'   cross-validation (refits the whole model `cv_folds` times per
#'   lambda) and selects by held-out predictive log-likelihood; this is
#'   far more expensive than the closed-form criteria.
#' @param ebic_gamma Extended-BIC tuning constant in `[0,1]` (Chen and
#'   Chen, 2008); `0` reduces EBIC to ordinary BIC. Only used when
#'   `criterion = "ebic"`.
#' @param cv_folds Number of folds (by group, not by observation) for
#'   `criterion = "cv"`.
#' @param cv_seed Optional integer seed for the fold assignment.
#' @param parallel_over `"lambda"` (default: one full ECME fit per
#'   worker, embarrassingly parallel -- usually the better choice for
#'   `nlambda > workers`) or `"groups"` (parallelize the E-step across
#'   groups within each, sequential, lambda -- better for a single/small
#'   grid on data with many groups). Only one level of parallelism is
#'   ever used at a time; `workers` are never oversubscribed by nesting
#'   both.
#' @param refit_at_best As in [pecme()]'s `refit`, but applied only to
#'   the final selected lambda (grid points themselves are fit with
#'   `refit = FALSE` to avoid `nlambda` extra unpenalized re-fits).
#' @param keep_path Keep the fitted `beta` at every grid point (for
#'   plotting a regularization path)? Default `TRUE`.
#'
#' @return An object of class `"pecme_grid"`: a list with `$best` (the
#'   selected `"pecme"` fit, refit per `refit_at_best`), `$lambda`,
#'   `$scores` (a data frame of one row per lambda with every
#'   criterion, `df`, `logLik`, and convergence status), `$best_index`,
#'   and `$criterion`.
#' @export
grid_pecme <- function(fixed_formula, random_var, data,
                        censor = "none", left = NULL, right = NULL,
                        penalty = c("lasso", "adaptive", "scad", "mcp", "elastic", "ridge", "none"),
                        lambda = NULL, nlambda = 50L, lambda_min_ratio = NULL,
                        alpha = 0.5, gamma = NULL, adaptive_gamma = 1,
                        criterion = c("bic", "ebic", "aic", "gcv", "cv"),
                        ebic_gamma = 0.5, cv_folds = 5L, cv_seed = NULL,
                        missing_y = FALSE, missing_x = c("error", "cca", "mean", "reg"),
                        standardize = TRUE, quadrature = 40L, workers = 1L,
                        parallel_over = c("lambda", "groups"),
                        variance_step = c("ecme", "em"), refit_at_best = TRUE,
                        init_method = c("lmer", "zero", "random"),
                        max_iter = 200L, tol = 1e-6, cd_max_iter = 1000L, cd_tol = 1e-8,
                        keep_path = TRUE, verbose = FALSE) {
  penalty <- match.arg(penalty)
  criterion <- match.arg(criterion)
  parallel_over <- match.arg(parallel_over)
  variance_step <- match.arg(variance_step)
  missing_x <- match.arg(missing_x)
  init_method <- match.arg(init_method)
  if (is.null(gamma)) gamma <- if (penalty == "scad") 3.7 else if (penalty == "mcp") 3 else 3.7
  call <- match.call()

  setup <- pecme_setup(fixed_formula, random_var, data, censor, left, right,
                        missing_y, missing_x, standardize, quadrature, init_method)

  ## Compute a single, lightly-ridge-penalized pilot fit under the
  ## ACTUAL censored/missing-data likelihood -- not lmer()'s raw,
  ## censoring-blind initial values -- and use its converged (sigma,
  ## sigma_b) to calibrate lambda_max. Using the uncalibrated initial
  ## sigma0/sigma_b0 here can leave lambda_max too small relative to
  ## where the model actually converges, in which case NO lambda in the
  ## grid reaches true sparsity for ANY penalty, and every information
  ## criterion ends up selecting essentially the same, least-shrunk
  ## point -- which looks like (and would otherwise be mistaken for)
  ## different penalties "all agreeing" by coincidence.
  pilot_workers <- if (workers > 1L && parallel_over == "groups") workers else 1L
  pilot <- pecme_compute_pilot(setup, ridge_lambda_ratio = 1e-3, variance_step = variance_step,
                                max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter,
                                cd_tol = cd_tol, workers = pilot_workers, cluster = NULL)
  setup$beta0 <- pilot$beta0; setup$sigma0 <- pilot$sigma0; setup$sigma_b0 <- pilot$sigma_b0

  ## Adaptive Lasso's weights are a simple transform of this SAME pilot
  ## fit's beta -- no separate/duplicate pilot fit needed.
  adaptive_pilot <- NULL
  if (penalty == "adaptive") {
    adaptive_pilot <- pecme_compute_adaptive_weights(
      setup, adaptive_gamma, variance_step, max_iter, tol, cd_max_iter, cd_tol,
      workers = pilot_workers, cluster = NULL, pilot = pilot
    )
  }

  if (is.null(lambda)) {
    lam_max <- pecme_lambda_max(setup$Xs, setup$y, setup$group, setup$sigma0, setup$sigma_b0,
                                 setup$censor, setup$left, setup$right,
                                 penalty = if (penalty == "none") "lasso" else penalty,
                                 alpha = alpha, weights = if (penalty == "adaptive") adaptive_pilot$weights else NULL,
                                 has_intercept = setup$has_intercept, quadrature = quadrature)
    lambda <- pecme_lambda_grid(lam_max, nlambda = nlambda, lambda_min_ratio = lambda_min_ratio,
                                 n = setup$n, p = setup$p)
  }
  lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  nl <- length(lambda)

  cl <- NULL
  if (workers > 1L) {
    cl <- pecme_make_cluster(workers)
    on.exit(pecme_stop_cluster(cl), add = TRUE)
  }

  if (!is.null(cl) && parallel_over == "lambda") {
    ## Embarrassingly parallel over lambda: each worker runs one
    ## complete, sequential ECME fit (workers = 1 inside the task).
    parallel::clusterExport(cl, varlist = c("setup", "penalty", "alpha", "gamma", "adaptive_gamma",
                                             "variance_step", "max_iter", "tol", "cd_max_iter", "cd_tol",
                                             "adaptive_pilot"),
                             envir = environment())
    fits <- parallel::parLapplyLB(cl, as.list(lambda), function(lam) {
      tryCatch(
        pecme_run(setup, penalty = penalty, lambda = lam, alpha = alpha, gamma = gamma,
                  adaptive_gamma = adaptive_gamma, variance_step = variance_step, refit = FALSE,
                  max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter, cd_tol = cd_tol,
                  workers = 1L, cluster = NULL, verbose = FALSE, precomputed_adaptive = adaptive_pilot),
        error = function(e) structure(list(error = conditionMessage(e), lambda = lam),
                                      class = "pecme_grid_error")
      )
    })
  } else {
    ## Sequential over lambda; if parallel_over == "groups" the SAME
    ## cluster (created once above) is reused for the E-step at every
    ## lambda instead of being torn down and rebuilt each iteration.
    fits <- vector("list", nl)
    for (i in seq_len(nl)) {
      if (verbose) message(sprintf("[%d/%d] lambda = %.6g", i, nl, lambda[i]))
      fits[[i]] <- tryCatch(
        pecme_run(setup, penalty = penalty, lambda = lambda[i], alpha = alpha, gamma = gamma,
                  adaptive_gamma = adaptive_gamma, variance_step = variance_step, refit = FALSE,
                  max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter, cd_tol = cd_tol,
                  workers = if (is.null(cl)) 1L else workers, cluster = cl, verbose = FALSE,
                  precomputed_adaptive = adaptive_pilot),
        error = function(e) structure(list(error = conditionMessage(e), lambda = lambda[i]),
                                      class = "pecme_grid_error")
      )
    }
  }

  get_num <- function(f, nm) {
    if (inherits(f, "pecme_grid_error") || is.null(f) || is.null(f[[nm]])) NA_real_
    else as.numeric(f[[nm]])
  }
  get_log <- function(f) {
    if (inherits(f, "pecme_grid_error") || is.null(f) || is.null(f$logLik)) NA_real_
    else as.numeric(f$logLik)
  }
  get_bool <- function(f) {
    if (inherits(f, "pecme_grid_error") || is.null(f) || is.null(f$converged)) FALSE
    else isTRUE(f$converged)
  }
  scores <- data.frame(
    lambda = lambda,
    df = vapply(fits, get_num, numeric(1), nm = "df"),
    n_selected = vapply(fits, function(f)
      if (inherits(f, "pecme_grid_error") || is.null(f) || is.null(f$selected)) NA_integer_
      else as.integer(length(f$selected)), integer(1)),
    logLik = vapply(fits, get_log, numeric(1)),
    AIC = vapply(fits, get_num, numeric(1), nm = "AIC"),
    BIC = vapply(fits, get_num, numeric(1), nm = "BIC"),
    EBIC = vapply(fits, get_num, numeric(1), nm = "EBIC"),
    GCV = vapply(fits, get_num, numeric(1), nm = "GCV"),
    converged = vapply(fits, get_bool, logical(1))
  )

  ## Recompute EBIC only when EBIC is the requested selection criterion.
  ## Failed grid fits can have n_selected = NA; nonfinite/invalid k must
  ## never enter the combinatorial term.
  p_total <- length(if (setup$has_intercept) seq_len(setup$p)[-1L] else seq_len(setup$p))
  if (criterion == "ebic") {
    n_sel <- scores$n_selected
    log_choose <- vapply(n_sel, function(k) {
      if (!is.finite(k) || k < 0 || k > p_total) return(NA_real_)
      lchoose(p_total, k)
    }, numeric(1))
    scores$EBIC <- scores$BIC + 2 * ebic_gamma * log_choose
  }

  if (penalty %in% c("lasso", "adaptive", "scad", "mcp", "elastic") &&
      length(unique(scores$n_selected)) == 1L) {
    warning(
      "Every lambda on the grid selected the same ", scores$n_selected[1L],
      " of ", p_total, " predictor(s): no sparsity variation was achieved ",
      "anywhere on this grid, so different penalties will tend to agree ",
      "with each other by default (nothing to disagree about) and the ",
      "selection criterion cannot meaningfully choose among them. This ",
      "usually means lambda_max is too small relative to where the model ",
      "converges. Try a larger explicit `lambda` vector (e.g. multiply ",
      "the current grid's `lambda` by 10-100), or re-run with a larger ",
      "`nlambda`/smaller `lambda_min_ratio` and inspect `$scores$n_selected` ",
      "directly.", call. = FALSE
    )
  }

  if (criterion == "cv") {
    cv_ll <- pecme_cv_score(setup, lambda, penalty, alpha, gamma, adaptive_gamma, variance_step,
                             max_iter, tol, cd_max_iter, cd_tol, cv_folds, cv_seed, verbose = verbose)
    scores$CV_logLik <- cv_ll
    score_col <- -scores$CV_logLik
  } else {
    score_col <- scores[[toupper(criterion)]]
  }

  valid <- is.finite(score_col) & scores$converged
  if (!any(valid)) {
    stop("No valid finite grid candidate was obtained.", call. = FALSE)
  }
  best_index <- which(valid)[which.min(score_col[valid])]

  best <- pecme_run(setup, penalty = penalty, lambda = lambda[best_index], alpha = alpha, gamma = gamma,
                     adaptive_gamma = adaptive_gamma, variance_step = variance_step,
                     refit = refit_at_best, max_iter = max_iter, tol = tol,
                     cd_max_iter = cd_max_iter, cd_tol = cd_tol,
                     workers = if (is.null(cl)) 1L else workers,
                     cluster = cl, verbose = verbose, call = call, precomputed_adaptive = adaptive_pilot)

  structure(list(
    best = best, best_index = best_index, lambda = lambda, criterion = criterion,
    scores = scores, penalty = penalty, alpha = alpha, gamma = gamma,
    path = if (keep_path) {
      do.call(rbind, lapply(fits, function(f)
        if (inherits(f, "pecme_grid_error") || is.null(f)) rep(NA_real_, setup$p)
        else as.numeric(f$coefficients)))
    } else NULL,
    coef_names = setup$coef_names, call = call
  ), class = "pecme_grid")
}

#' @keywords internal
pecme_cv_score <- function(setup, lambda, penalty, alpha, gamma, adaptive_gamma, variance_step,
                            max_iter, tol, cd_max_iter, cd_tol, cv_folds, cv_seed, verbose = FALSE) {
  group_ids <- sort(unique(setup$group))
  J <- length(group_ids)
  cv_folds <- max(2L, min(as.integer(cv_folds), J))

  fold_of_group <- if (is.null(cv_seed)) {
    sample(rep(seq_len(cv_folds), length.out = J))
  } else {
    withr::with_seed(
      cv_seed,
      sample(rep(seq_len(cv_folds), length.out = J))
    )
  }

  scores <- matrix(NA_real_, nrow = length(lambda), ncol = cv_folds)

  for (k in seq_len(cv_folds)) {
    train_groups <- group_ids[fold_of_group != k]
    test_groups  <- group_ids[fold_of_group == k]
    tr <- setup$group %in% train_groups
    te <- setup$group %in% test_groups
    if (!any(tr) || !any(te)) next

    setup_k <- setup
    setup_k$y <- setup$y[tr]
    ## Xs is a row subset of the FULL-data standardized design matrix
    ## (i.e. standardized using all groups' mean/sd, not train-fold-only
    ## statistics). This is a minor, deliberate simplification: it keeps
    ## lambda on the same scale across folds and against the final
    ## full-data fit, at the cost of a small amount of test-fold
    ## information leaking into the training fold's column scaling
    ## (not into the response or the fitted coefficients themselves).
    setup_k$Xs <- setup$Xs[tr, , drop = FALSE]
    setup_k$X_raw <- setup$X_raw[tr, , drop = FALSE]
    setup_k$group <- match(setup$group[tr], sort(unique(setup$group[tr])))
    setup_k$censor <- setup$censor[tr]; setup_k$left <- setup$left[tr]; setup_k$right <- setup$right[tr]
    setup_k$group_levels <- setup$group_levels[sort(unique(setup$group[tr]))]
    setup_k$n <- sum(tr)

    for (i in seq_along(lambda)) {
      if (verbose) message(sprintf("[cv fold %d/%d] lambda %d/%d", k, cv_folds, i, length(lambda)))
      fit_k <- tryCatch(
        pecme_run(setup_k, penalty = penalty, lambda = lambda[i], alpha = alpha, gamma = gamma,
                  adaptive_gamma = adaptive_gamma, variance_step = variance_step, refit = FALSE,
                  max_iter = max_iter, tol = tol, cd_max_iter = cd_max_iter, cd_tol = cd_tol,
                  workers = 1L, cluster = NULL, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit_k)) next

      X_te <- setup$X_raw[te, , drop = FALSE]
      group_te <- match(setup$group[te], sort(unique(setup$group[te])))
      ll <- tryCatch(
        pecme_marginal_loglik(setup$y[te], X_te, group_te, fit_k$coefficients,
                               fit_k$sigma, fit_k$sigma_b, setup$censor[te], setup$left[te], setup$right[te],
                               setup$quadrature),
        error = function(e) NA_real_
      )
      scores[i, k] <- ll
    }
  }

  out <- rep(NA_real_, length(lambda))
  for (i in seq_along(lambda)) {
    ok <- is.finite(scores[i, ])
    if (any(ok)) out[i] <- sum(scores[i, ok])
  }
  out
}
