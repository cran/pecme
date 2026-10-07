## -----------------------------------------------------------------------
## pecme_metrics(): the "give me everything for my thesis" function
##
## MSE/RMSE/MAE/MAPE/SMAPE and two flavors of R^2 are computed either
## in-sample (default) or out-of-sample (pass `newdata`/`newresponse`).
## AIC/BIC/EBIC/GCV/logLik are read off the fitted training-data model
## (they are not meaningful to recompute on new data). LOO is an
## optional, genuinely refit-based leave-one-group-out cross-validated
## log-likelihood -- NOT the PSIS-LOO of the `loo` package, which needs
## posterior draws pecme does not produce (see Details).
## -----------------------------------------------------------------------

#' Comprehensive fit-quality metrics for a pecme model
#'
#' Collects the quantities typically reported for a mixed-model fit:
#' pointwise prediction error (MSE/RMSE/MAE/MAPE/SMAPE), two kinds of
#' pseudo-R-squared, information criteria, and (optionally) a
#' cross-validated predictive log-likelihood.
#'
#' @param object A `"pecme"` fit (from [pecme()] or `grid_pecme()$best`).
#' @param newdata,newresponse Optional: evaluate prediction-error metrics
#'   out-of-sample instead of on the training data. `newdata` needs the
#'   predictor and grouping columns. The observed response is looked up
#'   automatically from `newdata` using the same column name as in the
#'   fitted model's formula (`response` has no independent role -- see
#'   [pecme()]); pass `newresponse` explicitly only if `newdata` does
#'   not itself contain that column (e.g. truly held-out rows with
#'   unknown outcomes, being scored on other rows' observed values).
#'   Rows with a missing/`NA` response (however obtained) are silently
#'   dropped from the error metrics, since there is no single observed
#'   value to compare a point prediction against.
#' @param loo If `TRUE`, additionally compute a cross-validated
#'   predictive log-likelihood by actually refitting the model with one
#'   group (`loo_method = "exact"`) or one fold of groups
#'   (`loo_method = "kfold"`) held out at a time, and evaluating the
#'   held-out marginal log-likelihood. This is a genuine, if expensive,
#'   out-of-sample criterion; see Details for why it is not the PSIS-LOO
#'   of the `loo` package.
#' @param loo_method `"exact"` (leave-one-group-out; `J` refits) or
#'   `"kfold"` (leave-one-fold-out; `loo_folds` refits).
#' @param loo_folds Number of folds when `loo_method = "kfold"`.
#' @param workers Parallel workers for the LOO refits (default: reuse
#'   `object$workers`).
#' @param ebic_gamma Extended-BIC tuning constant (see [grid_pecme()]);
#'   only affects the recomputed `EBIC` entry.
#'
#' @details
#' # Why there is no WAIC
#' WAIC (and DIC) are defined from the variance of the pointwise
#' log-likelihood *across posterior draws* of the parameters. `pecme`
#' is a frequentist penalized-likelihood estimator: it returns a single
#' point estimate (plus, for an unpenalized fit, an asymptotic Hessian-
#' based covariance), not a posterior sample. Reporting a WAIC number
#' here would therefore not be a real WAIC -- rather than approximate
#' it, `$WAIC` is returned as `NA` with an explanatory `attr`. For
#' out-of-sample model comparison, use `AIC`/`BIC`/`EBIC` (fast, from
#' the fitted likelihood) or `loo = TRUE` (slower, but a genuine
#' held-out predictive log-likelihood, which is the same quantity WAIC
#' and PSIS-LOO both aim to approximate).
#'
#' # Two R-squared values
#' `R2_marginal` and `R2_conditional` follow Nakagawa & Schielzeth
#' (2013): `R2_marginal = Var(fixed) / (Var(fixed) + sigma_b^2 +
#' sigma^2)` uses only the fixed-effects predictor's variance in the
#' numerator; `R2_conditional` additionally credits the random
#' intercept, `(Var(fixed) + sigma_b^2) / (Var(fixed) + sigma_b^2 +
#' sigma^2)`. These are the R-squared values normally reported for a
#' mixed model. `R2_predictive` = `1 - MSE / Var(y)` is the more
#' familiar OLS-style quantity computed directly from predictions; it
#' is provided for intuition/comparability but is not the standard
#' mixed-model R-squared. `R2_predictive_adj` heuristically applies the
#' usual OLS adjusted-R2 formula with `object$df`; this adjustment is
#' not a rigorously derived quantity for mixed models and should be
#' reported as a diagnostic, not a citation-grade statistic.
#'
#' @return An object of class `"pecme_metrics"` (a list; has a `print`
#'   method).
#' @export
pecme_metrics <- function(object, newdata = NULL, newresponse = NULL,
                           loo = FALSE, loo_method = c("exact", "kfold"),
                           loo_folds = 10L, workers = NULL, ebic_gamma = 0.5) {
  if (!inherits(object, "pecme")) stop("`object` must be a \"pecme\" fit.", call. = FALSE)
  loo_method <- match.arg(loo_method)
  if (is.null(workers)) workers <- object$workers

  ## --- prediction-error metrics (in-sample or out-of-sample) -------------
  if (is.null(newdata)) {
    y_true <- object$y
    y_pred <- fitted.pecme(object)
    ## Only rows with a genuinely observed response are comparable to a
    ## point prediction; censored/missing rows contribute no error term.
    ok <- object$censor == "none" & !is.na(y_true)
  } else {
    y_pred <- predict.pecme(object, newdata = newdata, level = "group")
    if (is.null(newresponse)) {
      ## `response` has no independent role: it is always the left-hand
      ## side of object$fixed_formula. If `newdata` already carries that
      ## column (the common case -- e.g. a held-out data frame built the
      ## same way as the training data), use it directly rather than
      ## requiring the caller to re-supply values that are already there.
      ## Evaluate the LHS through the original formula so transformed
      ## responses (e.g. log(y), sqrt(y), I(y^2)) are scored on the same
      ## response scale on which the model was fitted.
      mf_new <- tryCatch(
        stats::model.frame(object$fixed_formula, data = newdata,
                           xlev = object$xlevels, na.action = stats::na.pass),
        error = function(e) e
      )
      if (inherits(mf_new, "error")) {
        stop(
          "The response variable is not present in `newdata`. ",
          "Please supply it using `newresponse`.",
          call. = FALSE
        )
      }
      y_true <- as.numeric(stats::model.response(mf_new))
    } else {
      y_true <- as.numeric(newresponse)
    }
    ok <- !is.na(y_true)
  }
  pm <- pecme_prediction_metrics(y_true[ok], y_pred[ok])

  ## --- R-squared: predictive (OLS-style) ---------------------------------
  r2_pred <- pm$R2
  r2_pred_adj <- pecme_adjusted_r2(r2_pred, pm$n, object$df)

  ## --- R-squared: Nakagawa & Schielzeth (2013) marginal/conditional -----
  eta_fixed <- as.vector(object$X %*% object$coefficients)
  var_fixed <- stats::var(eta_fixed)
  var_b <- object$sigma_b^2
  var_e <- object$sigma^2
  var_total <- var_fixed + var_b + var_e
  r2_marginal <- if (var_total > 0) var_fixed / var_total else NA_real_
  r2_conditional <- if (var_total > 0) (var_fixed + var_b) / var_total else NA_real_

  out <- list(
    n = pm$n, df = object$df, nobs = object$nobs, n_groups = object$n_groups,
    logLik = object$logLik, AIC = object$AIC, BIC = object$BIC,
    EBIC = if (isTRUE(all.equal(ebic_gamma, 0.5))) object$EBIC else {
      n_sel <- length(object$selected)
      p_tot <- length(if (object$has_intercept) object$coef_names[-1L] else object$coef_names)
      lc <- if (is.finite(n_sel) && n_sel >= 0L && n_sel <= p_tot) {
        lchoose(p_tot, n_sel)
      } else NA_real_
      if (is.finite(lc)) object$BIC + 2 * ebic_gamma * lc else NA_real_
    },
    GCV = object$GCV,
    MSE = pm$MSE, RMSE = pm$RMSE, MAE = pm$MAE, MAPE = pm$MAPE, SMAPE = pm$SMAPE,
    R2_predictive = r2_pred, R2_predictive_adj = r2_pred_adj,
    R2_marginal = r2_marginal, R2_conditional = r2_conditional,
    sigma = object$sigma, sigma_b = object$sigma_b,
    LOO = NA_real_, LOO_SE = NA_real_, LOOIC = NA_real_,
    WAIC = structure(NA_real_, note = paste(
      "Not defined for a frequentist penalized-likelihood estimator",
      "(no posterior draws); see ?pecme_metrics Details. Use AIC/BIC/EBIC",
      "or loo = TRUE instead."
    )),
    out_of_sample = !is.null(newdata)
  )

  if (isTRUE(loo)) {
    loo_res <- pecme_loo(object, method = loo_method, folds = loo_folds, workers = workers)
    out$LOO <- loo_res$elpd_loo
    out$LOO_SE <- loo_res$se
    out$LOOIC <- if (is.finite(loo_res$elpd_loo)) -2 * loo_res$elpd_loo else NA_real_
    out$LOO_successful_folds <- loo_res$n_successful_folds
    out$loo_detail <- loo_res
  }

  structure(out, class = "pecme_metrics")
}

#' @export
print.pecme_metrics <- function(x, digits = 4L, ...) {
  cat("pecme fit metrics", if (x$out_of_sample) "(out-of-sample)" else "(in-sample)", "\n\n")
  cat("Error metrics (n =", x$n, "):\n")
  err <- c(MSE = x$MSE, RMSE = x$RMSE, MAE = x$MAE, MAPE = x$MAPE, SMAPE = x$SMAPE)
  print(round(err, digits))

  cat("\nR-squared:\n")
  r2 <- c(R2_predictive = x$R2_predictive, R2_predictive_adj = x$R2_predictive_adj,
          R2_marginal = x$R2_marginal, R2_conditional = x$R2_conditional)
  print(round(r2, digits))
  cat("(marginal/conditional follow Nakagawa & Schielzeth 2013; see ?pecme_metrics)\n")

  cat("\nInformation criteria (df =", x$df, "):\n")
  ic <- c(logLik = x$logLik, AIC = x$AIC, BIC = x$BIC, EBIC = x$EBIC, GCV = x$GCV)
  print(round(ic, digits))

  if (!is.na(x$LOO)) {
    cat(sprintf("\nLeave-one-group-out CV: elpd_loo = %.*f (SE %.*f), LOOIC = %.*f\n",
                digits, x$LOO, digits, x$LOO_SE, digits, x$LOOIC))
  } else {
    cat("\nLOO: not computed (set loo = TRUE in pecme_metrics(); refits the model once per group/fold).\n")
  }
  cat("WAIC: not applicable to this estimator (see ?pecme_metrics Details).\n")
  invisible(x)
}

#' Leave-one-group-out (or leave-one-fold-out) cross-validated log-likelihood
#'
#' Refits the model with each group (or fold of groups) held out in
#' turn and evaluates the held-out marginal log-likelihood at the
#' refit parameters -- a genuine, frequentist out-of-sample predictive
#' criterion for clustered data (leave-one-*observation*-out is not
#' meaningful here since observations within a group share a random
#' intercept).
#'
#' @inheritParams pecme_metrics
#' @param method Character string specifying the resampling method.
#' @param folds Integer specifying the number of folds for K-fold cross-validation.
#' @return A list with `elpd_loo` (summed held-out log-likelihood, on
#'   the same scale as `logLik`), `se` (a naive standard error across
#'   successful group/fold contributions), `n_successful_folds`, and
#'   `pointwise` (a data frame of one
#'   row per group/fold).
#' @export
pecme_loo <- function(object, method = c("exact", "kfold"), folds = 10L, workers = NULL) {
  if (!inherits(object, "pecme")) stop("`object` must be a \"pecme\" fit.", call. = FALSE)
  method <- match.arg(method)
  if (is.null(workers)) workers <- object$workers

  group_ids <- sort(unique(object$group))
  J <- length(group_ids)
  fold_of <- if (method == "exact") seq_len(J) else {
    k <- max(2L, min(as.integer(folds), J))
    sample(rep(seq_len(k), length.out = J))
  }
  n_folds <- length(unique(fold_of))

  refit_and_score <- function(f) {
    train_groups <- group_ids[fold_of != f]
    test_groups  <- group_ids[fold_of == f]
    tr <- object$group %in% train_groups
    te <- object$group %in% test_groups
    if (!any(tr) || !any(te)) return(NA_real_)

    Xtr_raw <- object$X[tr, , drop = FALSE]
    if (isTRUE(object$standardize)) {
      std <- pecme_standardize(Xtr_raw, object$has_intercept)
    } else {
      std <- list(X = Xtr_raw, center = rep(0, ncol(Xtr_raw)), scale = rep(1, ncol(Xtr_raw)))
    }
    Xs_tr <- std$X
    group_tr <- match(object$group[tr], sort(unique(object$group[tr])))

    beta0 <- object$coefficients * std$scale
    if (object$has_intercept) {
      beta0[1L] <- object$coefficients[1L] +
        sum(object$coefficients[-1L] * std$center[-1L])
    }
    
    setup_f <- object

    setup_f$X_raw <- Xtr_raw
    setup_f$Xs <- Xs_tr
    setup_f$y <- object$y[tr]
    setup_f$group <- group_tr
    setup_f$censor <- object$censor[tr]
    setup_f$left <- object$left[tr]
    setup_f$right <- object$right[tr]

    setup_f$group_levels <-
      object$group_levels[sort(unique(object$group[tr]))]

    setup_f$n <- sum(tr)
    setup_f$p <- ncol(Xtr_raw)

    setup_f$n_informative <-
      sum(
        setup_f$censor != "missing" &
          !is.na(setup_f$censor)
      )

    setup_f$center <- std$center
    setup_f$scale <- std$scale

    ## Warm starts expressed on the standardized training-fold scale.
    setup_f$beta0 <- beta0

    ## These fields are required by pecme_run(), but are not stored
    ## under the setup names in a fitted "pecme" object.
    setup_f$sigma0 <- object$sigma
    setup_f$sigma_b0 <- object$sigma_b

    setup_f$coef_names <- colnames(Xtr_raw)

    adaptive_f <- NULL
    if (object$penalty == "adaptive") {
      adaptive_f <- tryCatch(
        pecme_compute_adaptive_weights(
          setup_f, object$adaptive_gamma, object$variance_step,
          max_iter = 200L, tol = 1e-6, cd_max_iter = 1000L, cd_tol = 1e-8,
          workers = 1L, cluster = NULL, pilot = NULL
        ), error = function(e) NULL
      )
      if (is.null(adaptive_f)) return(NA_real_)
    }

    fit_f <- tryCatch(
      pecme_run(
        setup_f, penalty = object$penalty, lambda = object$lambda,
        alpha = object$alpha, gamma = object$gamma,
        adaptive_gamma = object$adaptive_gamma, variance_step = object$variance_step,
        refit = FALSE, workers = 1L, cluster = NULL,
        precomputed_adaptive = adaptive_f
      ), error = function(e) NULL
    )
    if (is.null(fit_f)) return(NA_real_)

    pars_ok <- 
      !is.null(fit_f$coefficients) &&
      all(is.finite(fit_f$coefficients)) &&
      is.finite(fit_f$sigma) &&
      fit_f$sigma > 0 &&
      is.finite(fit_f$sigma_b) &&
      fit_f$sigma_b >= 0

    if (!pars_ok) return(NA_real_)

    Xte <- object$X[te, , drop = FALSE]
    group_te <- match(object$group[te], sort(unique(object$group[te])))
    
    score_f <- tryCatch(
      pecme_marginal_loglik(
        object$y[te],
        Xte,
        group_te,
        fit_f$coefficients,
        fit_f$sigma,
        fit_f$sigma_b,
        object$censor[te],
        object$left[te],
        object$right[te],
        object$quadrature
      ),
      error = function(e) NA_real_
    )

    if (length(score_f) != 1L || !is.finite(score_f)) {
      return(NA_real_)
    }

    score_f
  }

  cl <- NULL
  if (!is.null(workers) && workers > 1L) {
    cl <- pecme_make_cluster(workers)
    on.exit(pecme_stop_cluster(cl), add = TRUE)
  }
  if (!is.null(cl)) {
    parallel::clusterExport(cl, varlist = c("object"), envir = environment())
    scores <- unlist(parallel::parLapplyLB(cl, as.list(seq_len(n_folds)), refit_and_score))
  } else {
    scores <- vapply(seq_len(n_folds), refit_and_score, numeric(1))
  }

  scores_ok <- scores[is.finite(scores)]
  n_successful <- length(scores_ok)
  elpd <- if (n_successful) sum(scores_ok) else NA_real_
  se <- if (n_successful > 1L) stats::sd(scores_ok) * sqrt(n_successful) else NA_real_

  list(elpd_loo = elpd, se = se,
       pointwise = data.frame(fold = seq_len(n_folds), logLik = scores),
       method = method, n_folds = n_folds, n_successful_folds = n_successful)
}
