## -----------------------------------------------------------------------
## pecme_bootstrap(): a group-level (cluster) nonparametric bootstrap
##
## Naive standard errors are not valid for a penalized (shrinkage)
## estimate, and the automatic post-selection "relaxed refit" ($refit)
## is a practical approximation, not valid unconditional post-selection
## inference either (see ?pecme's `refit` argument). This function
## provides a more defensible, if more expensive, resampling-based alternative: resample
## whole GROUPS with replacement (appropriate for clustered/longitudinal
## data -- resampling individual rows would break the within-group
## dependence structure the random intercept is modeling), refit at the
## SAME penalty/lambda/alpha/gamma every time (never a fresh grid
## search -- mixing lambda-selection variability into the bootstrap
## would conflate two different sources of uncertainty), and summarize
## the resulting sampling distribution of the coefficients.
## -----------------------------------------------------------------------

#' Cluster (group-level) bootstrap for a pecme fit
#'
#' Resamples whole groups with replacement, refits the identical model
#' (same formula, `penalty`, `lambda`, `alpha`, `gamma`) to each
#' resampled data set, and summarizes the resulting coefficient
#' distribution -- a valid source of standard errors and confidence
#' intervals for a penalized fit, where naive SEs do not apply.
#'
#' @param object A `"pecme"` fit.
#' @param R Number of bootstrap replicates (default `200`; increase for
#'   final reporting, e.g. `1000+`).
#' @param workers Parallel workers across replicates (default: reuse
#'   `object$workers`).
#' @param ci_level Confidence level for the percentile interval.
#' @param seed Optional integer seed (locally preserved; see
#'   [simulate_pecme_data()]).
#' @param verbose Print replicate progress.
#'
#' @details
#' Because each replicate is refit at a FIXED penalty/lambda, this
#' captures sampling variability in the coefficient estimates GIVEN
#' that regularization strength, but not uncertainty in the choice of
#' lambda itself, and not formal selection-event-conditional inference
#' for which predictors end up nonzero. A replicate in which a
#' coefficient is shrunk to exactly zero contributes 0 to that
#' coefficient's bootstrap distribution (not `NA`), so the reported
#' bootstrap SE reflects both estimation and selection variability
#' together for that coefficient.
#'
#' @return An object of class `"pecme_bootstrap"`: a list with `table`
#'   (a data frame of `estimate`, `boot_mean`, `boot_se`, and percentile
#'   `ci_lower`/`ci_upper` per coefficient), `replicates` (the raw `R x
#'   p` matrix of bootstrap coefficient estimates), `R`, `n_failed`
#'   (replicates that errored or failed to converge, excluded), and
#'   `ci_level`.
#' @export
pecme_bootstrap <- function(object, R = 200L, workers = NULL, ci_level = 0.95,
                             seed = NULL, verbose = FALSE) {
  if (!inherits(object, "pecme")) stop("`object` must be a \"pecme\" fit.", call. = FALSE)
  if (is.null(workers)) workers <- object$workers
  R <- as.integer(R)

  if (!is.null(seed)) {
    return(withr::with_seed(
      seed,
      pecme_bootstrap(
        object = object,
        R = R,
        workers = workers,
        ci_level = ci_level,
        seed = NULL,
        verbose = verbose
      )
    ))
  }

  group_ids <- sort(unique(object$group))
  J <- length(group_ids)
  p <- length(object$coefficients)
  coef_names <- names(object$coefficients)
  if (is.null(object$analysis_data)) {
    stop("The fitted object does not contain the original analysis data required for formula-preserving bootstrap refits.",
         call. = FALSE)
  }

  boot_samples <- replicate(R, sample(group_ids, J, replace = TRUE), simplify = FALSE)

  one_rep <- function(sampled_groups) {
    row_idx <- lapply(sampled_groups, function(g) which(object$group == g))
    new_group <- unlist(Map(function(idx, k) rep(k, length(idx)), row_idx, seq_along(row_idx)))
    keep <- unlist(row_idx)

    ## Reuse the original analysis data rather than reconstructing it from X.
    ## This preserves factors, transformations, interactions, splines, etc.
    d <- object$analysis_data[keep, , drop = FALSE]
    d[[object$random_var]] <- factor(new_group)

    f <- tryCatch(
      pecme(object$fixed_formula, random_var = object$random_var, data = d,
            censor = object$censor[keep], left = object$left[keep], right = object$right[keep],
            penalty = object$penalty, lambda = object$lambda, alpha = object$alpha,
            gamma = object$gamma, adaptive_gamma = object$adaptive_gamma,
            missing_y = object$missing_y, missing_x = object$missing_x,
            standardize = object$standardize,
            quadrature = object$quadrature, variance_step = object$variance_step,
            refit = FALSE, init_method = "zero", verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(f) || !isTRUE(f$converged)) return(rep(NA_real_, p))
    out <- stats::setNames(rep(0, p), coef_names)
    out[names(f$coefficients)] <- f$coefficients
    out
  }

  cl <- NULL
  if (!is.null(workers) && workers > 1L) {
    cl <- pecme_make_cluster(workers)
    on.exit(pecme_stop_cluster(cl), add = TRUE)
  }

  if (!is.null(cl)) {
    parallel::clusterExport(cl, varlist = c("object", "p", "coef_names"),
                             envir = environment())
    reps <- parallel::parLapplyLB(cl, boot_samples, one_rep)
  } else {
    reps <- vector("list", R)
    for (r in seq_len(R)) {
      if (verbose && r %% max(1L, R %/% 10L) == 0L) message(sprintf("bootstrap replicate %d/%d", r, R))
      reps[[r]] <- one_rep(boot_samples[[r]])
    }
  }

  B <- do.call(rbind, reps)
  ok <- stats::complete.cases(B)
  n_failed <- sum(!ok)
  Bok <- B[ok, , drop = FALSE]

  alpha_ci <- 1 - ci_level
  ci <- apply(Bok, 2L, stats::quantile, probs = c(alpha_ci / 2, 1 - alpha_ci / 2), na.rm = TRUE)

  table <- data.frame(
    term = coef_names,
    estimate = as.numeric(object$coefficients[coef_names]),
    boot_mean = colMeans(Bok),
    boot_se = apply(Bok, 2L, stats::sd),
    ci_lower = ci[1L, ],
    ci_upper = ci[2L, ],
    row.names = NULL
  )

  structure(list(table = table, replicates = B, R = R, n_failed = n_failed, ci_level = ci_level),
            class = "pecme_bootstrap")
}

#' @export
print.pecme_bootstrap <- function(x, digits = 4L, ...) {
  cat(sprintf("Cluster bootstrap (R = %d, %d failed/excluded, %.0f%% percentile CI)\n\n",
              x$R, x$n_failed, 100 * x$ci_level))
  tab <- x$table
  tab[-1] <- lapply(tab[-1], round, digits = digits)
  print(tab, row.names = FALSE)
  invisible(x)
}
