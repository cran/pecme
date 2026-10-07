## -----------------------------------------------------------------------
## Missing-data handling
##
## Response `y` missing at random: reduces to the "missing" censoring
## category in truncation.R / estep.R -- no extra theory needed, the
## E-step already returns the EM-optimal (posterior mean) IMPUTATION
## E[Y* | X, group, model] for these rows every iteration. This is
## prediction/imputation, not estimation information: a "missing" row
## contributes exactly log(1) = 0 to the observed-data log-likelihood
## (see truncation.R), so it does not directly inform (beta, sigma,
## sigma_b) -- those are identified entirely by the observed and
## censored rows. The imputed value only reflects, and moves with,
## whatever the current parameter estimates already say.
##
## Covariate `X` missing: NOT covered by the censored-response ECME
## likelihood (that would require a joint model for X, a materially
## different and larger theoretical problem). `pecme()` offers
## practical, clearly-labelled options via `missing_x`; for rigorous
## propagation of imputation uncertainty into standard errors, combine
## multiple `mice`-imputed data sets with [pecme_mi_combine()].
## -----------------------------------------------------------------------

#' Recode a response/censoring specification to include missing values
#'
#' Given raw vectors `y`, `left`, `right` and a `censor` label vector
#' (`"none"`, `"left"`, `"right"`, `"interval"`), returns a new `censor`
#' vector in which any row with `is.na(y)` (for `"none"`/`"left"`/
#' `"right"`/`"interval"` rows that have no usable censoring bound
#' either) is relabelled `"missing"`.
#'
#' @keywords internal
pecme_flag_missing_response <- function(y, left, right, censor) {
  censor <- as.character(censor)
  n <- length(censor)
  if (length(y) != n) stop("y/censor length mismatch.", call. = FALSE)
  if (length(left)  != n) left  <- rep(NA_real_, n)
  if (length(right) != n) right <- rep(NA_real_, n)

  is_missing <-
    (censor == "none" & is.na(y)) |
    (censor == "left"  & is.na(left)) |
    (censor == "right" & is.na(right)) |
    (censor == "interval" & (is.na(left) | is.na(right)))

  censor[is_missing] <- "missing"
  list(censor = censor, n_missing = sum(is_missing), which_missing = which(is_missing))
}

#' Handle missing values in the fixed-effects design matrix
#'
#' @param X Numeric design matrix (already including the intercept
#'   column, if any -- the intercept column is never touched).
#' @param method One of:
#'  * `"error"` (default): stop with an informative message.
#'  * `"cca"`: complete-case analysis -- rows with any missing
#'    covariate are dropped (with a warning); response/group/censoring
#'    vectors must be subset identically by the caller.
#'  * `"mean"`: unconditional column-mean imputation (deterministic,
#'    ignores correlation between predictors -- fast but understates
#'    uncertainty; only intended for quick exploration).
#'  * `"reg"`: iterative conditional-mean (regression) imputation:
#'    each incomplete column is regressed in turn on the current
#'    values of the other columns and imputed with the fitted mean,
#'    cycling until stable or `max_iter` sweeps. This is a single
#'    imputation and, like `"mean"`, does not propagate imputation
#'    uncertainty into the final standard errors.
#' @param has_intercept Whether column 1 of `X` is an intercept
#'   (excluded from imputation modeling on the right-hand side is not
#'   necessary, but it is never itself imputed since it has no NAs).
#' @return A list with the (possibly modified) `X`, a logical `keep`
#'   vector (all `TRUE` unless `method = "cca"`), and `n_missing`.
#' @keywords internal
pecme_handle_missing_x <- function(X, method = c("error", "cca", "mean", "reg"),
                                    has_intercept = TRUE, max_iter = 10L) {
  method <- match.arg(method)
  X <- as.matrix(X)
  na_mask <- is.na(X)
  n_missing <- sum(na_mask)

  if (n_missing == 0L) {
    return(list(X = X, keep = rep(TRUE, nrow(X)), n_missing = 0L, method = method))
  }

  if (method == "error") {
    stop(n_missing, " missing value(s) found in the fixed-effects design matrix. ",
         "Set `missing_x` to \"cca\", \"mean\", or \"reg\" in pecme(), or impute ",
         "beforehand (e.g. with `mice`, then combine fits with pecme_mi_combine()).",
         call. = FALSE)
  }

  if (method == "cca") {
    keep <- rowSums(na_mask) == 0L
    if (!any(keep)) stop("`missing_x = \"cca\"` removed every row.", call. = FALSE)
    warning(sum(!keep), " row(s) dropped by complete-case analysis (`missing_x = \"cca\"`).",
            call. = FALSE)
    return(list(X = X[keep, , drop = FALSE], keep = keep, n_missing = n_missing, method = method))
  }

  cols <- seq_len(ncol(X))
  if (has_intercept) cols <- cols[-1L]
  cols <- cols[colSums(na_mask[, cols, drop = FALSE]) > 0]

  ## Mean/regression imputation is deliberately performed on the numeric
  ## model matrix. It is therefore only appropriate for suitable numeric
  ## predictors. In particular, missing factor levels would become fractional
  ## dummy values, which have no substantive interpretation.
  if (method %in% c("mean", "reg") && anyNA(X[, cols, drop = FALSE])) {
    ## The model matrix itself cannot reliably recover raw-variable classes,
    ## so callers should use raw-variable imputation (e.g. mice) whenever
    ## factors, interactions, splines, or nonlinear terms are involved.
    warning("`missing_x = \"mean\"`/`\"reg\"` imputes the numeric model matrix. ",
            "Use it only for suitable numeric predictors; for factors, interactions, ",
            "splines, or nonlinear terms, impute raw variables before fitting (e.g. with `mice`).",
            call. = FALSE)
  }

  if (method == "mean") {
    for (j in cols) {
      m <- mean(X[, j], na.rm = TRUE)
      X[is.na(X[, j]), j] <- m
    }
    return(list(X = X, keep = rep(TRUE, nrow(X)), n_missing = n_missing, method = method))
  }

  ## method == "reg": iterative conditional-mean imputation
  for (j in cols) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)  # warm start

  for (sweep in seq_len(max_iter)) {
    max_change <- 0
    for (j in cols) {
      rows_j <- which(na_mask[, j])
      if (!length(rows_j)) next
      rhs <- setdiff(seq_len(ncol(X)), c(j, if (has_intercept) 1L else integer(0)))
      fit_data <- data.frame(X[-rows_j, rhs, drop = FALSE])
      colnames(fit_data) <- paste0("v", rhs)
      fit_data$.y <- X[-rows_j, j]
      new_data <- data.frame(X[rows_j, rhs, drop = FALSE])
      colnames(new_data) <- paste0("v", rhs)

      fm <- stats::lm(stats::reformulate(paste0("v", rhs), response = ".y"), data = fit_data)
      pred <- as.numeric(stats::predict(fm, newdata = new_data))
      pred[!is.finite(pred)] <- mean(X[-rows_j, j])

      max_change <- max(max_change, max(abs(pred - X[rows_j, j])))
      X[rows_j, j] <- pred
    }
    if (max_change < 1e-6) break
  }

  list(X = X, keep = rep(TRUE, nrow(X)), n_missing = n_missing, method = method)
}

#' Combine `pecme` fits across multiply-imputed data sets (Rubin's rules)
#'
#' A convenience helper for the recommended workflow when fixed-effects
#' covariates (not the response) have missing values: impute `m`
#' complete data sets externally (e.g. with the `mice` package), fit
#' the *same* `pecme()` model to each, and combine with Rubin's rules.
#'
#' Rubin's rules require a (classically) valid standard error for every
#' fit being combined. A penalized fit (`penalty != "none"` and
#' `lambda > 0`) deliberately stores `se = NA` (see [pecme()]): naive
#' standard errors are not valid for a shrinkage estimate, so silently
#' pooling `NA`s here would produce `NA` within-imputation variance and
#' a misleading pooled result. Because of this, `pecme_mi_combine()`:
#'
#' * by default, REQUIRES every element of `fits` to have valid SEs
#'   (`se_valid == TRUE`, i.e. `penalty == "none"` or `lambda == 0`) and
#'   errors otherwise, naming which elements fail and why;
#' * with `use_refit = TRUE`, instead pools each fit's `$refit` --  the
#'   automatic unpenalized "relaxed refit" on the selected predictors
#'   (see [pecme()]'s `refit` argument) -- which does have valid SEs.
#'   This is a practical compromise, not a rigorous solution: (a) the
#'   relaxed refit is not valid unconditional post-selection inference
#'   to begin with (its SEs ignore the selection step); (b) different
#'   imputed data sets can select different variable subsets, so a
#'   coefficient absent from one imputation's refit is treated as
#'   exactly 0 (with `se = 0`) in that imputation for pooling purposes,
#'   which folds selection instability into the between-imputation
#'   variance rather than separating it out. Report results obtained
#'   this way as approximate/exploratory.
#'
#' For fully rigorous inference on a penalized, multiply-imputed fit,
#' see [pecme_bootstrap()] (a cluster bootstrap) applied within each
#' imputed data set, or consult a specialist reference on selective/
#' post-selection inference under multiple imputation.
#'
#' @param fits A list of `m` objects of class `"pecme"`, all fit with
#'   identical `penalty`/`lambda`/`alpha`/`gamma` and the same fixed
#'   coefficient names/order.
#' @param use_refit If `TRUE`, pool each fit's `$refit` sub-fit instead
#'   of `fits` themselves (see Details); required whenever `fits` are
#'   penalized.
#' @return A data frame with the pooled estimate, within-imputation
#'   variance, between-imputation variance, total variance, and pooled
#'   standard error for each fixed-effect coefficient.
#' @export
pecme_mi_combine <- function(fits, use_refit = FALSE) {
  if (!is.list(fits) || length(fits) < 2L) {
    stop("`fits` must be a list of >= 2 pecme fits.", call. = FALSE)
  }
  if (!all(vapply(fits, inherits, logical(1), "pecme"))) {
    stop("All elements of `fits` must be objects of class \"pecme\".", call. = FALSE)
  }
  m <- length(fits)

  if (use_refit) {
    missing_refit <- vapply(fits, function(f) is.null(f$refit), logical(1))
    if (any(missing_refit)) {
      stop("`use_refit = TRUE` but fit(s) ", paste(which(missing_refit), collapse = ", "),
           " have no `$refit`; fit with `refit = TRUE` before pooling.", call. = FALSE)
    }
    bad_refit <- vapply(fits, function(f) {
      r <- f$refit
      !isTRUE(r$converged) || !isTRUE(r$se_valid) ||
        is.null(r$se) || any(!is.finite(r$se))
    }, logical(1))
    if (any(bad_refit)) {
      stop("`use_refit = TRUE` requires each `$refit` to have converged with valid finite standard errors; ",
           "invalid refit(s): ", paste(which(bad_refit), collapse = ", "), call. = FALSE)
    }
    src <- lapply(fits, function(f) f$refit)
  } else {
    bad <- !vapply(fits, function(f) isTRUE(f$se_valid), logical(1))
    if (any(bad)) {
      stop(
        "pecme_mi_combine() requires classically valid standard errors, but fit(s) ",
        paste(which(bad), collapse = ", "), " do not have them (penalized fit with ",
        "lambda > 0: naive SEs for a shrinkage estimate are not valid, so pecme() ",
        "stores se = NA by design). Either fit with penalty = \"none\"/lambda = 0, ",
        "or set `use_refit = TRUE` to pool each fit's relaxed post-selection refit ",
        "instead (see ?pecme_mi_combine Details for the caveats of that approach), ",
        "or use pecme_bootstrap() within each imputed data set for a fully valid ",
        "(if more expensive) alternative.",
        call. = FALSE
      )
    }
    src <- fits
  }

  ## Union of coefficient names across all m fits: with use_refit = TRUE,
  ## different imputations can select different predictors. A term
  ## absent from a given fit's refit is treated as exactly 0 (se = 0)
  ## in that fit for pooling purposes -- see Details above.
  nm <- Reduce(union, lapply(src, function(f) names(f$coefficients)))

  B <- vapply(src, function(f) {
    v <- stats::setNames(rep(0, length(nm)), nm)
    v[names(f$coefficients)] <- f$coefficients
    v
  }, numeric(length(nm)))
  U <- vapply(src, function(f) {
    v <- stats::setNames(rep(0, length(nm)), nm)
    se <- f$se
    if (!is.null(se)) v[names(se)] <- se
    v
  }, numeric(length(nm)))^2
  rownames(B) <- rownames(U) <- nm

  qbar <- rowMeans(B)
  ubar <- rowMeans(U, na.rm = TRUE)
  bvar <- apply(B, 1L, stats::var)
  total_var <- ubar + (1 + 1 / m) * bvar

  data.frame(
    term = nm, estimate = qbar,
    within_var = ubar, between_var = bvar,
    total_var = total_var, se = sqrt(total_var),
    row.names = NULL
  )
}
