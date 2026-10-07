## -----------------------------------------------------------------------
## S3 methods for class "pecme"
## -----------------------------------------------------------------------

#' @export
print.pecme <- function(x, digits = 4L, ...) {
  cat("Penalized ECME fit (class \"pecme\")\n")
  cat("Penalty:", x$penalty)
  if (x$penalty == "elastic") cat(sprintf(" (alpha = %.3g)", x$alpha))
  if (x$penalty %in% c("scad", "mcp")) cat(sprintf(" (gamma = %.3g)", x$gamma))
  cat(sprintf("  |  lambda = %.6g\n", x$lambda))
  cat(sprintf("Groups: %d   Observations: %d   Informative responses: %d   df: %d\n",
              x$n_groups, x$nobs, x$n_informative, x$df))
  if (x$missing_y || x$n_missing_x > 0) {
    cat(sprintf("Missing data: %d response value(s) treated as unobserved and imputed only for the E-step, ",
                x$n_missing_y),
        sprintf("%d covariate value(s) handled via missing_x = \"%s\"\n", x$n_missing_x, x$missing_x))
  }
  cat(sprintf("Converged: %s (%s) in %d iteration(s)\n", x$converged, x$reason, x$iterations))
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  cat(sprintf("\nsigma = %.*f   sigma_b = %.*f\n", digits, x$sigma, digits, x$sigma_b))
  cat(sprintf("logLik = %.*f   AIC = %.*f   BIC = %.*f\n",
              digits, x$logLik, digits, x$AIC, digits, x$BIC))
  if (x$penalty != "none" && x$lambda > 0) {
    cat("Selected (nonzero) predictors:", if (length(x$selected_names)) paste(x$selected_names, collapse = ", ") else "(none)", "\n")
  }
  invisible(x)
}

#' @export
summary.pecme <- function(object, ...) {
  structure(list(fit = object), class = "summary.pecme")
}

#' @export
print.summary.pecme <- function(x, digits = 4L, ...) {
  object <- x$fit
  cat("Call:\n")
  print(object$call)
  cat("\n")
  print.pecme(object, digits = digits)

  if (object$se_valid) {
    z <- object$coefficients / object$se
    pval <- 2 * stats::pnorm(-abs(z))
    tab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se,
                 `z value` = z, `Pr(>|z|)` = pval)
    cat("\nCoefficient table (conditional on sigma, sigma_b; unpenalized fit):\n")
    stats::printCoefmat(tab, digits = digits)
  } else {
    cat("\nStandard errors are not reported for the penalized fit ",
        "(naive SEs are not classically valid for a shrinkage estimate). ",
        "See `$refit` for a relaxed (unpenalized, selected-predictors-only) ",
        "refit, or set `refit = TRUE` in pecme()/grid_pecme().\n", sep = "")
    if (!is.null(object$refit)) {
      cat("\nRelaxed refit (unpenalized, on the selected predictors only --\n",
          "reduces shrinkage bias, but its SE/p-values do NOT account for\n",
          "the selection step and so are not valid unconditional\n",
          "post-selection inference; treat as approximate/exploratory):\n", sep = "")
      z <- object$refit$coefficients / object$refit$se
      pval <- 2 * stats::pnorm(-abs(z))
      tab <- cbind(Estimate = object$refit$coefficients, `Std. Error` = object$refit$se,
                   `z value` = z, `Pr(>|z|)` = pval)
      stats::printCoefmat(tab, digits = digits)
    }
  }
  invisible(x)
}

#' @export
coef.pecme <- function(object, use_refit = FALSE, ...) {
  if (use_refit && !is.null(object$refit)) {
    
    out <- stats::setNames(
      rep(0, length(object$coefficients)),
      names(object$coefficients)
    )
    
    refit_coef <- object$refit$coefficients
    
    ## Map refit coefficient names back to original predictor names
    if (!is.null(object$refit$refit_name_map)) {
      
      refit_name_map <- object$refit$refit_name_map
      
      mapped_names <- names(refit_coef)
      
      idx <- mapped_names %in% names(refit_name_map)
      
      mapped_names[idx] <-
        unname(refit_name_map[mapped_names[idx]])
      
      names(refit_coef) <- mapped_names
    }
    
    ## Insert refitted coefficients into the original coefficient vector
    common_names <- intersect(names(out), names(refit_coef))
    
    out[common_names] <- refit_coef[common_names]
    
    return(out)
  }
  
  object$coefficients
}

#' @export
fitted.pecme <- function(object, ...) {
  eta <- as.vector(object$X %*% object$coefficients)
  b <- object$random_coef[match(object$group_levels[object$group], names(object$random_coef))]
  eta + b
}

#' @export
residuals.pecme <- function(object, ...) {
  object$y - fitted.pecme(object)
}

#' @export
logLik.pecme <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$df
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' @export
nobs.pecme <- function(object, ...) object$nobs

#' @export
AIC.pecme <- function(object, ..., k = 2) -2 * object$logLik + k * object$df

#' @export
BIC.pecme <- function(object, ...) object$BIC

#' Predict from a fitted pecme model
#'
#' @param object A `"pecme"` fit.
#' @param newdata Optional data frame with the predictor columns (and,
#'   for `level = "group"`, the grouping column) used to fit `object`.
#'   If omitted, fitted values on the training data are returned.
#' @param level `"population"` (fixed effects only, `X beta`) or
#'   `"group"` (`X beta + b_j` for known groups; groups absent from the
#'   training data effectively use `b_j = 0`, i.e. the population mean).
#' @param use_refit Use the post-selection refit coefficients (see
#'   `object$refit`) instead of the penalized estimates, if available.
#' @param se_fit Also return standard errors (only available when
#'   `object$se_valid` is `TRUE`, i.e. an unpenalized/lambda=0 fit; a
#'   naive delta-method approximation that ignores `sigma_b`
#'   uncertainty).
#' @param ... Additional arguments passed to the prediction method.
#' @return If `se_fit = FALSE`, a numeric vector of predictions.
#'   If `se_fit = TRUE`, a list with components `fit` and `se.fit`,
#'   containing the predicted values and their estimated standard errors,
#'   respectively.
#' @export
predict.pecme <- function(object, newdata = NULL,
                           level = c("population", "group"),
                           use_refit = FALSE, se_fit = FALSE, ...) {
  level <- match.arg(level)
  beta <- coef.pecme(object, use_refit = use_refit)
  

  if (is.null(newdata)) {
    X <- object$X
    group_id <- object$group_levels[object$group]
  } else if (isTRUE(object$refit_is_relaxed) &&
             !is.null(object$refit_parent_formula)) {
    ## Relaxed refits are internally built on synthetic v1/v2/... columns.
    ## For user-facing prediction, reconstruct the ORIGINAL model matrix
    ## from the parent formula/data semantics and map selected columns.
    parent_terms <- object$refit_parent_terms
    mf_parent <- stats::model.frame(stats::delete.response(parent_terms),
                                    newdata,
                                    xlev = object$refit_parent_xlevels,
                                    na.action = stats::na.pass)
    X_parent <- stats::model.matrix(stats::delete.response(parent_terms),
                                    mf_parent,
                                    contrasts.arg = object$refit_parent_contrasts)
    parent_cols <- object$refit_parent_model_matrix_cols
    selected_cols <- names(object$coefficients)
    ## refit coefficients are named with original selected coefficient names.
    X <- matrix(0, nrow(X_parent), length(selected_cols),
                dimnames = list(rownames(X_parent), selected_cols))
    common <- intersect(colnames(X_parent), selected_cols)
    if (length(common)) X[, common] <- X_parent[, common, drop = FALSE]
    beta <- beta[selected_cols]
    group_id <- if (level == "group") {
      rv <- object$refit_parent_random_var
      if (is.null(rv) || !rv %in% names(newdata)) {
        stop("`newdata` must contain the original grouping column ", rv,
             " when level = 'group'.", call. = FALSE)
      }
      as.character(newdata[[rv]])
    } else NULL
  } else {
    mf <- stats::model.frame(stats::delete.response(object$terms), newdata,
                              xlev = object$xlevels, na.action = stats::na.pass)
    X_new <- stats::model.matrix(stats::delete.response(object$terms), mf, contrasts.arg = object$contrasts)
    
    model_matrix_cols <- object$model_matrix_cols
    
    if (is.null(model_matrix_cols)) {
      ## Backward compatibility for objects created before this field existed
      model_matrix_cols <- names(beta)
    }
    
    X <- matrix(
      0,
      nrow = nrow(X_new),
      ncol = length(model_matrix_cols),
      dimnames = list(
        rownames(X_new),
        model_matrix_cols
      )
    )
    
    common <- intersect(
      colnames(X_new),
      model_matrix_cols
    )
    
    if (length(common)) {
      X[, common] <- X_new[, common, drop = FALSE]
    }
    
    beta <- beta[model_matrix_cols]
    group_id <- if (level == "group") {
      if (is.null(object$random_var) || is.null(newdata[[object$random_var]])) {
        stop("`newdata` must contain the grouping column \"", object$random_var,
             "\" when level = \"group\".", call. = FALSE)
      }
      as.character(newdata[[object$random_var]])
    } else {
      NULL
    }
  }

  eta <- as.vector(X %*% beta)

  if (level == "group") {
    b <- object$random_coef[group_id]
    b[is.na(b)] <- 0  # unseen groups: predict at the population mean
    eta <- eta + b
  }

  if (!se_fit) return(eta)

  if (!isTRUE(object$se_valid)) {
    warning("Standard errors are not available for a penalized fit; returning NA. ",
            "Fit with penalty = \"none\" or lambda = 0 for valid SEs, or use `object$refit`.",
            call. = FALSE)
    return(list(fit = eta, se.fit = rep(NA_real_, length(eta))))
  }
  if (is.null(object$covariance)) {
    warning("The fitted object has no coefficient covariance matrix; returning NA prediction SEs.", call. = FALSE)
    return(list(fit = eta, se.fit = rep(NA_real_, length(eta))))
  }
  V <- object$covariance[colnames(X), colnames(X), drop = FALSE]
  var_eta <- rowSums((X %*% V) * X)
  var_eta[!is.finite(var_eta) | var_eta < 0] <- NA_real_
  se_eta <- sqrt(var_eta)
  list(fit = eta, se.fit = se_eta)
}

#' @export
plot.pecme <- function(x, ...) {
  fv <- fitted.pecme(x)
  rs <- residuals.pecme(x)
  op <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(fv, rs, xlab = "Fitted", ylab = "Residuals",
                 main = "Residuals vs Fitted", pch = 20)
  graphics::abline(h = 0, lty = 2, col = "grey50")

  stats::qqnorm(rs, main = "Normal Q-Q (residuals)")
  stats::qqline(rs, col = "grey50")

  b <- x$random_coef
  stats::qqnorm(b, main = "Normal Q-Q (random intercepts)")
  stats::qqline(b, col = "grey50")

  if (length(x$objective_history) > 1L) {
    graphics::plot(x$objective_history, type = "b", pch = 20,
                   xlab = "ECME iteration", ylab = "Penalized objective",
                   main = "Convergence")
  } else {
    graphics::plot.new()
  }
  invisible(x)
}

## --- pecme_grid ---------------------------------------------------------

#' @export
print.pecme_grid <- function(x, digits = 4L, ...) {
  cat("Penalized ECME grid search (class \"pecme_grid\")\n")
  cat("Penalty:", x$penalty, " | Criterion:", x$criterion,
      " | Grid size:", length(x$lambda), "\n")
  cat(sprintf("Selected lambda = %.6g (index %d of %d)\n",
              x$lambda[x$best_index], x$best_index, length(x$lambda)))
  cat("\n")
  print.pecme(x$best, digits = digits)
  invisible(x)
}

#' @export
summary.pecme_grid <- function(object, ...) {
  structure(list(grid = object), class = "summary.pecme_grid")
}

#' @export
print.summary.pecme_grid <- function(x, digits = 4L, ...) {
  print.pecme_grid(x$grid, digits = digits)
  cat("\nScores by lambda:\n")
  sc <- x$grid$scores
  sc[-1] <- lapply(sc[-1], function(col) if (is.numeric(col)) round(col, digits) else col)
  print(sc)
  invisible(x)
}

#' Plot the regularization path and/or selection-criterion curve
#'
#' @param x A `"pecme_grid"` object.
#' @param which `"path"` (coefficient paths vs `log(lambda)`),
#'   `"criterion"` (the selection criterion vs `log(lambda)`), or
#'   `"both"` (default).
#' @param ... Additional graphical arguments. Currently ignored.
#' @return Invisibly returns the input `x`, an object of class
#'   `"pecme_grid"`.
#' @export
plot.pecme_grid <- function(x, which = c("both", "path", "criterion"), ...) {
  which <- match.arg(which)
  log_lambda <- log(x$lambda)

  draw_path <- function() {
    if (is.null(x$path)) { graphics::plot.new(); return(invisible()) }
    graphics::matplot(log_lambda, x$path, type = "l", lty = 1,
                       xlab = expression(log(lambda)), ylab = "Coefficient",
                       main = "Regularization path")
    graphics::abline(v = log(x$lambda[x$best_index]), lty = 2, col = "grey40")
  }
  draw_criterion <- function() {
    crit_col <- if (x$criterion == "cv") "CV_logLik" else toupper(x$criterion)
    y <- x$scores[[crit_col]]
    graphics::plot(log_lambda, y, type = "b", pch = 20,
                   xlab = expression(log(lambda)), ylab = toupper(x$criterion),
                   main = paste("Selection criterion:", toupper(x$criterion)))
    graphics::abline(v = log(x$lambda[x$best_index]), lty = 2, col = "grey40")
  }

  if (which == "path") { draw_path(); return(invisible(x)) }
  if (which == "criterion") { draw_criterion(); return(invisible(x)) }

  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  draw_path(); draw_criterion()
  invisible(x)
}
