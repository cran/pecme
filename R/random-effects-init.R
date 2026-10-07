## -----------------------------------------------------------------------
## Starting values for (beta, sigma, sigma_b)
##
## `lme4::lmer()` is used ONLY to obtain numerically stable starting
## values for the ECME iterations. It never enters the ECME likelihood
## or E-step itself (lmer's own random-intercept model ignores
## censoring, so its estimates are simply a fast, well-conditioned
## initialization, not a fitted censored model).
## -----------------------------------------------------------------------

#' @keywords internal
pecme_init_random_effects <- function(data, response, fixed_formula, random_var,
                                       has_intercept = TRUE, REML = FALSE,
                                       init_method = c("lmer", "zero", "random")) {
  init_method <- match.arg(init_method)
  fit <- NULL
  beta <- NULL

  if (init_method == "zero") {
    rand_eff <- rep(0, nrow(data))
    sigma_b <- 1e-4
    sigma <- stats::sd(data[[response]], na.rm = TRUE)

  } else if (init_method == "random") {
    lev <- levels(factor(data[[random_var]]))
    G <- length(lev)
    sigma_b <- 0.10
    b_group <- stats::rnorm(G, mean = 0, sd = sigma_b)
    names(b_group) <- lev
    rand_eff <- b_group[as.character(data[[random_var]])]
    sigma <- stats::sd(data[[response]], na.rm = TRUE)

  } else {
    term_labels <- attr(stats::terms(fixed_formula), "term.labels")
    rhs <- if (length(term_labels)) paste(term_labels, collapse = " + ") else "1"
    if (!has_intercept) rhs <- paste(rhs, "- 1")
    form <- stats::as.formula(paste(response, "~", rhs, "+ (1|", random_var, ")"))
    ## A singular-fit ("isSingular") warning here just means lmer's own
    ## (censoring-blind) random-intercept variance landed near the
    ## boundary -- harmless for a starting-value fit, since sigma_b is
    ## floored below regardless and the ECME iterations re-estimate it
    ## properly under the actual censored/missing-data likelihood. The
    ## check is disabled at the source (rather than via a blanket
    ## suppressWarnings()) so genuine optimizer non-convergence warnings
    ## would still surface.
    ctrl <- lme4::lmerControl(check.conv.singular = "ignore")
    fit <- lme4::lmer(form, data = data, REML = REML, control = ctrl)
    re_df <- lme4::ranef(fit)[[random_var]]
    rand_eff <- re_df[match(as.character(data[[random_var]]), rownames(re_df)), 1]
    sigma_b <- as.numeric(attr(lme4::VarCorr(fit)[[random_var]], "stddev"))
    sigma <- as.numeric(stats::sigma(fit))
    beta <- as.numeric(lme4::fixef(fit))
  }

  list(fit = fit, beta = beta, rand_eff = as.numeric(rand_eff),
       sigma_b = max(sigma_b, 1e-8), sigma = max(sigma, 1e-8),
       init_method = init_method)
}
