## -----------------------------------------------------------------------
## Variance-component updates
##
## `pecme_variance_em()` is the ordinary EM closed-form CM-step, using
## only the E-step sufficient statistics (fast, always available as a
## fallback/starting value).
##
## `pecme_variance_ecme()` is the actual ECME step: with beta held
## fixed, it maximizes the REAL observed-data marginal likelihood
## (random intercept integrated out) over (sigma, sigma_b) by direct
## numerical optimization on the log scale. This is what distinguishes
## ECME from plain ECM/EM and is known to accelerate convergence of the
## variance components in censored/incomplete-data Gaussian mixed
## models (Liu & Rubin, 1994; Vaida & Liu, 2009).
## -----------------------------------------------------------------------

#' @keywords internal
pecme_variance_em <- function(moments, X, beta, group) {
  X <- as.matrix(X); beta <- as.numeric(beta)
  eta <- as.vector(X %*% beta)
  n <- nrow(X)
  if (length(moments$Ey) != n) stop("E-step moments have incorrect length.", call. = FALSE)

  var_y <- pmax(moments$Ey2 - moments$Ey^2, 0)
  var_b <- pmax(moments$Ebb - moments$Eb^2, 0)
  cov_yb <- moments$Eyb - moments$Ey * moments$Eb
  mean_resid <- moments$Ey - eta - moments$Eb

  ## E[(Y* - X beta - b)^2 | data] via the law of total variance:
  ## Var(Y*-b|data) + (E[Y*-b|data] - X beta)^2
  residual_sq <- pmax(var_y + var_b - 2 * cov_yb + mean_resid^2, 0)

  sigma2_new   <- max(mean(residual_sq, na.rm = TRUE), 1e-12)
  sigma_b2_new <- max(mean(moments$group_Ebb, na.rm = TRUE), 1e-12)

  list(sigma = sqrt(sigma2_new), sigma_b = sqrt(sigma_b2_new),
       residual_sq = residual_sq, var_y = var_y, var_b = var_b,
       cov_yb = cov_yb, mean_resid = mean_resid, mu_fixed = eta)
}

#' @keywords internal
pecme_variance_ecme <- function(y, X, group, beta, sigma, sigma_b,
                                 censor, left = NULL, right = NULL,
                                 quadrature = 40L, maxit = 60L,
                                 method = c("L-BFGS-B", "Nelder-Mead"),
                                 reltol = 1e-7) {

  method <- match.arg(method)

  n <- length(y)

  if (is.null(left)) {
    left <- rep(NA_real_, n)
  }

  if (is.null(right)) {
    right <- rep(NA_real_, n)
  }

  objective <- function(log_par) {

    s <- exp(log_par[1L])
    sb <- exp(log_par[2L])

    ll <- pecme_marginal_loglik(
      y = y,
      X = X,
      group = group,
      beta = beta,
      sigma = s,
      sigma_b = sb,
      censor = censor,
      left = left,
      right = right,
      quadrature = quadrature
    )

    if (!is.finite(ll)) {
      return(1e100)
    }

    -ll
  }

  ## ---------------------------------------------------------------
  ## Initial values on log scale
  ## ---------------------------------------------------------------

  start <- log(
    c(
      max(sigma, 1e-8),
      max(sigma_b, 1e-8)
    )
  )

  ## ---------------------------------------------------------------
  ## Bounds for L-BFGS-B
  ## ---------------------------------------------------------------

  y_scale <- max(
    stats::sd(y, na.rm = TRUE),
    1,
    na.rm = TRUE
  )

  lower <- log(1e-8)
  upper <- log(max(y_scale * 100, 1e2))

  start <- pmin(
    pmax(start, lower + 1e-8),
    upper - 1e-8
  )

  ## ---------------------------------------------------------------
  ## First optimization
  ## ---------------------------------------------------------------

  if (method == "L-BFGS-B") {

    opt <- stats::optim(
      par = start,
      fn = objective,
      method = "L-BFGS-B",
      lower = c(lower, lower),
      upper = c(upper, upper),
      control = list(
        maxit = as.integer(maxit),
        factr = 1e7,
        pgtol = reltol
      )
    )

    ## -------------------------------------------------------------
    ## Retry with a larger iteration budget if needed
    ## -------------------------------------------------------------

    if (
      is.finite(opt$value) &&
      any(is.finite(opt$par)) &&
      !identical(as.integer(opt$convergence), 0L)
    ) {

      retry_maxit <- max(
        as.integer(maxit) * 5L,
        300L
      )

      opt_retry <- stats::optim(
        par = opt$par,
        fn = objective,
        method = "L-BFGS-B",
        lower = c(lower, lower),
        upper = c(upper, upper),
        control = list(
          maxit = retry_maxit,
          factr = 1e7,
          pgtol = reltol
        )
      )

      ## Keep the better finite solution
      if (
        is.finite(opt_retry$value) &&
        (
          !is.finite(opt$value) ||
          opt_retry$value <= opt$value
        )
      ) {
        opt <- opt_retry
      }
    }

    ## -------------------------------------------------------------
    ## Final fallback: Nelder-Mead
    ## -------------------------------------------------------------

    if (
      !identical(as.integer(opt$convergence), 0L)
    ) {

      opt_nm <- stats::optim(
        par = opt$par,
        fn = objective,
        method = "Nelder-Mead",
        control = list(
          maxit = max(
            as.integer(maxit) * 5L,
            300L
          ),
          reltol = reltol
        )
      )

      if (
        is.finite(opt_nm$value) &&
        (
          !is.finite(opt$value) ||
          opt_nm$value <= opt$value
        )
      ) {
        opt <- opt_nm
      }
    }

  } else {

    ## -------------------------------------------------------------
    ## Direct Nelder-Mead
    ## -------------------------------------------------------------

    opt <- stats::optim(
      par = start,
      fn = objective,
      method = "Nelder-Mead",
      control = list(
        maxit = as.integer(maxit),
        reltol = reltol
      )
    )
  }

  ## ---------------------------------------------------------------
  ## Validate optimization result
  ## ---------------------------------------------------------------

  if (
    !is.finite(opt$value) ||
    any(!is.finite(opt$par))
  ) {

    return(
      list(
        sigma = sigma,
        sigma_b = sigma_b,
        convergence = 1L,
        converged = FALSE,
        iterations = 0L,
        observed_loglik = NA_real_,
        optimizer = opt
      )
    )
  }

  ## ---------------------------------------------------------------
  ## Final result
  ## ---------------------------------------------------------------

  list(
    sigma = exp(opt$par[1L]),
    sigma_b = exp(opt$par[2L]),
    convergence = as.integer(opt$convergence),

    converged = identical(
      as.integer(opt$convergence),
      0L
    ),

    iterations = as.integer(opt$counts[1L]),

    observed_loglik = -opt$value,

    optimizer = opt
  )
}
