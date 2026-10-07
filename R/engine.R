## -----------------------------------------------------------------------
## Core ECME iteration for one (penalty, lambda, alpha, gamma)
##
## Operates entirely in the (already standardized) design-matrix scale;
## the user-facing pecme() wrapper in fit.R handles formula parsing,
## standardization, and back-transforming the fitted coefficients.
##
## One outer iteration =
##   E-step                      (estep.R;      GH quadrature; numerical)
##   CM-step 1: beta             (coordinate-descent.R; penalized Q-step)
##   CM-step 2: sigma, sigma_b   (variance-update.R;     "ecme" or "em")
## -----------------------------------------------------------------------

#' @keywords internal
pecme_fit_engine <- function(X, y, group, censor, left, right,
                              beta_init, sigma_init, sigma_b_init,
                              penalty = .pecme_penalties, lambda = 0,
                              alpha = 0.5, gamma = 3.7, weights = NULL,
                              has_intercept = TRUE, quadrature = 40L,
                              max_iter = 200L, tol = 1e-6,
                              cd_max_iter = 1000L, cd_tol = 1e-8,
                              variance_step = c("ecme", "em"),
                              variance_maxit = 60L,
                              workers = 1L, cluster = NULL,
                              check_monotonicity = FALSE, monotonicity_tol = 1e-8,
                             verbose = FALSE) {
  penalty <- match.arg(penalty)
  variance_step <- match.arg(variance_step)

  beta <- as.numeric(beta_init)
  sigma <- as.numeric(sigma_init)
  sigma_b <- as.numeric(sigma_b_init)

  ## NA initialization is essential: an optimizer failure must never be
  ## represented by the artificial value 0 in the returned likelihood or
  ## objective history.
  obj_hist <- rep(NA_real_, max_iter)
  ll_hist  <- rep(NA_real_, max_iter)
  monotonicity_violation <- rep(FALSE, max_iter)
  
  converged <- FALSE
  outer_converged <- FALSE
  
  ## Track convergence of the inner coordinate-descent steps
  inner_converged <- logical(max_iter)
  inner_iterations <- integer(max_iter)
  inner_failed <- logical(max_iter)
  
  variance_converged <- logical(max_iter)
  variance_iterations <- integer(max_iter)
  variance_failed <- logical(max_iter)
  
  reason <- "max_iterations"
  it_used <- 0L

  for (iter in seq_len(max_iter)) {
    beta_prev <- beta; sigma_prev <- sigma; sigma_b_prev <- sigma_b

    ## --- E-step -----------------------------------------------------
    moments <- pecme_e_step_dispatch(
      y = y, X = X, beta = beta, group = group, sigma = sigma, sigma_b = sigma_b,
      censor = censor, left = left, right = right, quadrature = quadrature,
      workers = workers, cluster = cluster
    )
    y_adj <- moments$Ey - moments$Eb

    ## --- CM-step 1: beta (penalized Q-step) --------------------------
    cd <- pecme_coordinate_descent(
      X = X, y = y_adj, beta = beta, sigma = sigma, lambda = lambda,
      penalty = penalty, weights = weights, alpha = alpha, gamma = gamma,
      max_iter = cd_max_iter, tol = cd_tol, has_intercept = has_intercept,
      check_monotonicity = check_monotonicity
    )
    
    ## Record inner coordinate-descent convergence
    cd_ok <- isTRUE(cd$converged)
    
    inner_converged[iter] <- cd_ok
    
    ## coordinate-descent returns `iter` (not `iterations`).
    if (!is.null(cd$iter)) {
      inner_iterations[iter] <- as.integer(cd$iter)
    }
    
    inner_failed[iter] <- !cd_ok
    
    if (!cd_ok) {
      warning(
        "Coordinate-descent failed to converge at outer iteration ",
        iter,
        " (maximum inner iterations = ",
        cd_max_iter,
        ").",
        call. = FALSE
      )
    }
    
    beta <- cd$beta

    ## --- CM-step 2: sigma, sigma_b -----------------------------------
    if (variance_step == "ecme") {
      
      vout <- pecme_variance_ecme(
        y = y,
        X = X,
        group = group,
        beta = beta,
        sigma = sigma,
        sigma_b = sigma_b,
        censor = censor,
        left = left,
        right = right,
        quadrature = quadrature,
        maxit = variance_maxit
      )
      
      ## Check variance-optimizer convergence
      variance_ok <- isTRUE(vout$converged)
      
      variance_converged[iter] <- variance_ok
      
      if (!is.null(vout$iterations)) {
        variance_iterations[iter] <- as.integer(vout$iterations)
      }
      
      if (!variance_ok) {
        
        variance_failed[iter] <- TRUE
        
        warning(
          "Variance optimizer failed to converge at outer iteration ",
          iter,
          " after ",
          variance_maxit,
          " iterations. ",
          "The outer ECME iteration is terminated.",
          call. = FALSE
        )
        
        ## Retain the returned values for diagnostics,
        ## but do not regard this ECME iteration as successful.
        sigma <- vout$sigma
        sigma_b <- vout$sigma_b
        
        it_used <- iter
        reason <- "variance_optimizer_failure"
        
        break
      }
      
      sigma <- vout$sigma
      sigma_b <- vout$sigma_b
    } else {
      ## Re-run the E-step conditional on the NEW beta so the classical
      ## EM closed-form CM-step is evaluated at the current iterate.
      moments2 <- pecme_e_step_dispatch(
        y = y, X = X, beta = beta, group = group, sigma = sigma, sigma_b = sigma_b,
        censor = censor, left = left, right = right, quadrature = quadrature,
        workers = workers, cluster = cluster
      )
      vout <- pecme_variance_em(moments2, X, beta, group)
      sigma <- vout$sigma
      sigma_b <- vout$sigma_b
      ## The EM variance CM-step is closed form; reaching this line means
      ## that its update succeeded.
      variance_converged[iter] <- TRUE
      variance_iterations[iter] <- 1L
      variance_failed[iter] <- FALSE
    }

    ## --- bookkeeping / convergence ------------------------------------
    obj <- pecme_objective(
      beta, y, X, sigma, sigma_b, group, censor, left, right,
      penalty, lambda, alpha, gamma, weights, has_intercept, quadrature
    )
    
    obj_hist[iter] <- obj$objective
    ll_hist[iter]  <- obj$logLik
    it_used <- iter
    
    if (iter > 1L) {
      
      objective_change <- obj_hist[iter] - obj_hist[iter - 1L]
      
      allowed_decrease <-
        monotonicity_tol *
        max(1, abs(obj_hist[iter - 1L]))
      
      ## Q(θ) = -logLik(θ) + penalty(β) is minimized. Therefore
      ## monotonicity is violated only by a positive increase beyond tolerance.
      if (objective_change > allowed_decrease) {
        monotonicity_violation[iter] <- TRUE
        if (check_monotonicity) {
          warning(
            "Penalized objective increased at outer iteration ", iter,
            ": previous = ", format(obj_hist[iter - 1L], digits = 8),
            ", current = ", format(obj_hist[iter], digits = 8),
            ", increase = ", format(objective_change, digits = 8), ".",
            call. = FALSE
          )
        }
      }
    }
    it_used <- iter

    if (verbose) {
      message(sprintf("iter %3d | objective %.6f | logLik %.6f | sigma %.5f | sigma_b %.5f",
                       iter, obj$objective, obj$logLik, sigma, sigma_b))
    }

    beta_diff    <- max(abs(beta - beta_prev))
    sigma_diff   <- abs(sigma - sigma_prev)
    sigma_b_diff <- abs(sigma_b - sigma_b_prev)
    obj_diff <- if (iter > 1L) abs(obj_hist[iter] - obj_hist[iter - 1L]) / max(abs(obj_hist[iter - 1L]), 1) else Inf

    if (!is.finite(obj$objective)) {
      reason <- "nonfinite_objective"
      break
    }

    if (iter > 1L &&
        beta_diff < tol &&
        sigma_diff < tol &&
        sigma_b_diff < tol &&
        obj_diff < tol) {
      
      outer_converged <- TRUE
      
      ## Overall convergence requires successful inner optimization
      if (all(inner_converged[seq_len(iter)])) {
        converged <- TRUE
        reason <- "converged"
        break
      } else {
        converged <- FALSE
        reason <- "outer_converged_inner_failure"
        break
      }
    }
  }

  list(
    beta = beta,
    sigma = sigma,
    sigma_b = sigma_b,
    
    iterations = it_used,
    
    ## Convergence diagnostics
    converged = converged,
    outer_converged = outer_converged,
    inner_converged = if (it_used > 0L) {
      all(inner_converged[seq_len(it_used)])
    } else {
      FALSE
    },
    
    inner_failed = if (it_used > 0L) {
      any(inner_failed[seq_len(it_used)])
    } else {
      FALSE
    },
    
    inner_failure_count = if (it_used > 0L) {
      sum(inner_failed[seq_len(it_used)])
    } else {
      0L
    },
    
    inner_converged_history =
      inner_converged[seq_len(it_used)],
    
    inner_iterations =
      inner_iterations[seq_len(it_used)],
    
    reason = reason,
    
    objective = obj_hist[it_used],
    logLik = ll_hist[it_used],
    
    objective_history = obj_hist[seq_len(it_used)],
    logLik_history = ll_hist[seq_len(it_used)],
    
    penalty = penalty,
    lambda = lambda,
    alpha = alpha,
    gamma = gamma,
    weights = weights,
    variance_step = variance_step,
    variance_converged = if (it_used > 0L) {
      all(variance_converged[seq_len(it_used)])
    } else {
      FALSE
    },
    
    variance_failed = if (it_used > 0L) {
      any(variance_failed[seq_len(it_used)])
    } else {
      FALSE
    },
    
    variance_failure_count = if (it_used > 0L) {
      sum(variance_failed[seq_len(it_used)])
    } else {
      0L
    },
    
    variance_converged_history =
      variance_converged[seq_len(it_used)],
    
    variance_iterations =
      variance_iterations[seq_len(it_used)],
    monotonicity_violations = sum(
      monotonicity_violation[seq_len(it_used)]
    ),
    
    monotonicity_violation_history =
      monotonicity_violation[seq_len(it_used)],
    
    monotonicity_checked = check_monotonicity,
    
    monotonicity_tol = monotonicity_tol
  )
}
