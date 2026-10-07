## -----------------------------------------------------------------------
## Conditional moments of Y* | B = b under the pecme observation model
##
## Convention (fixed relative to the original draft, where left/right
## were swapped):
##   "none"     Y* = y                          (fully observed)
##   "left"     Y* <= left                       (left-censored)
##   "right"    Y* >= right                      (right-censored)
##   "interval" left < Y* <= right                (interval-censored)
##   "missing"  Y* completely unobserved, MAR given (X, group, model)
##
## "missing" is handled as its own case rather than as an
## interval with (-Inf, Inf) bounds: R evaluates (-Inf) * 0 as NaN,
## which the generic truncated-interval algebra would hit whenever a
## boundary term such as a * dnorm(a) is evaluated at a = -Inf.
## -----------------------------------------------------------------------

#' Conditional first/second moments of (possibly censored or missing) Y*
#'
#' Vectorized over observations and, optionally, over a second dimension
#' (quadrature nodes for the random intercept): `mu` may be a matrix of
#' shape `n x K`, in which case moments are returned as `n x K` matrices.
#'
#' @keywords internal
pecme_truncated_moments <- function(mu, sigma, censor, left, right, y) {
  mu_dim <- dim(mu)
  if (is.null(mu_dim)) {
    mu <- matrix(mu, nrow = length(mu), ncol = 1L)
    mu_dim <- dim(mu)
  }
  n <- mu_dim[1L]; K <- mu_dim[2L]

  censor <- as.character(censor)
  if (length(censor) != n) stop("censor length mismatch.", call. = FALSE)
  if (length(left)  != n) left  <- rep(NA_real_, n)
  if (length(right) != n) right <- rep(NA_real_, n)
  if (length(y)     != n) y     <- rep(NA_real_, n)

  m1 <- matrix(as.numeric(y), nrow = n, ncol = K)
  m2 <- m1^2

  ## --- fully observed --------------------------------------------------
  ## nothing to do: m1/m2 already hold y and y^2.

  ## --- missing: unconstrained conditional mean/variance -----------------
  ii <- which(censor == "missing")
  if (length(ii)) {
    m1[ii, ] <- mu[ii, , drop = FALSE]
    m2[ii, ] <- mu[ii, , drop = FALSE]^2 + sigma^2
  }

  ## --- left censoring: Y* <= L -------------------------------------------
  ii <- which(censor == "left")
  if (length(ii)) {
    a <- (matrix(left[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    logZ <- stats::pnorm(a, log.p = TRUE)
    lambda <- exp(stats::dnorm(a, log = TRUE) - logZ)
    mm <- mu[ii, , drop = FALSE] - sigma * lambda
    vv <- pmax(sigma^2 * (1 - a * lambda - lambda^2), 0)
    m1[ii, ] <- mm; m2[ii, ] <- vv + mm^2
  }

  ## --- right censoring: Y* >= R -------------------------------------------
  ii <- which(censor == "right")
  if (length(ii)) {
    a <- (matrix(right[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    logZ <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
    lambda <- exp(stats::dnorm(a, log = TRUE) - logZ)
    mm <- mu[ii, , drop = FALSE] + sigma * lambda
    vv <- pmax(sigma^2 * (1 + a * lambda - lambda^2), 0)
    m1[ii, ] <- mm; m2[ii, ] <- vv + mm^2
  }

  ## --- interval censoring: L < Y* <= R -------------------------------------
  ii <- which(censor == "interval")
  if (length(ii)) {
    a <- (matrix(left[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    b <- (matrix(right[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    logZ <- matrix(.pecme_interval_logprob(as.vector(a), as.vector(b)), nrow = length(ii))

    ## Stable log-scale Mills-ratio form: each of phi(a)/Z and phi(b)/Z
    ## is computed as a single exp(bounded log-difference), exactly as
    ## for left/right censoring above. The previous implementation
    ## instead formed Zinv = exp(-logZ) (which can overflow to Inf for
    ## a narrow/deep-tail interval, where logZ is very negative) and
    ## multiplied it by ordinary-scale dnorm(a)/dnorm(b) (which can
    ## underflow to exactly 0 for extreme a/b) -- an Inf * 0 = NaN
    ## failure mode well before either factor individually overflows.
    ## phi(a)/Z and phi(b)/Z are each mathematically bounded (they are
    ## themselves valid truncated-normal-type ratios), so computing them
    ## this way keeps the whole calculation finite even at z = +-30 or
    ## beyond, unlike the raw-scale product.
    log_phi_a <- stats::dnorm(a, log = TRUE)
    log_phi_b <- stats::dnorm(b, log = TRUE)
    ratio_a <- exp(log_phi_a - logZ)
    ratio_b <- exp(log_phi_b - logZ)

    delta <- ratio_a - ratio_b
    mm <- mu[ii, , drop = FALSE] + sigma * delta
    vv <- pmax(sigma^2 * (1 + (a * ratio_a - b * ratio_b) - delta^2), 0)
    m1[ii, ] <- mm; m2[ii, ] <- vv + mm^2
  }

  bad <- !(censor %in% c("none", "left", "right", "interval", "missing"))
  if (any(bad)) {
    stop("Unsupported censoring label(s): ",
         paste(unique(censor[bad]), collapse = ", "), call. = FALSE)
  }

  list(mean = m1, second = m2)
}

#' Conditional log-likelihood contribution log P(observation | B = b)
#'
#' Vectorized over observations; `mu` may be a matrix (`n x K`) for
#' simultaneous evaluation at several quadrature nodes.
#'
#' @keywords internal
pecme_conditional_loglik <- function(y, mu, sigma, censor, left, right) {
  mu_dim <- dim(mu)
  if (is.null(mu_dim)) {
    mu <- matrix(mu, nrow = length(mu), ncol = 1L)
    mu_dim <- dim(mu)
  }
  n <- mu_dim[1L]; K <- mu_dim[2L]
  censor <- as.character(censor)

  logL <- matrix(0, n, K)

  none <- which(censor == "none")
  if (length(none)) {
    logL[none, ] <- stats::dnorm(
      matrix(y[none], length(none), K),
      mean = mu[none, , drop = FALSE], sd = sigma, log = TRUE
    )
  }

  ## A completely missing response integrates to probability 1: log(1) = 0.
  ## logL already initialized to 0 for these rows/"missing" label. This
  ## means a "missing" row contributes NO information to (beta, sigma,
  ## sigma_b): it never nudges the posterior over quadrature nodes, and
  ## its E-step "moments" (see pecme_truncated_moments()) are a model-
  ## implied prediction/imputation for that row, not evidence used to
  ## estimate the parameters that produced the prediction.

  li <- which(censor == "left")
  if (length(li)) {
    zL <- (matrix(left[li], length(li), K) - mu[li, , drop = FALSE]) / sigma
    logL[li, ] <- stats::pnorm(zL, log.p = TRUE)
  }

  ri <- which(censor == "right")
  if (length(ri)) {
    zR <- (matrix(right[ri], length(ri), K) - mu[ri, , drop = FALSE]) / sigma
    logL[ri, ] <- stats::pnorm(zR, lower.tail = FALSE, log.p = TRUE)
  }

  ii <- which(censor == "interval")
  if (length(ii)) {
    zL <- (matrix(left[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    zR <- (matrix(right[ii], length(ii), K) - mu[ii, , drop = FALSE]) / sigma
    logL[ii, ] <- matrix(.pecme_interval_logprob(as.vector(zL), as.vector(zR)), nrow = length(ii))
  }

  bad <- !(censor %in% c("none", "left", "right", "interval", "missing"))
  if (any(bad)) {
    stop("Unsupported censoring label(s): ",
         paste(unique(censor[bad]), collapse = ", "), call. = FALSE)
  }

  logL
}
