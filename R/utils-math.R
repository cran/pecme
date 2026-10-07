## -----------------------------------------------------------------------
## Low-level numerical utilities shared across the package.
## -----------------------------------------------------------------------

.pecme_GH_CACHE <- new.env(parent = emptyenv())

#' Gauss-Hermite nodes and weights for N(0, 1) integration
#'
#' Returns nodes \code{z} and weights \code{w} such that, for a smooth
#' function \code{f}, \eqn{E[f(Z)] \approx \sum_k w_k f(z_k)} when
#' \eqn{Z \sim N(0,1)}. Results are cached per quadrature order because
#' the Golub-Welsch eigendecomposition is by far the most expensive part
#' of the E-step if it is repeated every ECME iteration.
#'
#' @param n Integer number of quadrature nodes (>= 5).
#' @return A list with \code{nodes} and \code{weights}, both of length
#'   \code{n}, ordered by increasing node value, calibrated so that they
#'   directly integrate against the standard normal density (i.e. the
#'   physicists' Hermite rule has already been rescaled).
#' @keywords internal
gauss_hermite_normal <- function(n = 40L) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 5L) {
    stop("`n` must be a single integer >= 5.", call. = FALSE)
  }

  key <- as.character(n)
  if (exists(key, envir = .pecme_GH_CACHE, inherits = FALSE)) {
    return(get(key, envir = .pecme_GH_CACHE, inherits = FALSE))
  }

  ## Golub-Welsch construction of the physicists' Hermite rule:
  ##   integral exp(-x^2) f(x) dx  ~=  sum_k w_k f(x_k)
  k <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  off <- sqrt(k / 2)
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off

  ee <- eigen(J, symmetric = TRUE, only.values = FALSE)
  nodes_h <- ee$values
  weights_h <- sqrt(pi) * ee$vectors[1L, ]^2
  ord <- order(nodes_h)

  ## Rescale to integrate directly against the standard normal density:
  ## if Z ~ N(0,1) and x_k, w_k/sqrt(pi) are the (probabilists'-adjusted)
  ## Hermite nodes/weights, then Z = sqrt(2) x_k with probability
  ## w_k / sqrt(pi).
  nodes <- sqrt(2) * nodes_h[ord]
  weights <- weights_h[ord] / sqrt(pi)

  ans <- list(nodes = nodes, weights = weights)
  assign(key, ans, envir = .pecme_GH_CACHE)
  ans
}

#' Numerically stable log(sum(exp(x)))
#' @keywords internal
logsumexp <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' log(exp(log_a) - exp(log_b)) for log_a > log_b, stable near 0
#' @keywords internal
.pecme_logspace_sub <- function(log_a, log_b) {
  out <- rep(-Inf, length(log_a))
  ok <- is.finite(log_a) & is.finite(log_b) & (log_b < log_a)
  if (any(ok)) {
    out[ok] <- log_a[ok] + log1p(-exp(log_b[ok] - log_a[ok]))
  }
  zero_b <- is.finite(log_a) & is.infinite(log_b) & log_b < 0
  out[zero_b] <- log_a[zero_b]
  out
}

#' Stable log(Phi(b) - Phi(a)) for possibly-vectorized a < b
#' @keywords internal
.pecme_interval_logprob <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  n <- max(length(a), length(b))
  if (length(a) != n) a <- rep_len(a, n)
  if (length(b) != n) b <- rep_len(b, n)

  out <- rep(-Inf, n)
  ok <- is.finite(a) & is.finite(b) & (a < b)
  if (!any(ok)) return(out)

  aa <- a[ok]; bb <- b[ok]
  la <- stats::pnorm(aa, log.p = TRUE)
  lb <- stats::pnorm(bb, log.p = TRUE)
  direct <- .pecme_logspace_sub(lb, la)

  ## Right-tail identity: Phi(b) - Phi(a) = Phi(-a) - Phi(-b), evaluated
  ## with survival probabilities to keep precision in the far right tail.
  lsa <- stats::pnorm(aa, lower.tail = FALSE, log.p = TRUE)
  lsb <- stats::pnorm(bb, lower.tail = FALSE, log.p = TRUE)
  survival <- .pecme_logspace_sub(lsa, lsb)

  use_survival <- is.finite(survival) & (!is.finite(direct) | aa > 0)
  ans <- direct
  ans[use_survival] <- survival[use_survival]
  out[ok] <- ans
  out
}
