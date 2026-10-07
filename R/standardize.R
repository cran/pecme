## -----------------------------------------------------------------------
## Standardize the fixed-effects design matrix (mean 0 / sd 1 per
## column), leaving the intercept column (if present) untouched.
## -----------------------------------------------------------------------

#' @keywords internal
pecme_standardize <- function(X, has_intercept = TRUE) {
  X <- as.matrix(X)
  p <- ncol(X)
  center <- rep(0, p)
  scale  <- rep(1, p)
  start <- if (has_intercept) 2L else 1L

  if (start <= p) {
    for (j in start:p) {
      m <- mean(X[, j])
      s <- stats::sd(X[, j])
      center[j] <- m
      if (is.na(s) || s <= .Machine$double.eps) {
        scale[j] <- 1
      } else {
        scale[j] <- s
        X[, j] <- (X[, j] - m) / s
      }
    }
  }

  list(X = X, center = center, scale = scale)
}
