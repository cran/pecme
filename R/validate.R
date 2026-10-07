## -----------------------------------------------------------------------
## Input validation
## -----------------------------------------------------------------------

#' @keywords internal
pecme_validate_inputs <- function(y, X, group, censor, left, right,
                                   penalty, lambda, alpha, gamma,
                                   quadrature, workers, allow_missing_response) {
  n <- length(y)
  if (!is.matrix(X) || nrow(X) != n) stop("X must be a matrix with nrow(X) == length(y).", call. = FALSE)
  if (length(group) != n) stop("length(group) must equal length(y).", call. = FALSE)
  if (length(censor) != n) stop("length(censor) must equal length(y).", call. = FALSE)
  if (length(left) != n || length(right) != n) stop("length(left)/length(right) must equal length(y).", call. = FALSE)

  bad_lab <- !(censor %in% c("none", "left", "right", "interval", "missing"))
  if (any(bad_lab)) {
    stop("`censor` must be one of \"none\",\"left\",\"right\",\"interval\",\"missing\"; ",
         "found: ", paste(unique(censor[bad_lab]), collapse = ", "), call. = FALSE)
  }

  if (any(is.na(group))) stop("`group` cannot contain NA.", call. = FALSE)
  if (any(!is.finite(X))) stop("X contains non-finite values (NA/NaN/Inf are not allowed here; ",
                                "use `missing_x` in pecme() for missing covariates).", call. = FALSE)

  if (!allow_missing_response) {
    if (any(censor == "none" & is.na(y))) {
      stop("`y` has NA where `censor == \"none\"`. Either relabel as censored/missing, ",
           "or set `missing_y = TRUE` in pecme() to estimate it via the E-step.", call. = FALSE)
    }
  }
  if (any(censor == "left"  & is.na(left)))  stop("`left` is NA for a \"left\"-censored row.", call. = FALSE)
  if (any(censor == "right" & is.na(right))) stop("`right` is NA for a \"right\"-censored row.", call. = FALSE)
  if (any(censor == "interval" & (is.na(left) | is.na(right)))) {
    stop("`left`/`right` must both be non-NA for \"interval\"-censored rows.", call. = FALSE)
  }
  if (any(censor == "interval" & left >= right, na.rm = TRUE)) {
    stop("`left` must be strictly less than `right` for \"interval\"-censored rows.", call. = FALSE)
  }

  penalty <- match.arg(penalty, .pecme_penalties)
  if (!is.finite(lambda) || lambda < 0) stop("`lambda` must be finite and >= 0.", call. = FALSE)
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("`alpha` must be in [0, 1].", call. = FALSE)

  if (penalty == "scad" && (!is.finite(gamma) || gamma <= 2)) {
    stop("`gamma` must be > 2 for the SCAD penalty (a common default is 3.7).", call. = FALSE)
  }
  if (penalty == "mcp" && (!is.finite(gamma) || gamma <= 1)) {
    stop("`gamma` must be > 1 for the MCP penalty (a common default is 3).", call. = FALSE)
  }

  if (!is.finite(quadrature) || quadrature < 5L) stop("`quadrature` must be an integer >= 5.", call. = FALSE)
  if (!is.finite(workers) || workers < 1L) stop("`workers` must be an integer >= 1.", call. = FALSE)

  invisible(TRUE)
}
