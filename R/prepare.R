## -----------------------------------------------------------------------
## Turn (formula, data, group, censor spec) into the numeric objects the
## engine needs: y, X (raw, not yet standardized), group (integer),
## censor/left/right. Unlike a plain `model.frame(..., na.action =
## na.omit)` pipeline, missing values are NOT silently dropped here:
## they are handed to the caller so pecme() can route them through the
## E-step ("missing" response) or the `missing_x` option (covariates).
## -----------------------------------------------------------------------

.pecme_has_intercept <- function(Terms) identical(as.integer(attr(Terms, "intercept")), 1L)

#' @keywords internal
pecme_prepare_data <- function(data, fixed_formula, random_var,
                                censor, left = NULL, right = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!inherits(fixed_formula, "formula")) stop("`fixed_formula` must be a formula.", call. = FALSE)
  ## `response` is always the left-hand side of `fixed_formula`; there
  ## is no separate `response` argument to keep in sync with it.
  ## Use the evaluated model-frame response. This is essential for
  ## transformed responses such as log(y), sqrt(y), I(y^2), etc.;
  ## all.vars() would return the underlying variable rather than the
  ## response actually specified on the left-hand side.
  response_expression <- paste(deparse(fixed_formula[[2L]]), collapse = "")


  resolve_vec <- function(x, name) {
    if (is.null(x)) return(rep(NA_real_, nrow(data)))
    if (is.character(x) && length(x) == 1L && x %in% names(data)) return(as.numeric(data[[x]]))
    if (length(x) != nrow(data)) stop("`", name, "` must have length nrow(data) or name a column of `data`.", call. = FALSE)
    as.numeric(x)
  }
  left_vec  <- resolve_vec(left,  "left")
  right_vec <- resolve_vec(right, "right")

  censor_labels <- c("none", "left", "right", "interval", "missing")
  if (is.character(censor) && length(censor) == 1L && censor %in% names(data)) {
    censor_vec <- as.character(data[[censor]])
  } else if (is.character(censor) && length(censor) == 1L && censor %in% censor_labels) {
    ## Convenience: a single recognized label (e.g. the default "none")
    ## applied to every row, rather than a column name.
    censor_vec <- rep(censor, nrow(data))
  } else if (length(censor) == nrow(data)) {
    censor_vec <- as.character(censor)
  } else {
    stop("`censor` must be \"none\"/\"left\"/\"right\"/\"interval\", name a column of `data`, ",
         "or be a vector of length nrow(data).", call. = FALSE)
  }

  needed <- unique(c(all.vars(fixed_formula), random_var))

  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  Terms <- stats::terms(fixed_formula)
  has_intercept <- .pecme_has_intercept(Terms)

  ## na.action = na.pass: keep every row so missingness can be handled
  ## explicitly rather than silently dropped by model.frame().
  mf <- stats::model.frame(fixed_formula, data = data, na.action = stats::na.pass)
  X_raw <- stats::model.matrix(fixed_formula, data = mf)

  if (nrow(X_raw) != nrow(data)) {
    stop("model.matrix() dropped rows unexpectedly (check factor levels in `data`).", call. = FALSE)
  }
  if (ncol(X_raw) == 0L) stop("The fixed-effects design matrix has zero columns.", call. = FALSE)

  y <- as.numeric(stats::model.response(mf))
  group_raw <- data[[random_var]]
  if (any(is.na(group_raw))) stop("`random_var` (grouping column) cannot contain NA.", call. = FALSE)
  group_levels <- sort(unique(as.character(group_raw)))
  group <- match(as.character(group_raw), group_levels)

  list(
    X = X_raw, y = y, group = group, group_levels = group_levels,
    censor = censor_vec, left = left_vec, right = right_vec,
    has_intercept = has_intercept, terms = Terms,
    response_expression = response_expression,
    xlevels = stats::.getXlevels(Terms, mf),
    contrasts = attr(X_raw, "contrasts")
  )
}
