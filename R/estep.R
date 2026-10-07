## -----------------------------------------------------------------------
## E-step for the censored/missing-response random-intercept LMM
##
##   Y*_ij = x_ij'beta + b_j + e_ij,   e_ij ~ N(0, sigma^2)
##   b_j   ~ N(0, sigma_b^2)
##
## The random intercept is integrated out NUMERICALLY via Gauss-Hermite
## quadrature (cached, see utils-math.R) -- a fixed-order approximation
## to the continuous integral, not an exact evaluation of it (see
## pecme_quadrature_check() for a sensitivity diagnostic). Conditional
## on that discretization, the posterior weights across the quadrature
## nodes are themselves computed exactly. Conditional on a quadrature
## node b_j = z_k, the observations within group j are independent, so
## the truncated/missing-response moments in truncation.R can be
## evaluated once per (observation, node) pair and then averaged over
## those posterior weights.
## -----------------------------------------------------------------------

#' Sequential (single-core) E-step
#'
#' @param y,X,beta,group,sigma,sigma_b,censor,left,right Model
#'   quantities; see [pecme()].
#' @param quadrature Number of Gauss-Hermite nodes for the random
#'   intercept.
#' @return A list of posterior moments: `Ey`, `Ey2` (per observation),
#'   `Eb`, `Ebb` (per observation, broadcast from the group), `Eyb`
#'   (per observation), and `group_Eb`, `group_Ebb` (per group, in the
#'   order of `sort(unique(group))`).
#' @keywords internal
pecme_e_step <- function(y, X, beta, group, sigma, sigma_b,
                          censor, left = NULL, right = NULL,
                          quadrature = 40L) {
  X <- as.matrix(X); beta <- as.numeric(beta)
  group <- as.integer(group); censor <- as.character(censor)
  n <- nrow(X)

  if (length(y) != n || length(group) != n || length(censor) != n) {
    stop("y, group, and censor must have length nrow(X).", call. = FALSE)
  }
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(sigma_b) || sigma_b <= 0) {
    stop("sigma and sigma_b must be positive and finite.", call. = FALSE)
  }
  if (is.null(left))  left  <- rep(NA_real_, n)
  if (is.null(right)) right <- rep(NA_real_, n)

  gh <- gauss_hermite_normal(quadrature)
  b_nodes <- sigma_b * gh$nodes
  log_w <- log(gh$weights)

  eta <- as.vector(X %*% beta)

  Ey <- Ey2 <- Eb <- Ebb <- Eyb <- numeric(n)
  group_ids <- sort(unique(group))
  J <- length(group_ids)
  group_Eb <- group_Ebb <- numeric(J)

  for (gidx in seq_along(group_ids)) {
    idx <- which(group == group_ids[gidx])
    mu <- outer(eta[idx], b_nodes, "+")           # obs x nodes

    logL <- pecme_conditional_loglik(
      y = y[idx], mu = mu, sigma = sigma,
      censor = censor[idx], left = left[idx], right = right[idx]
    )

    log_post_group <- colSums(logL) + log_w
    log_norm <- logsumexp(log_post_group)
    if (!is.finite(log_norm)) {
      stop("E-step produced a non-finite posterior normalization for group ",
           group_ids[gidx], ". Check censoring limits and variance values.",
           call. = FALSE)
    }
    post_w <- exp(log_post_group - log_norm)

    Eb_g <- sum(post_w * b_nodes)
    Ebb_g <- sum(post_w * b_nodes^2)
    group_Eb[gidx] <- Eb_g
    group_Ebb[gidx] <- Ebb_g

    cm <- pecme_truncated_moments(
      mu = mu, sigma = sigma, censor = censor[idx],
      left = left[idx], right = right[idx], y = y[idx]
    )

    Ey[idx]  <- as.vector(cm$mean %*% post_w)
    Ey2[idx] <- as.vector(cm$second %*% post_w)
    Eb[idx]  <- Eb_g
    Ebb[idx] <- Ebb_g
    Eyb[idx] <- as.vector(sweep(cm$mean, 2L, b_nodes, "*") %*% post_w)
  }

  list(Ey = Ey, Ey2 = Ey2, Eb = Eb, Ebb = Ebb, Eyb = Eyb,
       group_Eb = group_Eb, group_Ebb = group_Ebb,
       group_ids = group_ids, quadrature = quadrature)
}

#' E-step dispatcher: sequential or cluster-parallel over groups
#'
#' With `workers <= 1`, or too few groups to make parallelism
#' worthwhile, this simply calls [pecme_e_step()]. Otherwise the groups
#' are partitioned across a pre-existing PSOCK/FORK cluster (created
#' once per [pecme()] call and reused across iterations).
#'
#' @param cluster A `"cluster"` object from `parallel::makeCluster()`,
#'   or `NULL` for sequential execution.
#' @keywords internal
pecme_e_step_dispatch <- function(y, X, beta, group, sigma, sigma_b,
                                   censor, left = NULL, right = NULL,
                                   quadrature = 40L, workers = 1L,
                                   cluster = NULL) {
  n <- length(y)
  group_levels <- sort(unique(group))
  J <- length(group_levels)

  workers <- as.integer(workers)
  if (!is.finite(workers) || workers < 1L) workers <- 1L
  workers <- min(workers, J)

  if (workers <= 1L || J <= 2L || is.null(cluster)) {
    return(pecme_e_step(y, X, beta, group, sigma, sigma_b,
                         censor, left, right, quadrature))
  }

  cluster_size <- length(cluster)
  workers <- min(workers, cluster_size, J)
  if (workers <= 1L) {
    return(pecme_e_step(y, X, beta, group, sigma, sigma_b,
                         censor, left, right, quadrature))
  }

  idx_list <- lapply(group_levels, function(g) which(group == g))

  worker_fun <- function(ii) {
    res <- pecme_e_step(
      y = y[ii], X = X[ii, , drop = FALSE], beta = beta,
      group = rep(1L, length(ii)), sigma = sigma, sigma_b = sigma_b,
      censor = censor[ii],
      left  = if (is.null(left))  NULL else left[ii],
      right = if (is.null(right)) NULL else right[ii],
      quadrature = quadrature
    )
    list(index = ii, result = res)
  }

  parallel::clusterExport(
    cluster,
    varlist = c("y", "X", "beta", "sigma", "sigma_b", "censor", "left",
                "right", "quadrature", "worker_fun"),
    envir = environment()
  )

  pieces <- parallel::parLapplyLB(cluster, idx_list, worker_fun)

  Ey <- Eb <- Ey2 <- Ebb_obs <- Eyb <- numeric(n)
  group_Eb <- group_Ebb <- numeric(J)

  for (j in seq_along(pieces)) {
    ii <- pieces[[j]]$index
    r  <- pieces[[j]]$result
    Ey[ii] <- r$Ey; Eb[ii] <- r$Eb
    Ey2[ii] <- r$Ey2; Ebb_obs[ii] <- r$Ebb; Eyb[ii] <- r$Eyb
    group_Eb[j]  <- as.numeric(r$group_Eb[1L])
    group_Ebb[j] <- as.numeric(r$group_Ebb[1L])
  }

  list(Ey = Ey, Ey2 = Ey2, Eb = Eb, Ebb = Ebb_obs, Eyb = Eyb,
       group_Eb = group_Eb, group_Ebb = group_Ebb,
       group_ids = group_levels, quadrature = quadrature)
}
