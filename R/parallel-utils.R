## -----------------------------------------------------------------------
## Parallel cluster helpers
##
## A single cluster is created once per pecme()/grid_pecme() call and
## reused across ECME iterations (E-step) and/or lambda values (grid
## search), rather than spun up and torn down every iteration -- that
## per-iteration overhead is what usually makes naive parallel EM
## slower than sequential EM on typical cluster/longitudinal data set
## sizes.
## -----------------------------------------------------------------------

#' @keywords internal
pecme_make_cluster <- function(workers) {
  workers <- as.integer(workers)
  if (!is.finite(workers) || workers <= 1L) return(NULL)
  avail <- parallel::detectCores(logical = TRUE)
  if (is.na(avail)) avail <- workers
  workers <- max(1L, min(workers, avail))
  if (workers <= 1L) return(NULL)

  cl <- parallel::makeCluster(workers, type = "PSOCK")

  ## Export every object of the pecme namespace (not just the exported
  ## API) so worker processes can call internal helpers such as
  ## pecme_e_step() regardless of whether pecme was loaded with
  ## library() or devtools::load_all(). Done once at cluster creation,
  ## not per iteration.
  ns <- tryCatch(asNamespace("pecme"), error = function(e) NULL)
  if (!is.null(ns)) {
    parallel::clusterExport(cl, varlist = ls(ns, all.names = TRUE), envir = ns)
  }
  cl
}

#' @keywords internal
pecme_stop_cluster <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  invisible(NULL)
}
