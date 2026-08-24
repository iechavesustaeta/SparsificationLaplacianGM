# Spectral-LCGGM and Spectral-HR pipelines with BIC selection over eta.

suppressPackageStartupMessages({
  library(parallel)
})
if (!exists("hr_golazo_cd"))    source("golazo_cd.R")
if (!exists("lcggm_mle_dense")) source("spectral_lcggm.R")
if (!exists("hr_mle_dense"))    source("spectral_hr.R")

spectral_lcggm_pipeline <- function(S, n, eta_grid,
                                    n_cores       = 1L,
                                    step1_tol     = 1e-5,
                                    step1_maxiter = 5000,
                                    step3_tol     = 1e-8,
                                    step3_maxiter = 5000,
                                    return_path   = FALSE) {
  mle <- tryCatch(
    lcggm_mle_dense(S, tol = step1_tol, maxiter = step1_maxiter),
    error = function(e) NULL
  )
  if (is.null(mle)) {
    return(list(K = NULL, param = NA_real_, bic = Inf,
                path = if (return_path) data.frame() else NULL))
  }

  eta_worker <- function(eta) {
    res <- tryCatch(
      spectral_lcggm_from_mle(mle, S, eta,
                              tol = step3_tol, maxiter = step3_maxiter),
      error = function(e) NULL
    )
    if (is.null(res) || is.null(res$K_refit)) {
      return(list(K = NULL, bic = Inf, eta = eta))
    }
    list(K = res$K_refit, bic = bic_lcggm(res$K_refit, S, n), eta = eta)
  }

  results <- if (n_cores > 1L && length(eta_grid) > 1L) {
    nc <- min(as.integer(n_cores), length(eta_grid))
    cl <- makeCluster(nc)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, source("experiment_helpers.R"))
    clusterExport(cl, varlist = c("mle", "S", "n",
                                  "step3_tol", "step3_maxiter"),
                  envir = environment())
    parLapplyLB(cl, eta_grid, eta_worker)
  } else {
    lapply(eta_grid, eta_worker)
  }

  bics <- vapply(results, function(r) r$bic, numeric(1))
  best <- which.min(bics)
  out <- list(K     = results[[best]]$K,
              param = results[[best]]$eta,
              bic   = results[[best]]$bic)
  if (return_path) {
    out$path <- data.frame(
      eta   = vapply(results, function(r) r$eta, numeric(1)),
      bic   = bics,
      edges = vapply(results, function(r)
                     if (is.null(r$K)) NA_integer_ else count_edges(r$K),
                     integer(1)),
      stringsAsFactors = FALSE
    )
  }
  out
}

spectral_hr_pipeline <- function(Gamma_hat, n, eta_grid,
                                 n_cores      = 1L,
                                 tol          = 1e-8,
                                 maxiter      = 3000,
                                 return_path  = FALSE) {
  mle <- tryCatch(hr_mle_dense(Gamma_hat, tol = tol, maxiter = maxiter),
                  error = function(e) NULL)
  if (is.null(mle)) {
    return(list(K = NULL, param = NA_real_, bic = Inf,
                path = if (return_path) data.frame() else NULL))
  }

  eta_worker <- function(eta) {
    res <- tryCatch(
      spectral_hr_from_mle(mle, eta, tol = tol, maxiter = maxiter),
      error = function(e) NULL
    )
    if (is.null(res) || is.null(res$Theta_refit)) {
      return(list(K = NULL, bic = Inf, eta = eta))
    }
    list(K = res$Theta_refit, bic = bic_hr(res$Theta_refit, Gamma_hat, n),
         eta = eta)
  }

  results <- if (n_cores > 1L && length(eta_grid) > 1L) {
    nc <- min(as.integer(n_cores), length(eta_grid))
    cl <- makeCluster(nc)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, source("experiment_helpers.R"))
    clusterExport(cl, varlist = c("mle", "Gamma_hat", "n",
                                  "tol", "maxiter"),
                  envir = environment())
    parLapplyLB(cl, eta_grid, eta_worker)
  } else {
    lapply(eta_grid, eta_worker)
  }

  bics <- vapply(results, function(r) r$bic, numeric(1))
  best <- which.min(bics)
  out <- list(K     = results[[best]]$K,
              param = results[[best]]$eta,
              bic   = results[[best]]$bic)
  if (return_path) {
    out$path <- data.frame(
      eta   = vapply(results, function(r) r$eta, numeric(1)),
      bic   = bics,
      edges = vapply(results, function(r)
                     if (is.null(r$K)) NA_integer_ else count_edges(r$K),
                     integer(1)),
      stringsAsFactors = FALSE
    )
  }
  out
}
