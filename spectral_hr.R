# Husler-Reiss surrogate maximum likelihood and its spectral sparsification.

suppressPackageStartupMessages({
  library(MASS); library(igraph)
})
source("golazo_cd.R")
source("spectral_lcggm.R")

hr_mle_dense <- function(Gamma_hat, tol = 1e-9, maxiter = 5000) {
  d <- nrow(Gamma_hat)
  Lp <- matrix(0,   d, d); Up <- matrix(Inf, d, d)
  diag(Lp) <- diag(Up) <- 0
  r <- hr_golazo_cd(Gamma_hat, lambda = 1, Lpen = Lp, Upen = Up,
                    tol = tol, maxiter = maxiter)
  list(Theta_hat = r$Theta, Gamma_hat = Gamma_hat,
       converged = r$converged, iters = r$iters)
}

hr_mle_refit <- function(Gamma_hat, edge_mask,
                         tol = 1e-9, maxiter = 5000,
                         w_init = NULL) {
  d <- nrow(Gamma_hat)
  Lp <- matrix(-Inf, d, d); Up <- matrix(Inf, d, d)
  diag(Lp) <- diag(Up) <- 0
  off <- (row(edge_mask) != col(edge_mask))
  Lp[edge_mask & off] <- 0
  r <- hr_golazo_cd(Gamma_hat, lambda = 1, Lpen = Lp, Upen = Up,
                    w_init = w_init,
                    tol = tol, maxiter = maxiter)
  r$Theta
}

spectral_hr_from_mle <- function(mle, eta,
                                 tol = 1e-9, maxiter = 5000,
                                 verbose = FALSE) {
  Theta_hat <- mle$Theta_hat; Gamma_hat <- mle$Gamma_hat
  d <- nrow(Theta_hat)

  if (verbose) cat("[Step 2] BSS on Theta_hat, eta =", eta, "...\n")
  sk <- bss_skeleton(Theta_hat, eta = eta)
  if (is.null(sk))
    return(list(Theta_hat = Theta_hat, Theta_tilde = NULL,
                Theta_refit = NULL, E_tilde = NULL, eta = eta, sk = NULL))

  Theta_tilde <- bss_laplacian(Theta_hat, sk)
  edge_mask <- matrix(FALSE, d, d)
  for (k in sk$keep) {
    i <- sk$int_idx[k, 1]; j <- sk$int_idx[k, 2]
    edge_mask[i, j] <- edge_mask[j, i] <- TRUE
  }

  if (verbose) cat("[Step 3] HR-MLE refit on", length(sk$keep), "edges...\n")
  Theta_init <- Theta_hat
  Theta_init[!edge_mask & (row(Theta_init) != col(Theta_init))] <- 0
  diag(Theta_init) <- 0
  diag(Theta_init) <- -rowSums(Theta_init)
  Theta_refit <- hr_mle_refit(Gamma_hat, edge_mask,
                              tol = tol, maxiter = maxiter,
                              w_init = Theta_init)

  list(Theta_hat = Theta_hat, Theta_tilde = Theta_tilde,
       Theta_refit = Theta_refit, E_tilde = edge_mask,
       eta = eta, sk = sk)
}

spectral_hr <- function(Gamma_hat, eta,
                        tol = 1e-9, maxiter = 5000, verbose = FALSE) {
  if (verbose) cat("[Step 1] dense HR-MLE...\n")
  mle <- hr_mle_dense(Gamma_hat, tol = tol, maxiter = maxiter)
  spectral_hr_from_mle(mle, eta, tol = tol, maxiter = maxiter, verbose = verbose)
}
