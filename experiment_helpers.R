# Shared solvers, baseline wrappers, evaluation metrics and tuning grids.

suppressPackageStartupMessages({
  library(MASS); library(igraph)
  library(glassoFast); library(graphicalExtremes)
  if (requireNamespace("spectralGraphTopology", quietly = TRUE))
    library(spectralGraphTopology)
})
source("golazo_cd.R")
source("spectral_lcggm.R")
source("spectral_hr.R")
source("spectral_pipeline.R")

edge_tol <- 1e-3

count_edges <- function(K, tol = edge_tol)
  if (is.null(K)) NA_integer_ else sum(abs(K[upper.tri(K)]) > tol)

edge_mask_of <- function(K, tol = edge_tol) {
  if (is.null(K)) return(NULL)
  m <- (abs(K) > tol) & (row(K) != col(K))
  m
}

edge_metrics <- function(K, E_true, tol = edge_tol) {
  if (is.null(K)) return(list(precision = NA, recall = NA, f1 = NA, n_edges = NA))
  est <- edge_mask_of(K, tol)
  TP <- sum(est & E_true) / 2
  FP <- sum(est & !E_true) / 2
  FN <- sum(!est & E_true) / 2
  prec <- if (TP + FP > 0) TP / (TP + FP) else 0
  rec  <- if (TP + FN > 0) TP / (TP + FN) else 0
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  list(precision = prec, recall = rec, f1 = f1, n_edges = TP + FP)
}

loglik_lcggm <- function(K, S_test) {
  if (is.null(K)) return(NA_real_)
  d <- nrow(K); P <- diag(d) - matrix(1, d, d) / d
  K_p <- P %*% K %*% P; S_p <- P %*% S_test %*% P
  ev  <- eigen(K_p, symmetric = TRUE, only.values = TRUE)$values
  evp <- ev[ev > 1e-10]
  if (length(evp) < d - 1) return(NA_real_)
  sum(log(evp)) - sum(K_p * S_p)
}

bic_lcggm <- function(K, S, n) {
  if (is.null(K)) return(Inf)
  ll <- loglik_lcggm(K, S)
  if (is.na(ll)) return(Inf)
  -n * ll + count_edges(K) * log(n)
}

loglik_hr <- function(Theta, Gamma_test) {
  if (is.null(Theta)) return(NA_real_)
  d <- nrow(Theta)
  P <- diag(d) - matrix(1, d, d) / d
  Th_p <- P %*% Theta %*% P
  ev  <- eigen(Th_p, symmetric = TRUE, only.values = TRUE)$values
  evp <- ev[ev > 1e-10]
  if (length(evp) < d - 1) return(NA_real_)
  sum(log(evp)) + 0.5 * sum(Theta * Gamma_test)
}

bic_hr <- function(Theta, Gamma, n) {
  if (is.null(Theta)) return(Inf)
  ll <- loglik_hr(Theta, Gamma)
  if (is.na(ll)) return(Inf)
  -n * ll + count_edges(Theta) * log(n)
}

fit_glasso <- function(S, lambda) {
  res <- tryCatch(glassoFast(S, rho = lambda, thr = 1e-7, maxIt = 10000),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  K <- res$wi; (K + t(K)) / 2
}

fit_cgl <- function(S, alpha = 0, tol = 1e-9, maxiter = 5000) {
  mle <- tryCatch(lcggm_mle_dense(S, tol = tol, maxiter = maxiter),
                  error = function(e) NULL)
  if (is.null(mle)) NULL else mle$K_hat
}

fit_gleadmm <- function(S, alpha, maxiter = 10000) {
  tryCatch(learn_laplacian_gle_admm(
             S, alpha = alpha, reltol = 1e-5,
             maxiter = maxiter, verbose = FALSE)$laplacian,
           error = function(e) NULL)
}

if (!exists("L", mode = "function") &&
    requireNamespace("spectralGraphTopology", quietly = TRUE))
  L <- spectralGraphTopology::L

fit_sgl_scad <- function(S, alpha, gamma = 2.001, maxiter = 2000) {
  d <- nrow(S)
  w0 <- rep(1 / d, d * (d - 1) / 2)
  res <- tryCatch(sparseGraph::learn_laplacian_pgd_connected(
                    S, w0 = w0, alpha = alpha,
                    sparsity_type = "scad", gamma = gamma,
                    backtrack = TRUE,
                    maxiter = maxiter, reltol = 1e-5,
                    verbose = FALSE),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  if (!is.null(res$laplacian))      res$laplacian
  else if (!is.null(res$Laplacian)) res$Laplacian
  else NULL
}

fit_emtp2 <- function(Gamma_hat, tol = 1e-8, maxiter = 3000) {
  mle <- tryCatch(hr_mle_dense(Gamma_hat, tol = tol, maxiter = maxiter),
                  error = function(e) NULL)
  if (is.null(mle)) NULL else mle$Theta_hat
}

fit_eglearn_pkg <- function(X, p, rholist, n_eff = NULL) {
  if (is.null(n_eff)) n_eff <- nrow(X)
  fit <- tryCatch(graphicalExtremes::eglearn(
                    data = X, p = p, rholist = rholist,
                    reg_method = "ns", complete_Gamma = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(Theta = NULL, rho = NA_real_))

  best_bic <- Inf; best_Th <- NULL; best_rho <- NA_real_
  Gamma_tr <- graphicalExtremes::emp_vario(data = X, p = p)
  for (j in seq_along(rholist)) {
    G_j <- fit$Gamma[[j]]
    if (is.null(G_j) || anyNA(G_j)) next
    Th <- tryCatch(graphicalExtremes::Gamma2Theta(G_j),
                   error = function(e) NULL)
    if (is.null(Th)) next
    b <- bic_hr(Th, Gamma_tr, n_eff)
    if (b < best_bic) { best_bic <- b; best_Th <- Th; best_rho <- rholist[j] }
  }
  list(Theta = best_Th, rho = best_rho)
}

bic_select_alpha <- function(S, n, grid, fit_fn,
                             scorer = bic_lcggm) {
  best_bic <- Inf; best_K <- NULL; best_param <- NA_real_
  for (param in grid) {
    K <- tryCatch(fit_fn(S, param), error = function(e) NULL)
    if (is.null(K)) next
    b <- scorer(K, S, n)
    if (b < best_bic) { best_bic <- b; best_K <- K; best_param <- param }
  }
  list(K = best_K, param = best_param, bic = best_bic)
}

bic_select_spectral_lcggm <- function(S, n, eta_grid,
                                      n_cores       = 1L,
                                      step1_tol     = 1e-5,
                                      step1_maxiter = 5000,
                                      step3_tol     = 1e-8,
                                      step3_maxiter = 5000) {
  spectral_lcggm_pipeline(S, n, eta_grid,
                          n_cores       = n_cores,
                          step1_tol     = step1_tol,
                          step1_maxiter = step1_maxiter,
                          step3_tol     = step3_tol,
                          step3_maxiter = step3_maxiter)
}

bic_select_spectral_hr <- function(Gamma_hat, n, eta_grid,
                                   n_cores = 1L,
                                   tol     = 1e-8,
                                   maxiter = 3000) {
  spectral_hr_pipeline(Gamma_hat, n, eta_grid,
                       n_cores = n_cores,
                       tol     = tol,
                       maxiter = maxiter)
}

alpha_grid_default  <- 10^seq(-4, 0, length.out = 8)
eta_grid_default    <- c(1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
eta_grid_lcggm_sim  <- c(1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
                         10.0, 12.0, 16.0, 20.0, 26.0, 32.0)
glasso_grid_default <- 10^seq(-3, 0, length.out = 8)
rho_grid_default    <- seq(0.05, 0.50, length.out = 12)

method_labels_lcggm <- c(
  glasso   = "GLasso",
  cgl      = "CGL (Egilmez 2017)",
  gleadmm  = "GLE-ADMM (Ying 2020)",
  sgl_scad = "NGL-SCAD (Ying 2020 NeurIPS)",
  spectral = "Spectral-LCGGM (ours)"
)
method_order_lcggm <- names(method_labels_lcggm)

method_labels_hr <- c(
  emtp2    = "EMTP2",
  eglearn  = "eglearn-NS (Engelke et al)",
  spectral = "Spectral-HR (ours)"
)
method_order_hr <- names(method_labels_hr)
