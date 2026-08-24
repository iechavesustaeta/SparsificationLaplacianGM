# Dense LCGGM maximum likelihood, BSS sparsification and refit on the retained support.

suppressPackageStartupMessages({
  library(MASS); library(igraph)
})
source("golazo_cd.R")

vario <- function(S) {
  d <- nrow(S); dS <- diag(S)
  G <- outer(dS, rep(1, d)) + outer(rep(1, d), dS) - 2 * S
  diag(G) <- 0; G
}
logdet_lap <- function(Th) as.numeric(determinant(Th[-1, -1], log = TRUE)$modulus)
lcggm_loglik <- function(K, S) logdet_lap(K) - sum(K * S)
count_edges  <- function(K, tol = 1e-3) sum(abs(K[upper.tri(K)]) > tol)

bss_skeleton <- function(L_B, eta = 2) {
  p_ <- nrow(L_B)
  int_idx <- which(upper.tri(L_B) & L_B < -1e-8, arr.ind = TRUE)
  m  <- nrow(int_idx); if (m == 0L) return(NULL)
  w_int <- -L_B[int_idx]

  eig <- eigen(L_B, symmetric = TRUE)
  ev  <- eig$values; tol_e <- 1e-10 * max(abs(ev))
  r   <- sum(ev > tol_e); if (r == 0L) return(NULL)
  Q   <- eig$vectors[, seq_len(r), drop = FALSE]; Lam <- ev[seq_len(r)]

  QT_b <- t(Q[int_idx[, 1], , drop = FALSE]) - t(Q[int_idx[, 2], , drop = FALSE])
  V    <- sweep(QT_b / sqrt(Lam), 2, sqrt(w_int), `*`)

  sqd <- sqrt(eta); delta_L <- 1; delta_U <- (sqd + 1) / (sqd - 1)
  ell <- -sqd * r;  u <- sqd * r * delta_U; q_iter <- ceiling(eta * r)
  A <- matrix(0, r, r); s <- rep(0, m); Ir <- diag(r)
  Phi_u_curr <- r / u; Phi_l_curr <- r / (-ell)

  for (t in seq_len(q_iter)) {
    u_new <- u + delta_U; ell_new <- ell + delta_L
    M_U <- tryCatch(solve(u_new * Ir - A), error = function(e) NULL)
    M_L <- tryCatch(solve(A - ell_new * Ir), error = function(e) NULL)
    if (is.null(M_U) || is.null(M_L)) break
    Phi_u_new <- sum(diag(M_U)); Phi_l_new <- sum(diag(M_L))
    dPhi_u <- Phi_u_curr - Phi_u_new; dPhi_l <- Phi_l_new - Phi_l_curr
    if (dPhi_u <= 0 || dPhi_l <= 0) break
    MUV <- M_U %*% V; MU2V <- M_U %*% MUV
    MLV <- M_L %*% V; ML2V <- M_L %*% MLV
    U_star <- colSums(V * MU2V) / dPhi_u + colSums(V * MUV)
    L_star <- colSums(V * ML2V) / dPhi_l - colSums(V * MLV)
    slack  <- L_star - U_star; feas <- which(slack > 0)
    if (length(feas) == 0L) break
    e_pick <- feas[which.max(slack[feas])]
    alpha  <- 2 / (U_star[e_pick] + L_star[e_pick])
    v <- V[, e_pick]; A <- A + alpha * tcrossprod(v); s[e_pick] <- s[e_pick] + alpha
    ell <- ell_new; u <- u_new
    Phi_u_curr <- tryCatch(sum(diag(solve(u * Ir - A))), error = function(e) Inf)
    Phi_l_curr <- tryCatch(sum(diag(solve(A - ell * Ir))), error = function(e) Inf)
    if (!is.finite(Phi_u_curr) || !is.finite(Phi_l_curr)) break
  }
  keep <- which(s > 0)
  if (length(keep) == 0L) return(NULL)
  list(keep = keep, int_idx = int_idx, s = s, p = p_)
}

mle_refit_lcggm <- function(G_S, keep, int_idx, p_, tol = 1e-9, maxiter = 5000,
                            w_init = NULL) {
  Lp <- matrix(-Inf, p_, p_); Up <- matrix(Inf, p_, p_)
  diag(Lp) <- 0; diag(Up) <- 0
  for (k in keep) {
    i <- int_idx[k, 1]; j <- int_idx[k, 2]
    Lp[i, j] <- 0; Lp[j, i] <- 0
  }
  res <- hr_golazo_cd(G_S, lambda = 1, Lpen = Lp, Upen = Up,
                      w_init = w_init,
                      tol = tol, maxiter = maxiter)
  res$Theta
}

bss_laplacian <- function(L_B, sk) {
  p_ <- sk$p
  out <- matrix(0, p_, p_)
  for (k in sk$keep) {
    i <- sk$int_idx[k, 1]; j <- sk$int_idx[k, 2]
    w <- sk$s[k] * (-L_B[i, j])
    out[i, i] <- out[i, i] + w; out[j, j] <- out[j, j] + w
    out[i, j] <- out[i, j] - w; out[j, i] <- out[j, i] - w
  }
  out
}

lcggm_mle_dense <- function(S, tol = 1e-9, maxiter = 5000) {
  d <- nrow(S); G_S <- vario(S)
  Lp <- matrix(0, d, d); Up <- matrix(Inf, d, d)
  diag(Lp) <- diag(Up) <- 0
  r <- hr_golazo_cd(G_S, lambda = 1, Lpen = Lp, Upen = Up,
                    tol = tol, maxiter = maxiter)
  list(K_hat = r$Theta, G_S = G_S,
       converged = r$converged, iters = r$iters)
}

spectral_lcggm_from_mle <- function(mle, S, eta,
                                    tol = 1e-9, maxiter = 5000,
                                    verbose = FALSE) {
  K_hat <- mle$K_hat; G_S <- mle$G_S; d <- nrow(K_hat)

  if (verbose) cat("[Step 2] BSS sparsification, eta =", eta, "...\n")
  sk <- bss_skeleton(K_hat, eta = eta)
  if (is.null(sk)) {
    warning("BSS returned NULL (degenerate Laplacian or no candidate edges)")
    return(list(K_hat = K_hat, K_tilde = NULL, K_refit = NULL,
                E_tilde = NULL, eta = eta, sk = NULL))
  }
  K_tilde <- bss_laplacian(K_hat, sk)
  g_t <- graph_from_adjacency_matrix(abs(K_tilde) > 1e-10,
                                     mode = "undirected", diag = FALSE)
  if (!is_connected(g_t)) {
    warning("BSS skeleton is disconnected; refit may fail")
  }

  edge_mask <- matrix(FALSE, d, d)
  for (k in sk$keep) {
    i <- sk$int_idx[k, 1]; j <- sk$int_idx[k, 2]
    edge_mask[i, j] <- edge_mask[j, i] <- TRUE
  }

  if (verbose) cat("[Step 3] MLE refit on", length(sk$keep), "edges...\n")
  K_init <- K_hat
  K_init[!edge_mask & (row(K_init) != col(K_init))] <- 0
  diag(K_init) <- 0
  diag(K_init) <- -rowSums(K_init)
  K_refit <- mle_refit_lcggm(G_S, sk$keep, sk$int_idx, sk$p,
                             tol = tol, maxiter = maxiter,
                             w_init = K_init)
  list(
    K_hat   = K_hat,
    K_tilde = K_tilde,
    K_refit = K_refit,
    E_tilde = edge_mask,
    eta = eta, sk = sk,
    n_edges_full  = sum(abs(K_hat[upper.tri(K_hat)])     > 1e-8),
    n_edges_bss   = length(sk$keep),
    n_edges_refit = sum(abs(K_refit[upper.tri(K_refit)]) > 1e-8),
    loglik_full   = lcggm_loglik(K_hat,   S),
    loglik_bss    = lcggm_loglik(K_tilde, S),
    loglik_refit  = lcggm_loglik(K_refit, S)
  )
}

spectral_lcggm <- function(S, eta, tol = 1e-9, maxiter = 5000, verbose = FALSE) {
  if (verbose) cat("[Step 1] dense LCGGM-MLE...\n")
  mle <- lcggm_mle_dense(S, tol = tol, maxiter = maxiter)
  spectral_lcggm_from_mle(mle, S, eta, tol = tol, maxiter = maxiter,
                          verbose = verbose)
}

if (sys.nframe() == 0) {
  set.seed(2024)
  d <- 20; n <- 600
  g <- igraph::sample_pa(d, m = 2, directed = FALSE)
  A <- as.matrix(igraph::as_adjacency_matrix(g))
  W <- A * matrix(runif(d * d, 0.5, 1.5), d, d); W <- (W + t(W)) / 2
  K_star <- diag(rowSums(W)) - W
  E_star <- (abs(K_star) > 1e-10) & (row(K_star) != col(K_star))
  cat(sprintf("Truth: d=%d, n=%d, edges=%d\n",
              d, n, sum(E_star[upper.tri(E_star)])))

  Kp  <- ginv(K_star)
  P   <- diag(d) - matrix(1, d, d) / d
  Kpp <- P %*% Kp %*% P
  ev  <- eigen(Kpp, symmetric = TRUE)
  pos <- ev$values > 1e-10
  U_  <- ev$vectors[, pos, drop = FALSE]
  D_  <- diag(sqrt(ev$values[pos]), nrow = sum(pos))
  X   <- matrix(rnorm(n * sum(pos)), n, sum(pos)) %*% t(U_ %*% D_)
  S   <- crossprod(X) / n

  metrics <- function(K, label) {
    Edge <- (abs(K) > 1e-3) & (row(K) != col(K))
    tp <- sum(Edge & E_star) / 2
    fp <- sum(Edge & !E_star & row(Edge) != col(Edge)) / 2
    fn <- sum(!Edge & E_star) / 2
    prec <- if (tp + fp > 0) tp / (tp + fp) else NA
    rec  <- if (tp + fn > 0) tp / (tp + fn) else NA
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
              2 * prec * rec / (prec + rec) else NA
    cat(sprintf("  %-20s edges=%3d  TP=%3d FP=%3d FN=%3d  prec=%.3f rec=%.3f F1=%.3f  ||K-K*||F=%.3f\n",
                label, sum(Edge) / 2, tp, fp, fn, prec, rec, f1,
                norm(K - K_star, "F")))
  }

  for (eta in c(1.5, 2, 3, 5, 8)) {
    cat(sprintf("\n--- eta = %.1f ---\n", eta))
    res <- spectral_lcggm(S, eta = eta, verbose = TRUE)
    if (is.null(res$K_refit)) { cat("  pipeline failed\n"); next }
    metrics(res$K_hat,   "K_hat (dense)")
    metrics(res$K_tilde, "K_tilde (BSS)")
    metrics(res$K_refit, "K_refit (Spectral-MLE)")
    cat(sprintf("  log-lik full / BSS / refit = %.3f / %.3f / %.3f\n",
                res$loglik_full, res$loglik_bss, res$loglik_refit))

    K <- res$K_refit
    extra <- ((abs(K) > 1e-8) & (row(K) != col(K))) & !res$E_tilde
    cat(sprintf("  support sanity: refit has %d edges outside BSS support (should be 0)\n",
                sum(extra) / 2))
    cat(sprintf("  Laplacian sanity: max |row sum| = %.2e, max off-diag = %.2e (<= 0)\n",
                max(abs(rowSums(K))),
                max(K[upper.tri(K)])))
  }
}
