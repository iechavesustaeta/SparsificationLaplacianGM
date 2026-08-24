# Coordinate-descent solvers for GOLAZO-penalised Gaussian and Laplacian models.

library(MASS)

soft_shrink_LU <- function(Y, lambda, L, U) {
  pmin(Y - lambda * L, 0) + pmax(Y - lambda * U, 0)
}

.sym_pinv <- function(M, tol = NULL) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  vals <- ev$values
  if (is.null(tol)) tol <- length(vals) * .Machine$double.eps * max(abs(vals))
  inv_vals <- ifelse(abs(vals) > tol, 1 / vals, 0)
  ev$vectors %*% (inv_vals * t(ev$vectors))
}

.w_update <- function(W, i, j, delta) {
  wi <- W[, i]; wj <- W[, j]
  a  <- 1 + delta * W[i, j]
  D  <- a^2 - delta^2 * W[i, i] * W[j, j]
  if (abs(D) < 1e-14) return(NULL)
  W - (delta / D) * (a * (outer(wi, wj) + outer(wj, wi)) -
                     delta * W[j, j] * outer(wi, wi) -
                     delta * W[i, i] * outer(wj, wj))
}

.tplus_update <- function(Tplus, v, delta) {
  tv   <- as.numeric(Tplus %*% v)
  alpha <- as.numeric(t(v) %*% tv)
  denom <- 1 + delta * alpha
  if (abs(denom) < 1e-10) return(NULL)
  Tplus - (delta / denom) * outer(tv, tv)
}

.tplus_update_ij <- function(Tplus, i, j, delta, tv = NULL, alpha = NULL) {
  if (is.null(tv))    tv    <- Tplus[, i] - Tplus[, j]
  if (is.null(alpha)) alpha <- tv[i] - tv[j]
  denom <- 1 + delta * alpha
  if (abs(denom) < 1e-10) return(NULL)
  Tplus - (delta / denom) * tcrossprod(tv)
}

gauss_golazo_cd <- function(S, lambda, Lpen, Upen,
                             tol = 1e-8, maxiter = 1000, verbose = FALSE) {
  p <- nrow(S)
  K <- diag(p)
  W <- diag(p)
  converged <- FALSE

  for (iter in seq_len(maxiter)) {
    max_chg <- 0

    for (i in seq_len(p)) {
      g_diag <- W[i, i] - S[i, i]
      c_diag <- W[i, i]^2
      if (c_diag > 1e-14) {
        delta_d <- g_diag / c_diag
        if (abs(delta_d) >= 1e-16) {
          K[i, i] <- K[i, i] + delta_d
          wi <- W[, i]
          denom_d <- 1 + delta_d * W[i, i]
          if (abs(denom_d) > 1e-14) {
            W <- W - (delta_d / denom_d) * outer(wi, wi)
            W <- (W + t(W)) / 2
          } else {
            W <- solve(K)
          }
          max_chg <- max(max_chg, abs(delta_d))
        }
      }

      for (j in seq_len(i - 1L)) {
        g <- W[i, j] - S[i, j]
        c <- W[i, i] * W[j, j] + W[i, j]^2
        if (c < 1e-14) next

        t_old <- K[i, j]
        t_new <- soft_shrink_LU(t_old + g / c, lambda / c, Lpen[i, j], Upen[i, j])
        delta  <- t_new - t_old
        if (abs(delta) < 1e-16) next

        K[i, j] <- K[j, i] <- t_new

        W_new <- .w_update(W, i, j, delta)
        if (is.null(W_new)) {
          W <- solve(K)
        } else {
          W <- (W_new + t(W_new)) / 2
        }

        max_chg <- max(max_chg, abs(delta))
      }
    }

    if (verbose && iter %% 50 == 0)
      cat(sprintf("  [gauss_cd] iter %4d  max_chg=%.2e\n", iter, max_chg))
    if (max_chg < tol) { converged <- TRUE; break }
  }
  if (!converged)
    warning("gauss_golazo_cd: did not converge in ", maxiter, " sweeps")
  list(K = K, W = W, converged = converged, iters = iter)
}

hr_golazo_cd <- function(Gamma_hat, lambda, Lpen, Upen,
                          w_init = NULL,
                          tol = 1e-8, maxiter = 2000, verbose = FALSE,
                          active_rescan = 10L) {
  d <- nrow(Gamma_hat)

  if (is.null(w_init)) {
    Theta <- matrix(0, d, d)
    w0    <- 0.01
    for (i in seq_len(d)) for (j in seq_len(i - 1L)) {
      v <- rep(0, d); v[i] <- 1; v[j] <- -1
      Theta <- Theta + w0 * outer(v, v)
    }
  } else {
    Theta <- w_init
  }

  Tplus <- .sym_pinv(Theta)
  converged <- FALSE

  is_feasible <- (Lpen > -Inf) | (Upen < Inf)
  diag(is_feasible) <- FALSE
  feas_pairs <- which(is_feasible & lower.tri(is_feasible), arr.ind = TRUE)
  active <- rep(TRUE, nrow(feas_pairs))

  for (iter in seq_len(maxiter)) {
    max_chg <- 0

    do_full_sweep <- ((iter %% active_rescan) == 1L)
    sweep_idx <- if (do_full_sweep) seq_len(nrow(feas_pairs)) else which(active)

    for (k_ in sweep_idx) {
      i <- feas_pairs[k_, 1]
      j <- feas_pairs[k_, 2]
      tv    <- Tplus[, i] - Tplus[, j]
      alpha <- tv[i] - tv[j]
      if (alpha < 1e-14) { active[k_] <- FALSE; next }

      g <- alpha - Gamma_hat[i, j]

      c <- alpha^2

      Theta_old  <- Theta[i, j]
      Theta_new  <- soft_shrink_LU(Theta_old - g / c, lambda / c, Lpen[i, j], Upen[i, j])
      Theta_new  <- min(Theta_new, 0)
      w_new      <- -Theta_new
      delta_w    <- -(Theta_new - Theta_old)
      if (abs(delta_w) < 1e-16) {
        active[k_] <- FALSE
        next
      }

      Theta[i, j] <- Theta[j, i] <- Theta_new
      Theta[i, i] <- Theta[i, i] + delta_w
      Theta[j, j] <- Theta[j, j] + delta_w

      Tp_new <- .tplus_update_ij(Tplus, i, j, delta_w, tv = tv, alpha = alpha)
      if (is.null(Tp_new) || w_new == 0) {
        Tplus <- .sym_pinv(Theta)
      } else {
        Tplus <- (Tp_new + t(Tp_new)) / 2
      }

      active[k_] <- TRUE
      max_chg <- max(max_chg, abs(delta_w))
    }

    if (verbose && iter %% 50 == 0)
      cat(sprintf("  [hr_cd]   iter %4d  max_chg=%.2e  |active|=%d/%d\n",
                  iter, max_chg, sum(active), length(active)))
    if (max_chg < tol && do_full_sweep) { converged <- TRUE; break }
  }
  if (!converged)
    warning("hr_golazo_cd: did not converge in ", maxiter, " sweeps")
  list(Theta = Theta, Tplus = Tplus, converged = converged, iters = iter)
}
