# Upper Danube river-network discharge experiment.

suppressPackageStartupMessages({
  library(MASS); library(igraph); library(graphicalExtremes); library(digest)
})
source("experiment_helpers.R")

set.seed(20260816)

data(danube, package = "graphicalExtremes")
X <- danube$data_clustered
d <- ncol(X); n <- nrow(X)
cat(sprintf("Danube data: n=%d, d=%d\n", n, d))

p_thresh <- 0.9
Gamma_full <- emp_vario(data = X, p = p_thresh)
data_mpareto <- data2mpareto(data = X, p = p_thresh)
n_pareto <- nrow(data_mpareto)
cat(sprintf("Threshold p=%.2f, full multivariate-Pareto sample size: n_Pareto=%d\n",
            p_thresh, n_pareto))

.ll_cache <- new.env(parent = emptyenv())
loglik_hr_full <- function(Theta) {
  if (is.null(Theta)) return(NA_real_)
  key <- digest::digest(round(Theta, 10))
  hit <- .ll_cache[[key]]
  if (!is.null(hit)) return(hit)
  Gamma_hat <- tryCatch(Theta2Gamma(Theta), error = function(e) NULL)
  if (is.null(Gamma_hat)) return(NA_real_)
  ll <- tryCatch(as.numeric(loglik_HR(data = data_mpareto,
                                      Gamma = Gamma_hat)[["loglik"]]),
                 error = function(e) NA_real_)
  ll <- if (is.na(ll) || !is.finite(ll)) NA_real_ else ll
  .ll_cache[[key]] <- ll
  ll
}
bic_hr_pareto <- function(Theta) {
  ll <- loglik_hr_full(Theta)
  if (is.na(ll)) return(Inf)
  -2 * ll + count_edges(Theta) * log(n_pareto)
}
aic_hr_pareto <- function(Theta) {
  ll <- loglik_hr_full(Theta)
  if (is.na(ll)) return(Inf)
  -2 * ll + 2 * count_edges(Theta)
}

rho_grid    <- sort(unique(c(seq(0.005, 0.05, length.out = 6),
                             seq(0.05,  0.50, length.out = 12))))
rho_gl_grid <- seq(0.005, 0.750, length.out = 20)
eta_grid    <- sort(unique(c(c(1.1, 1.2, 1.3),
                             eta_grid_default,
                             c(10.0, 15.0, 20.0, 25.0,
                               30.0, 40.0, 50.0, 75.0, 100.0))))
cat(sprintf("Grids: |rho_ns|=%d, |rho_glasso|=%d, |eta|=%d\n",
            length(rho_grid), length(rho_gl_grid), length(eta_grid)))

select_best <- function(theta_list, params) {
  bics <- vapply(theta_list, bic_hr_pareto, numeric(1))
  aics <- vapply(theta_list, aic_hr_pareto, numeric(1))
  best_b <- which.min(bics)
  best_a <- which.min(aics)
  list(
    bic = list(Theta = theta_list[[best_b]], param = params[best_b],
               bic = bics[best_b], aic = aic_hr_pareto(theta_list[[best_b]])),
    aic = list(Theta = theta_list[[best_a]], param = params[best_a],
               bic = bic_hr_pareto(theta_list[[best_a]]), aic = aics[best_a])
  )
}

cat("\nFitting...\n")

t0 <- proc.time()
K_emtp2 <- fit_emtp2(Gamma_full)
time_emtp2 <- (proc.time() - t0)[["elapsed"]]
emtp2_fit <- list(Theta = K_emtp2, param = NA_real_,
                  bic = bic_hr_pareto(K_emtp2),
                  aic = aic_hr_pareto(K_emtp2))

cat("  eglearn-NS sweep ...\n")
t0 <- proc.time()
eg_fit <- graphicalExtremes::eglearn(
  data = X, p = p_thresh, rholist = rho_grid,
  reg_method = "ns", complete_Gamma = TRUE)
eg_thetas <- lapply(seq_along(rho_grid), function(j) {
  G_j <- eg_fit$Gamma[[j]]
  if (is.null(G_j) || anyNA(G_j)) return(NULL)
  tryCatch(graphicalExtremes::Gamma2Theta(G_j), error = function(e) NULL)
})
time_eg <- (proc.time() - t0)[["elapsed"]]
eg_valid <- !sapply(eg_thetas, is.null)
eg_sel <- select_best(eg_thetas[eg_valid], rho_grid[eg_valid])

cat("  eglearn-glasso sweep ...\n")
t0 <- proc.time()
egl_fit <- graphicalExtremes::eglearn(
  data = X, p = p_thresh, rholist = rho_gl_grid,
  reg_method = "glasso", complete_Gamma = TRUE)
egl_thetas <- lapply(seq_along(rho_gl_grid), function(j) {
  G_j <- egl_fit$Gamma[[j]]
  if (is.null(G_j) || anyNA(G_j)) return(NULL)
  tryCatch(graphicalExtremes::Gamma2Theta(G_j), error = function(e) NULL)
})
time_egl <- (proc.time() - t0)[["elapsed"]]
egl_valid <- !sapply(egl_thetas, is.null)
egl_sel <- select_best(egl_thetas[egl_valid], rho_gl_grid[egl_valid])

cat("  Spectral-HR sweep ...\n")
t0 <- proc.time()
mle <- hr_mle_dense(Gamma_full, tol = 1e-8, maxiter = 3000)
sp_thetas <- lapply(eta_grid, function(eta) {
  r <- tryCatch(spectral_hr_from_mle(mle, eta, tol = 1e-8, maxiter = 3000),
                error = function(e) NULL)
  if (is.null(r) || is.null(r$Theta_refit)) NULL else r$Theta_refit
})
time_sp <- (proc.time() - t0)[["elapsed"]]
sp_valid <- !sapply(sp_thetas, is.null)
sp_sel <- select_best(sp_thetas[sp_valid], eta_grid[sp_valid])

cat("\n--- HR-Pareto BIC/AIC sweep diagnostics ---\n")
print_sweep <- function(label, thetas, params, valid, sel_bic_p, sel_aic_p) {
  cat(sprintf("\n%s:\n", label))
  for (j in seq_along(thetas)) {
    if (!valid[j]) {
      cat(sprintf("  param=%8.4f  --- (no valid Gamma)\n", params[j])); next
    }
    Th <- thetas[[j]]
    mark <- ""
    if (isTRUE(all.equal(params[j], sel_bic_p))) mark <- paste0(mark, " <-- BIC")
    if (isTRUE(all.equal(params[j], sel_aic_p))) mark <- paste0(mark, " <-- AIC")
    cat(sprintf("  param=%8.4f  edges=%4d  BIC=%9.2f  AIC=%9.2f%s\n",
                params[j], count_edges(Th),
                bic_hr_pareto(Th), aic_hr_pareto(Th), mark))
  }
}
print_sweep("eglearn-NS", eg_thetas, rho_grid, eg_valid,
            eg_sel$bic$param, eg_sel$aic$param)
print_sweep("eglearn-glasso", egl_thetas, rho_gl_grid, egl_valid,
            egl_sel$bic$param, egl_sel$aic$param)
print_sweep("Spectral-HR", sp_thetas, eta_grid, sp_valid,
            sp_sel$bic$param, sp_sel$aic$param)
cat("--- end diagnostics ---\n\n")

make_row <- function(method, criterion, fit, time_sec) {
  data.frame(method = method, criterion = criterion,
             param = fit$param, edges = count_edges(fit$Theta),
             BIC = fit$bic, AIC = fit$aic, time_sec = round(time_sec, 2),
             stringsAsFactors = FALSE)
}
out <- rbind(
  make_row("EMTP2",          "(none)", emtp2_fit, time_emtp2),
  make_row("eglearn-NS",     "BIC",    eg_sel$bic, time_eg),
  make_row("eglearn-NS",     "AIC",    eg_sel$aic, time_eg),
  make_row("eglearn-glasso", "BIC",    egl_sel$bic, time_egl),
  make_row("eglearn-glasso", "AIC",    egl_sel$aic, time_egl),
  make_row("Spectral-HR",    "BIC",    sp_sel$bic, time_sp),
  make_row("Spectral-HR",    "AIC",    sp_sel$aic, time_sp)
)

if (!dir.exists("../results")) dir.create("../results")
write.csv(out, "../results/realdata_hr_danube_nosplit.csv", row.names = FALSE)

sp_full <- data.frame()
for (j in seq_along(eta_grid)) {
  if (!sp_valid[j]) {
    sp_full <- rbind(sp_full, data.frame(
      eta = eta_grid[j], edges = NA_integer_,
      BIC = NA_real_, AIC = NA_real_,
      stringsAsFactors = FALSE))
    next
  }
  Th <- sp_thetas[[j]]
  sp_full <- rbind(sp_full, data.frame(
    eta = eta_grid[j], edges = count_edges(Th),
    BIC = bic_hr_pareto(Th), AIC = aic_hr_pareto(Th),
    stringsAsFactors = FALSE))
}
write.csv(sp_full, "../results/realdata_hr_danube_nosplit_spectral_path.csv",
          row.names = FALSE)

eg_full <- data.frame()
for (j in seq_along(rho_grid)) {
  if (!eg_valid[j]) {
    eg_full <- rbind(eg_full, data.frame(
      rho = rho_grid[j], edges = NA_integer_,
      BIC = NA_real_, AIC = NA_real_,
      stringsAsFactors = FALSE))
    next
  }
  Th <- eg_thetas[[j]]
  eg_full <- rbind(eg_full, data.frame(
    rho = rho_grid[j], edges = count_edges(Th),
    BIC = bic_hr_pareto(Th), AIC = aic_hr_pareto(Th),
    stringsAsFactors = FALSE))
}
write.csv(eg_full, "../results/realdata_hr_danube_nosplit_eglearn_ns_path.csv",
          row.names = FALSE)

egl_full <- data.frame()
for (j in seq_along(rho_gl_grid)) {
  if (!egl_valid[j]) {
    egl_full <- rbind(egl_full, data.frame(
      rho = rho_gl_grid[j], edges = NA_integer_,
      BIC = NA_real_, AIC = NA_real_,
      stringsAsFactors = FALSE))
    next
  }
  Th <- egl_thetas[[j]]
  egl_full <- rbind(egl_full, data.frame(
    rho = rho_gl_grid[j], edges = count_edges(Th),
    BIC = bic_hr_pareto(Th), AIC = aic_hr_pareto(Th),
    stringsAsFactors = FALSE))
}
write.csv(egl_full, "../results/realdata_hr_danube_nosplit_eglearn_glasso_path.csv",
          row.names = FALSE)

cat(sprintf("=== Danube HR, no split (d=%d, n_obs=%d, n_Pareto=%d, p=%.2f) ===\n",
            d, n, n_pareto, p_thresh))
cat(sprintf("  %-16s %-8s %8s %6s %10s %10s\n",
            "method", "crit", "param", "edges", "BIC", "AIC"))
cat("  ", strrep("-", 64), "\n", sep = "")
for (i in seq_len(nrow(out))) {
  cat(sprintf("  %-16s %-8s %8.4g %6d %10.2f %10.2f\n",
              out$method[i], out$criterion[i], out$param[i], out$edges[i],
              out$BIC[i], out$AIC[i]))
}
cat("\nResults saved to ../results/realdata_hr_danube_nosplit.csv\n")

cat(sprintf("\n=== Spectral-HR full eta-path (all grid points) ===\n"))
cat(sprintf("  %5s %6s %10s %10s%s\n", "eta", "edges", "BIC", "AIC", ""))
cat("  ", strrep("-", 46), "\n", sep = "")
for (i in seq_len(nrow(sp_full))) {
  r <- sp_full[i, ]
  if (is.na(r$edges)) { cat(sprintf("  %5.2f  ---\n", r$eta)); next }
  mark <- ""
  if (isTRUE(all.equal(r$eta, sp_sel$bic$param))) mark <- paste0(mark, " <-- BIC")
  if (isTRUE(all.equal(r$eta, sp_sel$aic$param))) mark <- paste0(mark, " <-- AIC")
  cat(sprintf("  %5.2f %6d %10.2f %10.2f%s\n",
              r$eta, r$edges, r$BIC, r$AIC, mark))
}
cat("\nSpectral-HR full path saved to ../results/realdata_hr_danube_nosplit_spectral_path.csv\n")
