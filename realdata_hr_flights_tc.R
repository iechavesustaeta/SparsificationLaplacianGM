# Texas Cluster flight-delay experiment.

suppressPackageStartupMessages({
  library(MASS); library(igraph); library(graphicalExtremes); library(digest)
})
source("experiment_helpers.R")

set.seed(20260816)

iatas    <- getFlightDelayData(what = "IATAs",
                               airportFilter = "tcCluster")
X_tr     <- getFlightDelayData(what = "delays",
                               airportFilter = "tcCluster",
                               dateFilter    = "tcTrain",
                               delayFilter   = "totals")
X_te     <- getFlightDelayData(what = "delays",
                               airportFilter = "tcCluster",
                               dateFilter    = "tcTest",
                               delayFilter   = "totals")
if (length(dim(X_tr)) == 3) X_tr <- X_tr[ , , 1]
if (length(dim(X_te)) == 3) X_te <- X_te[ , , 1]

d <- ncol(X_tr); n_train <- nrow(X_tr); n_test <- nrow(X_te)
cat(sprintf("Texas Cluster flights: d=%d airports (%s)\n",
            d, paste(iatas, collapse = ", ")))
cat(sprintf("Pre-split:  n_train=%d, n_test=%d\n", n_train, n_test))

cat(sprintf("NA fraction train: %.3f, test: %.3f\n",
            mean(is.na(X_tr)), mean(is.na(X_te))))
keep_tr <- rowSums(!is.na(X_tr)) > 0
keep_te <- rowSums(!is.na(X_te)) > 0
X_tr <- X_tr[keep_tr, , drop = FALSE]
X_te <- X_te[keep_te, , drop = FALSE]
n_train <- nrow(X_tr); n_test <- nrow(X_te)
cat(sprintf("After dropping fully-NA rows: n_train=%d, n_test=%d\n",
            n_train, n_test))

p_thresh <- 0.85
Gamma_tr <- emp_vario(data = X_tr, p = p_thresh)
Gamma_te <- emp_vario(data = X_te, p = p_thresh)

data_te_mpareto <- data2mpareto(data = X_te, p = p_thresh)
n_pareto_te <- nrow(data_te_mpareto)
data_tr_mpareto <- data2mpareto(data = X_tr, p = p_thresh)
n_pareto_tr <- nrow(data_tr_mpareto)
cat(sprintf("Multivariate-Pareto exceedance counts: train=%d, test=%d  (p=%.2f)\n",
            n_pareto_tr, n_pareto_te, p_thresh))

.ll_cache <- new.env(parent = emptyenv())
loglik_hr_cached <- function(Theta, data_mpareto, tag) {
  if (is.null(Theta)) return(NA_real_)
  key <- paste0(tag, "-", digest::digest(round(Theta, 10)))
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

bic_hr_pareto <- function(Theta, data_mpareto) {
  ll <- loglik_hr_cached(Theta, data_mpareto, "tr")
  if (is.null(Theta) || is.na(ll)) return(Inf)
  -2 * ll + count_edges(Theta) * log(nrow(data_mpareto))
}

rho_grid    <- sort(unique(c(seq(0.005, 0.05, length.out = 6),
                             seq(0.05,  0.50, length.out = 12))))
rho_gl_grid <- sort(unique(c(seq(0.005, 0.05, length.out = 6),
                             seq(0.05,  0.50, length.out = 12))))
eta_grid    <- sort(unique(c(c(1.1, 1.2, 1.3),
                             eta_grid_default,
                             c(10.0, 15.0, 20.0, 25.0,
                               30.0, 40.0, 50.0, 75.0, 100.0),
                             c(150.0, 200.0, 300.0, 500.0))))
cat(sprintf("Grids: |rho_ns|=%d, |rho_glasso|=%d, |eta|=%d\n",
            length(rho_grid), length(rho_gl_grid), length(eta_grid)))

aic_hr_pareto <- function(Theta, data_mpareto) {
  ll <- loglik_hr_cached(Theta, data_mpareto, "tr")
  if (is.null(Theta) || is.na(ll)) return(Inf)
  -2 * ll + 2 * count_edges(Theta)
}

sweep_grid <- function(thetas, params) {
  do.call(rbind, lapply(seq_along(params), function(j) {
    Th <- thetas[[j]]
    if (is.null(Th)) {
      return(data.frame(param = params[j], edges = NA_integer_,
                        loglik = NA_real_, BIC = NA_real_, AIC = NA_real_))
    }
    ll <- loglik_hr_cached(Th, data_tr_mpareto, "tr")
    k  <- count_edges(Th)
    if (is.na(ll)) {
      return(data.frame(param = params[j], edges = k,
                        loglik = NA_real_, BIC = NA_real_, AIC = NA_real_))
    }
    data.frame(param = params[j], edges = k, loglik = ll,
               BIC = -2 * ll + k * log(n_pareto_tr),
               AIC = -2 * ll + 2 * k)
  }))
}

pick_by <- function(sweep, thetas, crit) {
  ok <- which(is.finite(sweep[[crit]]))
  if (!length(ok)) return(list(K = NULL, param = NA_real_))
  j <- ok[which.min(sweep[[crit]][ok])]
  list(K = thetas[[j]], param = sweep$param[j])
}

print_sweep <- function(tag, sweep, sel_bic, sel_aic, fmt = "%.4f") {
  cat(sprintf("\n%s sweep (train BIC / AIC):\n", tag))
  for (j in seq_len(nrow(sweep))) {
    mark <- paste0(
      if (isTRUE(all.equal(sweep$param[j], sel_bic))) "  <-- BIC" else "",
      if (isTRUE(all.equal(sweep$param[j], sel_aic))) "  <-- AIC" else "")
    if (is.na(sweep$BIC[j])) {
      cat(sprintf(paste0("  ", fmt, "  ---\n"), sweep$param[j]))
    } else {
      cat(sprintf(paste0("  ", fmt, "  edges=%3d  BIC=%10.2f  AIC=%10.2f%s\n"),
                  sweep$param[j], sweep$edges[j], sweep$BIC[j], sweep$AIC[j], mark))
    }
  }
}

fits <- list()

cat("\n[1/4] EMTP2 (GOLAZO, L=0/U=Inf) ... ")
t1 <- system.time({ K_emtp2 <- fit_emtp2(Gamma_tr) })["elapsed"]
fits[["EMTP2|(none)"]] <- list(K = K_emtp2, param = NA_real_, time = t1,
                               label = "EMTP2", criterion = "(none)")
cat(sprintf("done (%.1fs)\n", t1))

cat("[2/4] eglearn-NS (sweep over rho) ... ")
t2 <- system.time({
  eg_fit <- tryCatch(graphicalExtremes::eglearn(
    data = X_tr, p = p_thresh, rholist = rho_grid,
    reg_method = "ns", complete_Gamma = TRUE),
    error = function(e) NULL)
  eg_thetas <- lapply(seq_along(rho_grid), function(j) {
    if (is.null(eg_fit)) return(NULL)
    G_j <- eg_fit$Gamma[[j]]
    if (is.null(G_j) || anyNA(G_j)) return(NULL)
    tryCatch(graphicalExtremes::Gamma2Theta(G_j), error = function(e) NULL)
  })
  eg_sweep <- sweep_grid(eg_thetas, rho_grid)
})["elapsed"]
eg_bic <- pick_by(eg_sweep, eg_thetas, "BIC")
eg_aic <- pick_by(eg_sweep, eg_thetas, "AIC")
fits[["eglearn-NS|BIC"]] <- list(K = eg_bic$K, param = eg_bic$param, time = t2,
                                 label = "eglearn-NS (Engelke et al)", criterion = "BIC")
fits[["eglearn-NS|AIC"]] <- list(K = eg_aic$K, param = eg_aic$param, time = t2,
                                 label = "eglearn-NS (Engelke et al)", criterion = "AIC")
cat(sprintf("done (%.1fs, BIC rho=%.4f, AIC rho=%.4f)\n", t2, eg_bic$param, eg_aic$param))

cat("[3/4] eglearn-glasso (sweep over rho) ... ")
t3 <- system.time({
  egl_fit <- tryCatch(graphicalExtremes::eglearn(
    data = X_tr, p = p_thresh, rholist = rho_gl_grid,
    reg_method = "glasso", complete_Gamma = TRUE),
    error = function(e) NULL)
  egl_thetas <- lapply(seq_along(rho_gl_grid), function(j) {
    if (is.null(egl_fit)) return(NULL)
    G_j <- egl_fit$Gamma[[j]]
    if (is.null(G_j) || anyNA(G_j)) return(NULL)
    tryCatch(graphicalExtremes::Gamma2Theta(G_j), error = function(e) NULL)
  })
  egl_sweep <- sweep_grid(egl_thetas, rho_gl_grid)
})["elapsed"]
egl_bic <- pick_by(egl_sweep, egl_thetas, "BIC")
egl_aic <- pick_by(egl_sweep, egl_thetas, "AIC")
fits[["eglearn-glasso|BIC"]] <- list(K = egl_bic$K, param = egl_bic$param, time = t3,
                                     label = "eglearn-glasso (Engelke et al)", criterion = "BIC")
fits[["eglearn-glasso|AIC"]] <- list(K = egl_aic$K, param = egl_aic$param, time = t3,
                                     label = "eglearn-glasso (Engelke et al)", criterion = "AIC")
cat(sprintf("done (%.1fs, BIC rho=%.4f, AIC rho=%.4f)\n", t3, egl_bic$param, egl_aic$param))

cat("[4/4] Spectral-HR (sweep over eta) ... ")
t4 <- system.time({
  mle <- hr_mle_dense(Gamma_tr, tol = 1e-8, maxiter = 3000)
  sp_thetas <- lapply(eta_grid, function(eta) {
    res_eta <- tryCatch(spectral_hr_from_mle(mle, eta, tol = 1e-8, maxiter = 3000),
                        error = function(e) NULL)
    if (is.null(res_eta)) NULL else res_eta$Theta_refit
  })
  sp_sweep <- sweep_grid(sp_thetas, eta_grid)
})["elapsed"]
sp_bic <- pick_by(sp_sweep, sp_thetas, "BIC")
sp_aic <- pick_by(sp_sweep, sp_thetas, "AIC")
fits[["Spectral-HR|BIC"]] <- list(K = sp_bic$K, param = sp_bic$param, time = t4,
                                  label = "Spectral-HR (ours)", criterion = "BIC")
fits[["Spectral-HR|AIC"]] <- list(K = sp_aic$K, param = sp_aic$param, time = t4,
                                  label = "Spectral-HR (ours)", criterion = "AIC")
cat(sprintf("done (%.1fs, BIC eta=%.3g, AIC eta=%.3g)\n", t4, sp_bic$param, sp_aic$param))

print_sweep("eglearn-NS (rho)",     eg_sweep,  eg_bic$param,  eg_aic$param)
print_sweep("eglearn-glasso (rho)", egl_sweep, egl_bic$param, egl_aic$param)
print_sweep("Spectral-HR (eta)",    sp_sweep,  sp_bic$param,  sp_aic$param, fmt = "%6.2f")
cat("\n--- end diagnostics ---\n\n")

if (!dir.exists("../results")) dir.create("../results")
write.csv(setNames(eg_sweep,  c("rho", "edges", "loglik", "BIC", "AIC")),
          "../results/realdata_hr_flights_tc_eglearn_ns_path.csv", row.names = FALSE)
write.csv(setNames(egl_sweep, c("rho", "edges", "loglik", "BIC", "AIC")),
          "../results/realdata_hr_flights_tc_eglearn_glasso_path.csv", row.names = FALSE)
write.csv(setNames(sp_sweep,  c("eta", "edges", "loglik", "BIC", "AIC")),
          "../results/realdata_hr_flights_tc_spectral_path.csv", row.names = FALSE)

out <- data.frame()
for (nm in names(fits)) {
  f <- fits[[nm]]
  if (is.null(f) || is.null(f$K)) next
  test_ll   <- loglik_hr_cached(f$K, data_te_mpareto, "te")
  train_bic <- bic_hr_pareto(f$K, data_tr_mpareto)
  train_aic <- aic_hr_pareto(f$K, data_tr_mpareto)
  out <- rbind(out, data.frame(
    method = f$label, criterion = f$criterion, param = f$param,
    n_edges = count_edges(f$K), train_bic = train_bic, train_aic = train_aic,
    test_loglik = test_ll, time_sec = f$time, stringsAsFactors = FALSE
  ))
}

write.csv(out, "../results/realdata_hr_flights_tc.csv", row.names = FALSE)

cat(sprintf("\n=== Texas Cluster HR (d=%d, n_train=%d, n_test=%d, p=%.2f) ===\n",
            d, n_train, n_test, p_thresh))
cat(sprintf("Pareto exceedance counts: train=%d, test=%d.\n",
            n_pareto_tr, n_pareto_te))
cat(sprintf("\n  %-32s %6s %10s %6s %13s %13s %14s %8s\n",
            "method", "sel.", "param", "edges", "train_BIC", "train_AIC",
            "test_loglik", "time(s)"))
cat("  ", strrep("-", 112), "\n", sep = "")
for (i in seq_len(nrow(out))) {
  cat(sprintf("  %-32s %6s %10.3g %6d %13.2f %13.2f %14.2f %7.1f\n",
              out$method[i], out$criterion[i], out$param[i], out$n_edges[i],
              out$train_bic[i], out$train_aic[i], out$test_loglik[i], out$time_sec[i]))
}

gap <- with(out, (train_bic - train_aic) - n_edges * (log(n_pareto_tr) - 2))
cat(sprintf("\nmax |BIC - AIC - k(log n - 2)| = %.3e\n", max(abs(gap))))

cat("\nResults saved to ../results/realdata_hr_flights_tc.csv\n")
cat("Path sweeps saved to ../results/realdata_hr_flights_tc_*_path.csv\n")
