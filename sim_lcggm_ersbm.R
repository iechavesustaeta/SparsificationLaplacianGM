# LCGGM simulation study over Erdos-Renyi and stochastic block model graphs.

suppressPackageStartupMessages({
  library(MASS); library(igraph); library(parallel)
})
source("experiment_helpers.R")

n_reps  <- 10
PROGRESS_LOG <- "../results/sim_lcggm_ersbm_progress.log"
if (!dir.exists("../results")) dir.create("../results")
unlink(PROGRESS_LOG)

log_progress <- function(msg) {
  line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line)
  cat(line, file = PROGRESS_LOG, append = TRUE)
}

gen_topology <- function(topo, d, seed = 42) {
  set.seed(seed)
  if (topo == "ER") {
    target_edges <- 2 * (d - 1)
    p <- target_edges / (d * (d - 1) / 2)
    g <- sample_gnp(d, p, directed = FALSE, loops = FALSE)
    tries <- 0
    while (!is_connected(g) && tries < 50) {
      g <- sample_gnp(d, p, directed = FALSE, loops = FALSE)
      tries <- tries + 1
    }
    label <- sprintf("ER(d=%d, p=%.4f)", d, p)
  } else if (topo == "SBM") {
    K_blocks    <- 5
    block_size  <- d %/% K_blocks
    block_sizes <- rep(block_size, K_blocks)
    rem <- d - sum(block_sizes); if (rem > 0) block_sizes[1] <- block_sizes[1] + rem
    pref <- matrix(0.01, K_blocks, K_blocks); diag(pref) <- 0.12
    g <- sample_sbm(d, pref.matrix = pref, block.sizes = block_sizes,
                    directed = FALSE, loops = FALSE)
    tries <- 0
    while (!is_connected(g) && tries < 50) {
      g <- sample_sbm(d, pref.matrix = pref, block.sizes = block_sizes,
                      directed = FALSE, loops = FALSE)
      tries <- tries + 1
    }
    label <- sprintf("SBM(d=%d, %d blocks, p_in=0.12, p_out=0.01)",
                     d, K_blocks)
  } else stop("unknown topology: ", topo)
  A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  W <- A * matrix(runif(d * d, 0.5, 1.5), d, d); W <- (W + t(W)) / 2
  K <- diag(rowSums(W)) - W
  list(K_true = K, A = A, g = g, label = label)
}

sample_lcggm <- function(n, K) {
  d <- nrow(K); Kp <- ginv(K)
  P <- diag(d) - matrix(1, d, d) / d
  Kpp <- P %*% Kp %*% P
  ev <- eigen(Kpp, symmetric = TRUE); pos <- ev$values > 1e-10
  U <- ev$vectors[, pos, drop = FALSE]
  Dm <- diag(sqrt(ev$values[pos]), nrow = sum(pos))
  matrix(rnorm(n * sum(pos)), n, sum(pos)) %*% t(U %*% Dm)
}

run_one_task <- function(task_idx, tasks, rep_data, E_true_mat,
                         setting_tag, progress_log) {
  source("experiment_helpers.R")
  rep <- tasks$rep[task_idx]
  mth <- tasks$method[task_idx]
  rd  <- rep_data[[rep]]

  cat(sprintf("[%s] [start] %s  task %d/%d  rep=%d  method=%s\n",
              format(Sys.time(), "%H:%M:%S"), setting_tag,
              task_idx, nrow(tasks), rep, mth),
      file = progress_log, append = TRUE)
  t_start <- proc.time()

  res <- if (mth == "cgl") {
    list(K = fit_cgl(rd$S_tr), param = 0)
  } else if (mth == "sgl_scad") {
    sel <- bic_select_alpha(rd$S_tr, rd$n_train, alpha_grid_default, fit_sgl_scad)
    list(K = sel$K, param = sel$param)
  } else if (mth == "spectral") {
    sel <- bic_select_spectral_lcggm(rd$S_tr, rd$n_train, eta_grid_lcggm_sim)
    list(K = sel$K, param = sel$param)
  } else stop("unknown method: ", mth)

  t_elapsed <- (proc.time() - t_start)["elapsed"]
  em <- edge_metrics(res$K, E_true_mat)
  ll <- loglik_lcggm(res$K, rd$S_te)
  ll_total <- (rd$n_test / 2) * ll

  cat(sprintf("[%s] [ end ] %s  task %d/%d  rep=%d  method=%-10s  edges=%4d  F1=%.3f  time=%6.1fs\n",
              format(Sys.time(), "%H:%M:%S"), setting_tag,
              task_idx, nrow(tasks), rep, mth,
              ifelse(is.na(em$n_edges), -1L, as.integer(em$n_edges)),
              ifelse(is.na(em$f1), 0, em$f1), t_elapsed),
      file = progress_log, append = TRUE)

  data.frame(
    rep = rep, method = method_labels_lcggm[[mth]],
    param = res$param, n_edges = em$n_edges,
    precision = em$precision, recall = em$recall, f1 = em$f1,
    test_loglik = ll_total, time_sec = t_elapsed,
    stringsAsFactors = FALSE
  )
}

run_setting <- function(d, topo, seed) {
  n_train <- 4 * d
  n_test  <- 4000
  td <- gen_topology(topo, d, seed = seed)
  K_true <- td$K_true
  E_true <- (abs(K_true) > 1e-10) & (row(K_true) != col(K_true))
  true_edges <- sum(E_true[upper.tri(E_true)])

  setting_tag <- sprintf("d=%d/%s", d, topo)
  log_progress(sprintf("========== %s ==========", setting_tag))
  log_progress(sprintf("%s  -- edges=%d, density=%.3f",
                       td$label, true_edges, true_edges / (d * (d - 1) / 2)))
  log_progress(sprintf("n_train=%d, n_test=%d, B=%d",
                       n_train, n_test, n_reps))

  log_progress("pre-generating rep data ...")
  t_gen <- proc.time()
  rep_data <- vector("list", n_reps)
  for (rep in seq_len(n_reps)) {
    set.seed(1000 + rep)
    X_all <- sample_lcggm(n_train + n_test, K_true)
    X_tr  <- X_all[seq_len(n_train), ]
    X_te  <- X_all[(n_train + 1):(n_train + n_test), ]
    rep_data[[rep]] <- list(
      X_tr = X_tr, X_te = X_te,
      S_tr = crossprod(X_tr) / n_train,
      S_te = crossprod(X_te) / n_test,
      n_train = n_train, n_test = n_test
    )
  }
  log_progress(sprintf("rep data ready (%.1fs)",
                       (proc.time() - t_gen)["elapsed"]))

  methods_in_order <- c("spectral", "sgl_scad", "cgl")
  tasks <- expand.grid(method = methods_in_order, rep = seq_len(n_reps),
                       stringsAsFactors = FALSE)
  tasks <- tasks[, c("rep", "method")]
  log_progress(sprintf("%d tasks (= %d reps x %d methods)",
                       nrow(tasks), n_reps, length(methods_in_order)))

  n_cores <- max(1L, parallel::detectCores())
  log_progress(sprintf("dispatching to %d workers (detectCores=%d)",
                       min(n_cores, nrow(tasks)), parallel::detectCores()))
  cl <- makeCluster(min(n_cores, nrow(tasks)))
  clusterExport(cl, varlist = c("tasks", "rep_data", "E_true",
                                "method_labels_lcggm", "method_order_lcggm",
                                "alpha_grid_default", "eta_grid_lcggm_sim",
                                "glasso_grid_default", "run_one_task",
                                "setting_tag", "PROGRESS_LOG"),
                envir = environment())
  clusterSetRNGStream(cl, iseed = seed + 1000)

  t0 <- proc.time()
  rows <- parLapplyLB(cl, seq_len(nrow(tasks)), function(i)
    run_one_task(i, tasks, rep_data, E_true, setting_tag, PROGRESS_LOG))
  stopCluster(cl)
  log_progress(sprintf("parallel phase done in %.1fs",
                       (proc.time() - t0)["elapsed"]))

  all_df <- do.call(rbind, rows)

  dead <- names(which(tapply(all_df$n_edges, all_df$method,
                             function(x) all(is.na(x)))))
  if (length(dead)) {
    stop("methods produced only NA in every replication: ",
         paste(dead, collapse = ", "),
         " -- check that their packages are installed (sparseGraph, ",
         "spectralGraphTopology) before trusting this run.")
  }
  n_na <- sum(is.na(all_df$n_edges))
  if (n_na > 0) log_progress(sprintf("WARNING: %d of %d (rep, method) fits failed",
                                     n_na, nrow(all_df)))

  if (!dir.exists("../results")) dir.create("../results")
  tag <- sprintf("sim_lcggm_ersbm_d%d_%s", d, tolower(topo))
  write.csv(all_df, sprintf("../results/%s_results.csv", tag), row.names = FALSE)

  agg <- aggregate(cbind(n_edges, precision, recall, f1, test_loglik, time_sec)
                   ~ method, data = all_df,
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                       se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))))
  agg <- do.call(data.frame, agg)
  names(agg) <- c("method",
                  "edges_mean","edges_se","precision_mean","precision_se",
                  "recall_mean","recall_se","f1_mean","f1_se",
                  "test_loglik_mean","test_loglik_se","time_sec_mean","time_sec_se")
  ord <- unname(method_labels_lcggm[method_order_lcggm])
  agg <- agg[match(ord, agg$method), ]
  agg <- agg[!is.na(agg$method), ]
  write.csv(agg, sprintf("../results/%s_summary.csv", tag), row.names = FALSE)

  cat(sprintf("\n  === Summary: %s ===\n", td$label))
  cat(sprintf("  %-26s %12s %12s %12s %12s %16s %10s\n",
              "method","edges","precision","recall","F1","test_loglik","time(s)"))
  cat("  ", strrep("-", 108), "\n", sep = "")
  for (i in seq_len(nrow(agg))) {
    cat(sprintf("  %-26s %5.1f (%4.1f) %5.3f (%5.3f) %5.3f (%5.3f) %5.3f (%5.3f) %9.0f (%5.0f) %7.1f\n",
                agg$method[i],
                agg$edges_mean[i], agg$edges_se[i],
                agg$precision_mean[i], agg$precision_se[i],
                agg$recall_mean[i], agg$recall_se[i],
                agg$f1_mean[i], agg$f1_se[i],
                agg$test_loglik_mean[i], agg$test_loglik_se[i],
                agg$time_sec_mean[i]))
  }
  invisible(list(label = td$label, true_edges = true_edges, agg = agg))
}

for (d in c(100, 200)) {
  for (topo in c("ER", "SBM")) {
    seed_ <- 100 * d + (if (topo == "ER") 1 else 2)
    run_setting(d, topo, seed = seed_)
  }
}

log_progress("ALL DONE")
