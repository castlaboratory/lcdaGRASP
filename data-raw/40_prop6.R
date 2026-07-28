# data-raw/40_prop6.R
#
# Proposition 6 (empirically): how fast does the Reactive update concentrate
# p_k, and is the empirical winner k_hat actually k* at B = 150? We run
# LCDA-GR with B in {150, 300, 600, 1200} on Karate and a larger SBM, and
# record the Shannon entropy H(p_k) trajectory + end-state diagnostics.
#
# -> inst/extdata/prop6_summary.rds
# -> inst/extdata/prop6_trajectories.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

set.seed(2026)
g_small <- igraph::make_graph("Zachary")
g_med   <- igraph::sample_sbm(400,
            pref.matrix = { m <- matrix(0.02, 4, 4); diag(m) <- 0.15; m },
            block.sizes = rep(100, 4), directed = FALSE)
Bs <- c(150, 300, 600, 1200)

summary_rows <- list(); traj_rows <- list()
for (g_lbl in c("karate", "sbm400")) {
  g <- if (g_lbl == "karate") g_small else g_med
  for (B in Bs) {
    cli::cli_alert_info("prop6: {g_lbl} B={B}")
    r <- lcda_gr(g, variant = 1, B = B, m = 20, y = 60, seed = 1)
    H_pk <- apply(r$pk_history, 1, function(p) -sum(ifelse(p > 0, p * log(p), 0)))
    final_p <- r$pk_history[B, ]
    rho <- suppressWarnings(cor(rank(final_p), rank(r$mu_k), method = "spearman"))
    k_hat <- which.max(final_p)
    summary_rows[[paste0(g_lbl, B)]] <- tibble::tibble(
      graph = g_lbl, B = B, H_p_initial = H_pk[1], H_p_final = H_pk[B],
      H_max = log(r$params$m), rank_corr_p_mu = rho, k_hat = k_hat,
      mu_at_k_hat = r$mu_k[k_hat], max_mu_k = max(r$mu_k), best_Q = r$best$Q)
    traj_rows[[paste0(g_lbl, B)]] <- tibble::tibble(
      graph = g_lbl, B_max = B, iter = seq_len(B), H_pk = H_pk)
  }
}

meta <- make_meta(
  title = "Reactive update concentration (Proposition 6, empirical)",
  description = "Entropy of the selection distribution p_k over iterations for B in {150,300,600,1200}; initial value is log(m)=log(20). Lower = more concentrated.",
  network_source = paste(SRC_ZACHARY, "and", SRC_SBM, "-- SBM(400, 4 blocks, p_in=0.15, p_out=0.02)"),
  seed = 2026, m = 20, y = 60,
  limitations = c("Two graphs only (Karate + a 4-block SBM); concentration speed is graph-dependent."))
save_dataset("prop6_summary", dplyr::bind_rows(summary_rows), meta)
save_dataset("prop6_trajectories", dplyr::bind_rows(traj_rows), meta)
