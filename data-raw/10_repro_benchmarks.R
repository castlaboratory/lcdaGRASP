# data-raw/10_repro_benchmarks.R
#
# Reproduce the paper's headline modularity numbers on the five benchmark
# networks, keeping the FULL distribution over replications (the paper
# reports only means/maxima). Louvain and Leiden serve as baselines.
#
# -> inst/extdata/repro_benchmarks.rds  (per-run tibble)
# -> inst/extdata/repro_summary.rds     (aggregated + baselines + literature)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

set.seed(20260527)
NREP <- as.integer(Sys.getenv("REPRO_NREP", "30"))
algos <- c("LCDA-GRASP 1", "LCDA-GRASP 2", "LCDA-GR 1", "LCDA-GR 2")

run_one <- function(g, algo, B, seed) {
  t0 <- Sys.time()
  res <- switch(algo,
    "LCDA-GRASP 1" = lcda_grasp(g, variant = 1, B = B, seed = seed),
    "LCDA-GRASP 2" = lcda_grasp(g, variant = 2, B = B, seed = seed),
    "LCDA-GR 1"    = lcda_gr(g, variant = 1, B = B, seed = seed),
    "LCDA-GR 2"    = lcda_gr(g, variant = 2, B = B, seed = seed)
  )
  tibble::tibble(algo = algo, seed = seed, Q = res$best$Q, H = res$best$H,
                 g_comms = length(res$best$leaders),
                 secs = as.numeric(Sys.time() - t0, units = "secs"))
}

per_run <- list(); baselines <- list()
for (net in BENCHMARKS$network) {
  g <- load_benchmark(net)
  cli::cli_h2("{net} (n={igraph::vcount(g)}, m={igraph::ecount(g)}) - {NREP} reps")

  lq <- vapply(seq_len(30), function(s) { set.seed(1000 + s); igraph::modularity(igraph::cluster_louvain(g)) }, numeric(1))
  leid <- tryCatch(vapply(seq_len(30), function(s) {
    set.seed(2000 + s)
    cl <- igraph::cluster_leiden(g, objective_function = "modularity")
    igraph::modularity(g, igraph::membership(cl))
  }, numeric(1)), error = function(e) rep(NA_real_, 30))
  baselines[[net]] <- tibble::tibble(network = net,
    Q_louvain_mean = mean(lq), Q_louvain_max = max(lq), Q_louvain_sd = sd(lq),
    Q_leiden_mean = mean(leid), Q_leiden_max = max(leid))
  cli::cli_alert_info("Louvain max={round(max(lq),4)}  Leiden max={round(max(leid),4)}")

  for (algo in algos) {
    # match "GRASP" explicitly: grepl("GR", "LCDA-GRASP 1") is TRUE and would
    # wrongly run the fixed-parameter GRASP variants at B=150 too. (Fixed 2026-05-28.)
    B <- if (grepl("GRASP", algo)) 50 else 150
    df <- purrr::map_dfr(seq_len(NREP), function(r)
      tryCatch(run_one(g, algo, B, seed = r), error = function(e) NULL))
    df$network <- net
    per_run[[paste(net, algo)]] <- df
    cli::cli_alert("{algo}: Q mean={round(mean(df$Q),4)} sd={round(sd(df$Q),4)} max={round(max(df$Q),4)}")
  }
}

results <- dplyr::bind_rows(per_run)
baseline_tbl <- dplyr::bind_rows(baselines)

summary_tbl <- results |>
  dplyr::group_by(network, algo) |>
  dplyr::summarise(Q_mean = mean(Q), Q_sd = sd(Q), Q_min = min(Q), Q_max = max(Q),
                   g_comms = mean(g_comms), secs = mean(secs), .groups = "drop") |>
  dplyr::left_join(BENCHMARKS[, c("network", "Q_ref", "kind")], by = "network") |>
  dplyr::left_join(baseline_tbl, by = "network") |>
  dplyr::mutate(pct_vs_louvain = 100 * (Q_max - Q_louvain_max) / Q_louvain_max,
                gap_vs_ref = Q_ref - Q_max)

meta <- make_meta(
  title = "Benchmark reproduction (5 networks x 4 LCDA variants + Louvain/Leiden)",
  description = "Best/mean/sd modularity over reps, vs igraph Louvain/Leiden and literature best-known Q. Confirms E3 (Louvain on Karate reaches 0.4198, not 0.3715).",
  seed = 20260527, nrep = NREP,
  limitations = c(
    "alpha_c/alpha_s use package defaults (0.1, 0.3), not the paper's exact RCL parameters.",
    "PolBlogs collapsed directed->undirected keeping all 1490 nodes (giant component = 1222).",
    "Reference Q are literature best-known; only Karate is a certified optimum."))

save_dataset("repro_benchmarks", results, meta)
save_dataset("repro_summary", summary_tbl, meta)
