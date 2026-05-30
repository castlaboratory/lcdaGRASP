# data-raw/31_lfr_add_ecg.R
#
# Adds ECG (Ensemble Clustering for Graphs; Poulin & Theberge 2019) as a recent,
# runnable recovery baseline to the existing LFR sweep, WITHOUT re-running the
# expensive LCDA fits. Regenerates each LFR graph from the SAME seed used in
# 30_lfr_robustness.R, runs ECG, computes Q/ARI/NMI vs the planted truth, and
# appends method = "ECG" rows to the cached lfr_robustness.rds quality table.
#
# Literature (deep-research, verified): ECG reportedly beats Louvain/Leiden on
# LFR recovery and stays stable to higher mixing -- a fair recent comparator,
# and the consensus principle behind a possible LCDA improvement.

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages(library(dplyr))

# must match 30_lfr_robustness.R exactly
MUS  <- as.numeric(strsplit(Sys.getenv("LFR_MUS", "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"), ",")[[1]])
REPS <- as.integer(Sys.getenv("LFR_REPS", "10"))
N    <- as.integer(Sys.getenv("LFR_N", "500"))

prev <- readRDS(file.path(EXTDATA_DIR, "lfr_robustness.rds"))
stopifnot(!is.null(prev$results$quality))
grid <- tidyr::expand_grid(mu = MUS, rep = seq_len(REPS))
cli::cli_h1("Add ECG to LFR sweep: {nrow(grid)} cells")

ecg_rows <- list()
for (i in seq_len(nrow(grid))) {
  mu <- grid$mu[i]; rep <- grid$rep[i]
  gen <- lfr_generate(N, mu, seed = 1000 * rep + as.integer(mu * 100))  # SAME seed as 30_
  if (is.null(gen)) next
  g <- gen$graph; truth <- gen$truth
  t0 <- Sys.time()
  mem <- ecg_membership(g, ens = 16)
  if (is.null(mem)) { cli::cli_alert_warning("ECG failed (mu={mu}, rep={rep})"); next }
  ecg_rows[[length(ecg_rows) + 1]] <- tibble::tibble(
    mu = mu, rep = rep, method = "ECG",
    Q = igraph::modularity(g, mem), ARI = adjusted_rand(mem, truth),
    NMI = nmi(mem, truth), n_comms = length(unique(mem)),
    secs = as.numeric(Sys.time() - t0, units = "secs"))
  if (i %% 10 == 1) cli::cli_alert_info("cell {i}/{nrow(grid)} (mu={mu}) ECG NMI={round(nmi(mem,truth),3)}")
}

ecg_tbl <- dplyr::bind_rows(ecg_rows)
quality2 <- dplyr::bind_rows(prev$results$quality |> dplyr::filter(method != "ECG"), ecg_tbl)
prev$results$quality <- quality2
prev$meta$description <- paste(prev$meta$description, "[ECG baseline appended via 31_lfr_add_ecg.R.]")
prev$meta$generated_at <- as.character(Sys.time())
saveRDS(prev, file.path(EXTDATA_DIR, "lfr_robustness.rds"))
cli::cli_alert_success("ECG appended: {nrow(ecg_tbl)} rows. Mean NMI by mu:")
print(as.data.frame(ecg_tbl |> group_by(mu) |> summarise(NMI = mean(NMI), Q = mean(Q), .groups="drop")), digits = 3)
