# data-raw/11_blogs_timing.R
#
# E1 remediation: genuine execution-time measurement on Political Blogs for
# the four algorithms x 3 centralities x 3 similarities, with the paper's
# B = 50 (GRASP) vs B = 150 (GR). The paper's Table 9 (tab:sim_blogs) reports
# the LCDA-GR rows IDENTICAL to the LCDA-GRASP rows in both Q and time, which
# is impossible: with 3x the iterations the reactive variants must take longer.
# This script produces the companion-package C++ timings so the duplication
# can be replaced with reproducible numbers.
#
# Times are the companion package's optimized C++ kernels (so absolute values
# are smaller than an R-only implementation); the POINT is the B=50 vs B=150
# ratio and the non-degenerate per-configuration structure.
#
# -> inst/extdata/blogs_timing.rds  (per-run + summary tibbles)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

set.seed(20260528)
NREP <- as.integer(Sys.getenv("BLOGS_TIMING_NREP", "3"))  # timing mean is stable at low reps
cents <- c("betweenness", "closeness", "eigen")
sims  <- c("dice", "hpi", "jaccard")
algos <- c("LCDA-GRASP 1", "LCDA-GRASP 2", "LCDA-GR 1", "LCDA-GR 2")

g <- load_benchmark("polblogs")
cli::cli_h2("polblogs (n={igraph::vcount(g)}, m={igraph::ecount(g)}) - timing sweep, {NREP} reps")

run_one <- function(algo, cent, sim, B, seed) {
  t0 <- Sys.time()
  res <- switch(algo,
    "LCDA-GRASP 1" = lcda_grasp(g, variant = 1, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GRASP 2" = lcda_grasp(g, variant = 2, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GR 1"    = lcda_gr(g, variant = 1, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GR 2"    = lcda_gr(g, variant = 2, B = B, centrality = cent, similarity = sim, seed = seed)
  )
  tibble::tibble(algo = algo, centrality = cent, similarity = sim, B = B,
                 seed = seed, Q = res$best$Q,
                 g_comms = length(res$best$leaders),
                 secs = as.numeric(Sys.time() - t0, units = "secs"))
}

grid <- expand.grid(algo = algos, centrality = cents, similarity = sims,
                    stringsAsFactors = FALSE)
per_run <- list()
for (i in seq_len(nrow(grid))) {
  algo <- grid$algo[i]; cent <- grid$centrality[i]; sim <- grid$similarity[i]
  # NOTE: match "GRASP" explicitly -- grepl("GR", "LCDA-GRASP 1") is TRUE, which
  # would wrongly give every algorithm B=150. GRASP variants use B=50, GR use B=150.
  B <- if (grepl("GRASP", algo)) 50L else 150L
  df <- purrr::map_dfr(seq_len(NREP), function(r)
    tryCatch(run_one(algo, cent, sim, B, seed = r), error = function(e) NULL))
  per_run[[paste(algo, cent, sim)]] <- df
  cli::cli_alert("{algo} / {cent} / {sim} (B={B}): time mean={round(mean(df$secs),2)}s  Q mean={round(mean(df$Q),4)}")
}

results <- dplyr::bind_rows(per_run)
summary_tbl <- results |>
  dplyr::group_by(algo, centrality, similarity, B) |>
  dplyr::summarise(Q_mean = mean(Q), Q_sd = sd(Q),
                   secs_mean = mean(secs), secs_sd = sd(secs),
                   g_comms = mean(g_comms), .groups = "drop")

meta <- make_meta(
  title = "Political Blogs timing sweep (E1 remediation): 4 algorithms x 3 centralities x 3 similarities",
  description = paste("Companion-package C++ execution times and mean Q on Political Blogs,",
                      "B=50 (GRASP) vs B=150 (GR). Produced to replace the duplicated",
                      "LCDA-GR rows of the paper's Table 9 (tab:sim_blogs)."),
  network_source = src_newman("polblogs"),
  seed = 20260528, nrep = NREP,
  limitations = c(
    "Absolute times are the package's optimized C++ kernels; an R-only implementation (as in the paper) is far slower. The reproducible content is the B=50 vs B=150 ratio and the per-configuration structure.",
    "Closeness/betweenness on the disconnected PolBlogs (giant = 1222 of 1490) can degenerate: the non-finite-centrality fallback may collapse the partition (Q near 0). This is a separate package finding, flagged for the weighted/large-scale work.",
    "PolBlogs collapsed directed->undirected keeping all 1490 nodes."))

save_dataset("blogs_timing", list(per_run = results, summary = summary_tbl), meta)
cli::cli_alert_success("E1 timing sweep done.")
