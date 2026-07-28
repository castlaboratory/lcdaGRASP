# data-raw/50_pool_sensitivity.R
#
# The paper fixes m=20, y=3m=60, B=150 without testing alternatives. Sweep
#   m in {5,10,20,40,80}, y/m in {1,2,3,5,10}, B in {60,120,240}
# and record best Q, efficiency, and the lex-decisive fraction (how often the
# H tie-break actually bound).
#
# -> inst/extdata/pool_sensitivity.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages(library(purrr))

set.seed(2026)
g <- igraph::sample_sbm(300,
       pref.matrix = { m <- matrix(0.02, 5, 5); diag(m) <- 0.12; m },
       block.sizes = rep(60, 5))

REPS <- 5L   # publication-grade for this grid (375 cells total)
grid <- tidyr::expand_grid(m = c(5, 10, 20, 40, 80), y_ratio = c(1, 2, 3, 5, 10),
                           B = c(60, 120, 240), rep = seq_len(REPS))
cli::cli_h1("Pool sensitivity: {nrow(grid)} cells")

results <- grid |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  purrr::pmap_dfr(function(m, y_ratio, B, rep, .row) {
    if (.row %% 50 == 0) cli::cli_alert_info("pool {(.row)}/{nrow(grid)}")
    t0 <- Sys.time()
    r <- lcda_gr(g, variant = 1, B = B, m = m, y = y_ratio * m, seed = 1000 * rep + .row)
    tibble::tibble(m = m, y_ratio = y_ratio, B = B, rep = rep,
                   best_Q = r$best$Q, best_H = r$best$H, best_iter = r$best$iter,
                   H_decisive_pct = 100 * length(r$lex_decisive_iters) / B,
                   secs = as.numeric(Sys.time() - t0, units = "secs"))
  })

meta <- make_meta(
  title = "Pool/refresh hyperparameter sensitivity for LCDA-GR",
  description = "Best Q and lex-decisive fraction over a (m, y/m, B) grid on a 5-block SBM(300). Tests whether the paper's m=20, y=3m, B=150 is justified and whether the H tie-break is mostly ornamental.",
  network_source = paste(SRC_SBM, "-- SBM(300, 5 blocks, p_in=0.12, p_out=0.02)"),
  seed = 2026, reps = REPS, network = "SBM(300, 5 blocks, p_in=0.12, p_out=0.02)",
  limitations = c("Single synthetic SBM; conclusions about m/y may differ on real networks.",
                  "B capped at 240 for runtime."))
save_dataset("pool_sensitivity", results, meta)
