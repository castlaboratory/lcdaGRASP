# data-raw/12_blogs_table9.R
#
# E1 resolution (Option 1): regenerate Table 9 (tab:sim_blogs) with the companion
# package, replacing the duplicated LCDA-GR/LCDA-GRASP rows. The package is R with
# performance-critical kernels in C++ via Rcpp.
#
# To match the paper's Political Blogs (Table 2: n=1490, m=19090 directed edges),
# we collapse the directed graph to an undirected SIMPLE graph whose edge WEIGHTS
# are the multiplicities (sum of weights = 19090 minus 3 self-loops = 19087);
# weighted modularity on this graph equals modularity on the 19090-edge
# multigraph, so the Q scale matches the paper rather than the simple 16715-edge
# reduction. (This also exercises the weighted kernels.)
#
# -> inst/extdata/blogs_table9.rds  (per-run + summary tibbles)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
NREP  <- as.integer(Sys.getenv("T9_NREP", "5"))
cents <- c("betweenness", "closeness", "eigen")
sims  <- c("dice", "hpi", "jaccard")
algos <- c("LCDA-GRASP 1", "LCDA-GRASP 2", "LCDA-GR 1", "LCDA-GR 2")
CKPT  <- file.path(PKG_DIR, "data-raw", ".t9_ckpt"); dir.create(CKPT, showWarnings = FALSE)

# weighted Political Blogs matching the paper's edge count via multiplicities
g0 <- suppressWarnings(igraph::read_graph(file.path(BENCH_DIR, "polblogs.gml"), format = "gml"))
igraph::E(g0)$weight <- 1
g <- igraph::simplify(igraph::as_undirected(g0, mode = "each"),
                      edge.attr.comb = list(weight = "sum"))
cli::cli_h2("polblogs (weighted): n={igraph::vcount(g)} m={igraph::ecount(g)} sum(w)={sum(igraph::E(g)$weight)}")

run_one <- function(algo, cent, sim, B, seed) {
  t0 <- Sys.time()
  res <- switch(algo,
    "LCDA-GRASP 1" = lcda_grasp(g, variant = 1, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GRASP 2" = lcda_grasp(g, variant = 2, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GR 1"    = lcda_gr(g, variant = 1, B = B, centrality = cent, similarity = sim, seed = seed),
    "LCDA-GR 2"    = lcda_gr(g, variant = 2, B = B, centrality = cent, similarity = sim, seed = seed))
  tibble::tibble(algo = algo, centrality = cent, similarity = sim, B = B, seed = seed,
                 Q = res$best$Q, secs = as.numeric(Sys.time() - t0, units = "secs"))
}

grid <- expand.grid(algo = algos, centrality = cents, similarity = sims, stringsAsFactors = FALSE)
# Per-run results are fully determined by the per-rep seed (run_one(..., seed=r)),
# so cells are independent and can run in parallel without changing the output.
# T9_CORES>1 forks over the grid (mclapply); default 1 keeps the sequential path.
CORES <- as.integer(Sys.getenv("T9_CORES", "1"))
run_cell <- function(i) {
  algo <- grid$algo[i]; cent <- grid$centrality[i]; sim <- grid$similarity[i]
  cp <- file.path(CKPT, sprintf("c%02d.rds", i))
  if (file.exists(cp)) return(readRDS(cp))
  B <- if (grepl("GRASP", algo)) 50L else 150L
  df <- purrr::map_dfr(seq_len(NREP), function(r) tryCatch(run_one(algo, cent, sim, B, r), error = function(e) NULL))
  saveRDS(df, cp)
  cli::cli_alert("{algo}/{cent}/{sim} (B={B}): Q={round(mean(df$Q),4)} t={round(mean(df$secs),1)}s")
  df
}
per_run <- if (CORES > 1L) {
  parallel::mclapply(seq_len(nrow(grid)), run_cell, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), run_cell)
}
results <- dplyr::bind_rows(per_run)
summary_tbl <- results |> dplyr::group_by(algo, centrality, similarity, B) |>
  dplyr::summarise(Q_mean = mean(Q), Q_sd = sd(Q), secs_mean = mean(secs), .groups = "drop")

meta <- make_meta(
  title = "Table 9 regeneration (E1): mean Q and time on weighted Political Blogs",
  description = paste("Companion-package reproduction of Table 9 (tab:sim_blogs). Package is R with",
                      "C++/Rcpp kernels. Political Blogs collapsed to a weighted simple graph",
                      "(sum of weights = paper's directed edge count), so weighted modularity matches",
                      "the paper's scale. Replaces the duplicated LCDA-GR rows with genuine numbers."),
  seed = 20260529, nrep = NREP,
  limitations = c("Times are the package's C++/Rcpp kernels (faster than a pure-R implementation).",
                  "Closeness now uses harmonic centrality on the disconnected graph (robust fallback)."))
save_dataset("blogs_table9", list(per_run = results, summary = summary_tbl), meta)
unlink(CKPT, recursive = TRUE)
cli::cli_alert_success("Table 9 regeneration done.")
