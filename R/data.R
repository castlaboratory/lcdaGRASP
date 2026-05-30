# R/data.R
#
# Accessors for the precomputed simulation datasets shipped under
# inst/extdata/. These back the package vignettes so they render without
# re-running any simulation. Each .rds is a list(results = <tibble>,
# meta = <list>); `meta` records the date, seed, package versions, network
# source, and limitations (always state context and limits).

#' Precomputed simulation datasets
#'
#' The package ships the results of its publication-grade simulation panel as
#' compressed `.rds` files under `inst/extdata/`. The vignettes read these so
#' they never re-run a simulation at build time. Regenerate them with the
#' scripts in `data-raw/` (see `data-raw/run_all_data.R`).
#'
#' Available datasets (name -> contents):
#' \describe{
#'   \item{`repro_benchmarks`, `repro_summary`}{Benchmark reproduction on the
#'     five networks vs Louvain/Leiden and literature best-known Q.}
#'   \item{`doe_screening`, `doe_rsm`}{Design-of-experiments runs: factorial
#'     screening and a central composite design in `(alpha_c, alpha_s)`.}
#'   \item{`lfr_robustness`}{Q and ARI across the LFR-like mixing sweep.}
#'   \item{`prop6_summary`, `prop6_trajectories`}{Reactive-update concentration.}
#'   \item{`pool_sensitivity`}{`(m, y, B)` grid for LCDA-GR.}
#'   \item{`degeneracy`}{Near-best partition counts and pairwise NMI.}
#'   \item{`eda_replicates`}{Replicated `(Q, H)` for exploratory analysis.}
#'   \item{`nce_alternatives`}{Global vs community-conditioned NCE leaders.}
#'   \item{`lcda_ecg`, `overlap_lcda_ecg`, `stability_pool`,
#'     `consensus_leader_test`}{LCDA-ECG ensemble consensus: LFR recovery,
#'     overlapping/bridge nodes, pool-stability stopping rule, and the
#'     consensus-vs-central leader comparison.}
#'   \item{`largescale`}{Recovery and runtime up to n = 5e4.}
#'   \item{`gnn_baseline`}{Graph-auto-encoder baseline (true-k and auto-k) vs
#'     classical methods on LFR.}
#'   \item{`weighted_demo`}{Weighted-graph demonstration (weights aid recovery).}
#'   \item{`blogs_table9`, `blogs_timing`}{Weighted Political Blogs reproduction
#'     (mean/max Q, timings).}
#'   \item{`doe_lfr_recovery`}{DoE screening/RSM scored by LFR recovery.}
#'   \item{`leader_utility`, `leader_vs_twostage`,
#'     `leader_vs_twostage_lfr`}{Downstream leader utility (IC/LT spread,
#'     coverage) vs two-stage detect-then-centrality pipelines.}
#' }
#'
#' @param name dataset name (without the `.rds` extension). If `NULL`
#'   (default), returns a tibble listing the available datasets.
#' @return If `name` is `NULL`, a tibble of available datasets. Otherwise the
#'   `list(results, meta)` bundle for that dataset, or `NULL` with a message if
#'   it has not been generated yet.
#' @examples
#' lcda_data()                       # list what is available
#' \dontrun{
#' d <- lcda_data("repro_summary")
#' d$meta$limitations
#' head(d$results)
#' }
#' @name lcda_data
#' @export
lcda_data <- function(name = NULL) {
  dir <- system.file("extdata", package = "lcdaGRASP")
  files <- list.files(dir, pattern = "\\.rds$", full.names = FALSE)
  if (is.null(name)) {
    return(tibble::tibble(dataset = sub("\\.rds$", "", files)))
  }
  path <- system.file("extdata", paste0(name, ".rds"), package = "lcdaGRASP")
  if (!nzchar(path) || !file.exists(path)) {
    cli::cli_alert_warning(
      "Dataset {.val {name}} not found. Generate it with {.path data-raw/run_all_data.R}.")
    return(invisible(NULL))
  }
  readRDS(path)
}
