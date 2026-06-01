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
#'   \item{`vnmi_nprime`}{VNMI candidate-subset size `n'` sweep: modularity and
#'     local-search time vs `n'` across density regimes.}
#'   \item{`degeneracy`}{Near-best partition counts and pairwise NMI.}
#'   \item{`eda_replicates`}{Replicated `(Q, H)` for exploratory analysis.}
#'   \item{`nce_alternatives`}{Global vs community-conditioned NCE leaders.}
#'   \item{`lcda_ecg`, `overlap_lcda_ecg`, `stability_pool`,
#'     `consensus_leader_test`}{LCDA-ECG ensemble consensus: LFR recovery,
#'     overlapping/bridge nodes, pool-stability stopping rule, and the
#'     consensus-vs-central leader comparison.}
#'   \item{`largescale`}{Recovery and runtime up to n = 5e4 (synthetic LFR).}
#'   \item{`realnet_amazon`, `realnet_coauthor`}{Large real networks with
#'     ground truth: modularity Q and recovery (NMI/ARI vs labels) on
#'     Amazon-Computers (co-purchase) and Coauthor-Physics (co-authorship).}
#'   \item{`openalex_leaders`}{Leader validation against an external citation
#'     signal on an OpenAlex co-authorship graph.}
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

#' Provenance of a shipped dataset (version, date, checksum)
#'
#' Reports the reproducibility metadata for one of the [lcda_data()] datasets:
#' the package version that generated it, the generation date, and its SHA-256
#' checksum looked up from the `SHA256SUMS` manifest shipped alongside the data.
#' The manifest stores bare file names, so verify it from inside the data
#' directory: `cd inst/extdata && shasum -a 256 -c SHA256SUMS`. Use it in
#' analyses and vignettes to pin exactly which version of a dataset a result
#' came from. The version is also enforced by the package's own tests, so a
#' shipped dataset always matches the installed `packageVersion("lcdaGRASP")`.
#'
#' @param name dataset name (without the `.rds` extension), as in [lcda_data()].
#' @return a one-row data frame with columns `dataset`, `pkg_version`,
#'   `generated_on`, and `sha256`; or `NULL` (invisibly) if the dataset is not
#'   installed.
#' @seealso [lcda_data()]
#' @examples
#' lcda_provenance("repro_summary")
#' @export
lcda_provenance <- function(name) {
  d <- lcda_data(name)
  if (is.null(d)) return(invisible(NULL))
  m <- d$meta
  field <- function(x) if (is.null(x)) NA_character_ else as.character(x)[1]
  sha <- NA_character_
  sums <- system.file("extdata", "SHA256SUMS", package = "lcdaGRASP")
  if (nzchar(sums) && file.exists(sums)) {
    L <- readLines(sums, warn = FALSE)
    hit <- L[grepl(paste0("[ /]", name, "\\.rds$"), L)]
    if (length(hit)) sha <- sub("[[:space:]].*$", "", hit[1])
  }
  data.frame(
    dataset      = name,
    pkg_version  = field(m$pkg_version),
    generated_on = field(m$generated_on),
    sha256       = sha,
    stringsAsFactors = FALSE
  )
}
