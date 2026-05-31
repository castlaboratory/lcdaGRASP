# data-raw/00_helpers.R
#
# Shared helpers for the data-generation pipeline. Sourced by every
# data-raw/*.R script and by run_all_data.R. These scripts run the FULL
# (publication-grade) simulations once and cache the results as .rds under
# inst/extdata/, so the package vignettes can render without ever re-running
# a simulation at build time.
#
# Run from the repository root OR from the package directory; paths are
# auto-detected.

suppressPackageStartupMessages({
  library(lcdaGRASP)
  library(igraph)
  library(tibble)
  library(dplyr)
})

# ---- path discovery --------------------------------------------------------

.locate_pkg_dir <- function() {
  cands <- c("lcdaGRASP", ".", "..", "../lcdaGRASP")
  for (d in cands) {
    desc <- file.path(d, "DESCRIPTION")
    if (file.exists(desc) && any(grepl("^Package: lcdaGRASP", readLines(desc, warn = FALSE)))) {
      return(normalizePath(d))
    }
  }
  cli::cli_abort("Could not locate the lcdaGRASP package directory.")
}

PKG_DIR     <- .locate_pkg_dir()
EXTDATA_DIR <- file.path(PKG_DIR, "inst", "extdata")
dir.create(EXTDATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Benchmark GML cache: prefer the repo-root data/benchmarks, else download.
BENCH_DIR <- local({
  cands <- c(file.path(PKG_DIR, "data-raw", "benchmarks"),
             file.path(dirname(PKG_DIR), "data", "benchmarks"),
             "data-raw/benchmarks", "data/benchmarks", "../data/benchmarks")
  hit <- cands[dir.exists(cands)]
  if (length(hit)) normalizePath(hit[1]) else {
    d <- file.path(dirname(PKG_DIR), "data", "benchmarks")
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    d
  }
})

NEWMAN_URL <- "https://websites.umich.edu/~mejn/netdata"

# ---- network loaders -------------------------------------------------------

.fetch_gml <- function(name) {
  gml <- file.path(BENCH_DIR, paste0(name, ".gml"))
  if (file.exists(gml)) return(gml)
  zip <- file.path(BENCH_DIR, paste0(name, ".zip"))
  cli::cli_alert_info("Downloading {name}.zip from Newman's collection")
  utils::download.file(file.path(NEWMAN_URL, paste0(name, ".zip")),
                       zip, quiet = TRUE, mode = "wb")
  utils::unzip(zip, exdir = BENCH_DIR)
  if (!file.exists(gml)) cli::cli_abort("expected {.path {gml}} after unzip")
  gml
}

#' Load a benchmark network as an undirected simple igraph.
load_benchmark <- function(name) {
  if (name == "karate") return(igraph::make_graph("Zachary"))
  g <- suppressWarnings(igraph::read_graph(.fetch_gml(name), format = "gml"))
  g <- igraph::as_undirected(g, mode = "collapse")
  igraph::simplify(g)
}

# The five paper benchmarks, with literature reference modularity.
BENCHMARKS <- tibble::tribble(
  ~network,   ~n,   ~Q_ref,  ~kind,
  "karate",   34L,  0.4198,  "exact optimum",
  "dolphins", 62L,  0.5285,  "best known",
  "football", 115L, 0.6046,  "best known",
  "polbooks", 105L, 0.5272,  "best known",
  "polblogs", 1490L,0.4271,  "best known"
)

# ---- metadata (Box lens: always state context/source/sample/limitations) ---

make_meta <- function(title, description, seed, ..., limitations = character()) {
  list(
    title         = title,
    description   = description,
    generated_on  = as.character(Sys.Date()),
    generated_at  = as.character(Sys.time()),
    pkg_version   = as.character(utils::packageVersion("lcdaGRASP")),
    R_version     = R.version.string,
    igraph_version= as.character(utils::packageVersion("igraph")),
    seed          = seed,
    host          = unname(Sys.info()[["sysname"]]),
    network_source= "Newman netdata collection (websites.umich.edu/~mejn) + igraph::make_graph('Zachary')",
    reference     = "Ospina, Silva, Matos Junior, Leite & Ochi (2026), preprint submitted to Knowledge-Based Systems.",
    limitations   = limitations,
    extra         = list(...)
  )
}

#' Save a {results, meta} bundle to inst/extdata/<name>.rds.
save_dataset <- function(name, results, meta) {
  stopifnot(tibble::is_tibble(results) || is.data.frame(results) || is.list(results))
  path <- file.path(EXTDATA_DIR, paste0(name, ".rds"))
  saveRDS(list(results = results, meta = meta), path, compress = "xz")
  cli::cli_alert_success("Saved {.path {path}}  ({nrow_safe(results)} rows)")
  invisible(path)
}

nrow_safe <- function(x) if (is.null(nrow(x))) length(x) else nrow(x)

# ---- partition-comparison metrics (no external deps) -----------------------

adjusted_rand <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  cmb <- function(v) sum(choose(v, 2))
  a <- cmb(rowSums(tab)); b <- cmb(colSums(tab)); ab <- cmb(as.numeric(tab))
  expected <- (a * b) / choose(n, 2); maxr <- (a + b) / 2
  (ab - expected) / (maxr - expected)
}

nmi <- function(a, b) {
  tab <- table(a, b); N <- sum(tab)
  pa <- rowSums(tab) / N; pb <- colSums(tab) / N
  Ha <- -sum(ifelse(pa > 0, pa * log2(pa), 0))
  Hb <- -sum(ifelse(pb > 0, pb * log2(pb), 0))
  pab <- tab / N
  I <- sum(ifelse(pab > 0, pab * log2(pab / outer(pa, pb)), 0))
  if (Ha + Hb == 0) return(1)
  2 * I / (Ha + Hb)
}

# ---- diffusion simulators (shared by 30_lfr and 90_leader_utility) ---------
# Pure-R reference implementations of Independent Cascade. A Rcpp kernel is the
# productionisation (PLANO_PUBLICACAO.md, item 2.9); this is the single place to
# swap it in. Neighbour lists are 1-based igraph vertex ids.

nbr_list <- function(g) lapply(igraph::as_adj_list(g, mode = "all"), as.integer)

# One Independent-Cascade realisation from `seeds`; returns the number activated.
ic_once <- function(nb, seeds, p) {
  active <- logical(length(nb)); active[seeds] <- TRUE
  frontier <- seeds
  while (length(frontier)) {
    newly <- integer(0)
    for (u in frontier) {
      cand <- nb[[u]]
      if (!length(cand)) next
      cand <- cand[!active[cand]]
      if (!length(cand)) next
      hit <- cand[stats::runif(length(cand)) < p]
      if (length(hit)) { active[hit] <- TRUE; newly <- c(newly, hit) }
    }
    frontier <- unique(newly)
  }
  sum(active)
}

# Expected IC spread (fraction of nodes activated), Monte Carlo mean.
ic_spread <- function(nb, seeds, p, mc) {
  mean(vapply(seq_len(mc), function(i) ic_once(nb, seeds, p), numeric(1))) / length(nb)
}

# Locate the .venv-lfr Python interpreter (the venv lives at the repo root).
.venv_python <- function() {
  cands <- c(file.path(dirname(PKG_DIR), ".venv-lfr", "bin", "python"),
             file.path(PKG_DIR, ".venv-lfr", "bin", "python"),
             file.path(dirname(PKG_DIR), "lcdaGRASP", ".venv-lfr", "bin", "python"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

# Generate one canonical LFR graph via the .venv-lfr Python generator.
# Returns list(graph = igraph, truth = integer membership) or NULL on failure.
lfr_generate <- function(n, mu, seed, avg_degree = 12, min_comm = 20, max_comm = 60) {
  py  <- .venv_python(); if (is.na(py)) return(NULL)
  gen <- file.path(PKG_DIR, "data-raw", "lfr_generate.py")
  ef <- tempfile(fileext = ".edges"); tf <- tempfile(fileext = ".truth")
  on.exit(unlink(c(ef, tf)), add = TRUE)
  st <- suppressWarnings(system2(py, c(gen, "--n", n, "--mu", mu, "--seed", seed,
    "--avg-degree", avg_degree, "--min-comm", min_comm, "--max-comm", max_comm,
    "--edges", ef, "--truth", tf), stdout = FALSE, stderr = FALSE))
  if (st != 0 || !file.exists(ef) || !file.exists(tf)) return(NULL)
  el <- as.matrix(utils::read.table(ef))
  truth <- scan(tf, what = integer(), quiet = TRUE)
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE))
  # ensure all n nodes exist even if some are isolated after simplify
  if (igraph::vcount(g) < n) g <- igraph::add_vertices(g, n - igraph::vcount(g))
  list(graph = g, truth = truth)
}

# ECG (Ensemble Clustering for Graphs, Poulin & Theberge 2019) via partition-igraph
# in .venv-lfr. Returns a 1-based membership vector, or NULL on failure.
ecg_membership <- function(g, ens = 16) {
  py <- .venv_python(); if (is.na(py)) return(NULL)
  helper <- file.path(PKG_DIR, "data-raw", "ecg_membership.py")
  n <- igraph::vcount(g)
  ef <- tempfile(fileext = ".edges"); of <- tempfile(fileext = ".mem")
  on.exit(unlink(c(ef, of)), add = TRUE)
  el <- igraph::as_edgelist(g, names = FALSE)
  utils::write.table(el, ef, row.names = FALSE, col.names = FALSE)
  st <- suppressWarnings(system2(py, c(helper, "--n", n, "--edges", ef, "--out", of,
                                       "--ens", ens), stdout = FALSE, stderr = FALSE))
  if (st != 0 || !file.exists(of)) return(NULL)
  scan(of, what = integer(), quiet = TRUE)
}

cli::cli_alert_info("data-raw helpers loaded. PKG_DIR={.path {PKG_DIR}}")
cli::cli_alert_info("EXTDATA_DIR={.path {EXTDATA_DIR}}  BENCH_DIR={.path {BENCH_DIR}}")
