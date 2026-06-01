# data-raw/36_realnet_amazon.R
#
# Fase 1.4 (Wald), external validity: a LARGE REAL network with ground-truth
# communities, complementing the synthetic LFR large-scale study. Amazon-Computers
# (Shchur et al. 2018): a co-purchase graph, 13,752 nodes, ~246k edges, 10
# product-category labels -- undirected, denser (avg degree ~36) than every
# paper benchmark, so it also probes the method outside the sparse-LFR regime.
#
# Recovery (NMI/ARI vs the category labels), modularity Q, and wall-clock time
# for LCDA-ECG (eigenvector) vs Louvain / Leiden / ECG. Honest expectation:
# product categories are not pure modular communities, so NMI is moderate for
# every method -- the point is real-data external validity, not an easy win.
#
# -> inst/extdata/realnet_amazon.rds  (list: results, meta)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
B_ECG   <- as.integer(Sys.getenv("RA_B_ECG", "32"))
B_GR    <- as.integer(Sys.getenv("RA_B_GR", "30"))
ECG_ENS <- as.integer(Sys.getenv("RA_ECG_ENS", "16"))

timed <- function(expr) { t0 <- Sys.time(); m <- force(expr)
  list(mem = m, secs = as.numeric(Sys.time() - t0, units = "secs")) }

# ---- fetch + build the largest connected component -------------------------
py <- .venv_python()
if (is.na(py)) { cli::cli_alert_warning("No .venv-lfr Python; skipping Amazon-Computers."); }
helper <- file.path(PKG_DIR, "data-raw", "amazon_computers.py")
ef <- file.path(tempdir(), "amazon_edges.txt")
lf <- file.path(tempdir(), "amazon_labels.txt")
cache <- file.path(tempdir(), "amazon_computers.npz")   # never under the repo tree
st <- if (is.na(py)) 1L else suppressWarnings(system2(py,
  c(helper, "--edges", ef, "--labels", lf, "--cache", cache),
  stdout = FALSE, stderr = FALSE))

if (is.na(py) || st != 0 || !file.exists(ef) || !file.exists(lf)) {
  cli::cli_alert_warning("Amazon-Computers fetch failed; dataset not (re)generated.")
} else {
  el <- as.matrix(utils::read.table(ef))
  truth_all <- scan(lf, what = integer(), quiet = TRUE)
  n_all <- length(truth_all)
  g <- igraph::graph_from_edgelist(el, directed = FALSE)
  if (igraph::vcount(g) < n_all) g <- igraph::add_vertices(g, n_all - igraph::vcount(g))
  igraph::V(g)$truth <- truth_all
  g <- igraph::simplify(g)
  comp <- igraph::components(g)
  g <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  truth <- igraph::V(g)$truth
  cli::cli_h1("Amazon-Computers LCC: n={vcount(g)} m={ecount(g)} avg_deg={round(2*ecount(g)/vcount(g),1)} classes={length(unique(truth))}")

  # ---- run the recovery-figure method set (eigenvector centrality) ---------
  res <- list()
  res$Louvain    <- timed(as.integer(igraph::membership(igraph::cluster_louvain(g))))
  res$Leiden     <- timed(as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity"))))
  ecg_m          <- ecg_membership(g, ens = ECG_ENS)
  res$ECG        <- list(mem = if (is.null(ecg_m)) NA else ecg_m, secs = NA_real_)
  res$`LCDA-GR`  <- timed(lcda_gr(g, B = B_GR, centrality = "eigen", similarity = "hpi", seed = 1)$best$membership)
  res$`LCDA-ECG` <- timed(lcda_ecg(g, B = B_ECG, centrality = "eigen", similarity = "hpi", seed = 1)$membership)

  # Q computed uniformly via igraph modularity for an apples-to-apples comparison.
  results <- dplyr::bind_rows(lapply(names(res), function(m) {
    mem <- res[[m]]$mem
    ok <- length(mem) == igraph::vcount(g)
    tibble::tibble(
      method = m,
      Q   = if (ok) igraph::modularity(g, mem) else NA_real_,
      NMI = if (ok) nmi(mem, truth) else NA_real_,
      ARI = if (ok) adjusted_rand(mem, truth) else NA_real_,
      n_comms = if (ok) length(unique(mem)) else NA_integer_,
      secs = res[[m]]$secs)
  }))
  print(as.data.frame(results), row.names = FALSE, digits = 3)

  meta <- make_meta(
    title = "Large real network with ground truth: Amazon-Computers (recovery + runtime)",
    description = paste("Recovery (NMI/ARI vs 10 product-category labels), modularity Q and",
                        "wall-clock time for LCDA-ECG (eigenvector) vs Louvain/Leiden/ECG on the",
                        "largest connected component of Amazon-Computers (Shchur et al. 2018).",
                        "External-validity complement to the synthetic LFR large-scale study."),
    seed = 20260529, n = igraph::vcount(g), m_edges = igraph::ecount(g),
    classes = length(unique(truth)), b_ecg = B_ECG, b_gr = B_GR, ecg_ens = ECG_ENS,
    network_source = "Amazon-Computers, gnn-benchmark npz (Shchur, Mumme, Bojchevski & Gunnemann 2018)",
    limitations = c(
      "Product-category labels are not pure modular communities, so NMI/ARI are moderate for ALL methods (expected; this is real-data external validity, not an easy benchmark).",
      "Single graph (one real network), single seed for the metaheuristics.",
      "Eigenvector centrality only (betweenness/closeness are O(nm), impractical here); n above the dense affinity-matrix comfort zone is left to future work.",
      "Times are single-threaded R/C++; the B GRASP iterations are parallelisable in principle (not measured)."))
  save_dataset("realnet_amazon", list(eval = results), meta)
  cli::cli_alert_success("Amazon-Computers evaluation done.")
}
