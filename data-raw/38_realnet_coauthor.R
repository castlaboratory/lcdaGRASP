# data-raw/38_realnet_coauthor.R
#
# Real-network external validity on a LEADER-MEANINGFUL domain: Coauthor-Physics
# (Shchur et al. 2018) -- nodes are authors, edges are co-authorship, labels are
# research fields. A "central author within a field" is a far more natural
# reading of a community leader than a hub product, so this replaces the
# Amazon-Computers co-purchase graph as the paper's real community-recovery
# anchor (Amazon remains available as a denser stress test, lcda_data).
#
# LCC ~ 34k nodes, 5 fields. Recovery (NMI/ARI vs field labels), modularity Q,
# and wall-clock time for LCDA-ECG (eigenvector) vs Louvain/Leiden/ECG/LCDA-GR.
# As elsewhere, field labels are not pure modular communities, so NMI is moderate
# for every method -- this is real-data external validity, not an easy win.
#
# -> inst/extdata/realnet_coauthor.rds  (list: eval, meta)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
B_ECG   <- as.integer(Sys.getenv("RC_B_ECG", "32"))
B_GR    <- as.integer(Sys.getenv("RC_B_GR", "30"))
ECG_ENS <- as.integer(Sys.getenv("RC_ECG_ENS", "16"))
DATASET <- Sys.getenv("RC_DATASET", "coauthor_physics")

timed <- function(expr) { t0 <- Sys.time(); m <- force(expr)
  list(mem = m, secs = as.numeric(Sys.time() - t0, units = "secs")) }

py <- .venv_python()
helper <- file.path(PKG_DIR, "data-raw", "realnet_npz.py")
ef <- file.path(tempdir(), "co_edges.txt"); lf <- file.path(tempdir(), "co_labels.txt")
cache <- file.path(tempdir(), paste0(DATASET, ".npz"))
st <- if (is.na(py)) 1L else suppressWarnings(system2(py,
  c(helper, "--dataset", DATASET, "--edges", ef, "--labels", lf, "--cache", cache),
  stdout = FALSE, stderr = FALSE))

if (is.na(py) || st != 0 || !file.exists(ef) || !file.exists(lf)) {
  cli::cli_alert_warning("Coauthor fetch failed; dataset not (re)generated.")
} else {
  el <- as.matrix(utils::read.table(ef))
  truth_all <- scan(lf, what = integer(), quiet = TRUE)
  g <- igraph::graph_from_edgelist(el, directed = FALSE)
  if (igraph::vcount(g) < length(truth_all))
    g <- igraph::add_vertices(g, length(truth_all) - igraph::vcount(g))
  igraph::V(g)$truth <- truth_all
  g <- igraph::simplify(g)
  cm <- igraph::components(g)
  g <- igraph::induced_subgraph(g, which(cm$membership == which.max(cm$csize)))
  truth <- igraph::V(g)$truth
  cli::cli_h1("{DATASET} LCC: n={vcount(g)} m={ecount(g)} avg_deg={round(2*ecount(g)/vcount(g),1)} fields={length(unique(truth))}")

  res <- list()
  res$Louvain    <- timed(as.integer(igraph::membership(igraph::cluster_louvain(g))))
  res$Leiden     <- timed(as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity"))))
  ecg_m          <- ecg_membership(g, ens = ECG_ENS)
  res$ECG        <- list(mem = if (is.null(ecg_m)) NA else ecg_m, secs = NA_real_)
  res$`LCDA-GR`  <- timed(lcda_gr(g, B = B_GR, centrality = "eigen", similarity = "hpi", seed = 1)$best$membership)
  res$`LCDA-ECG` <- timed(lcda_ecg(g, B = B_ECG, centrality = "eigen", similarity = "hpi", seed = 1)$membership)

  results <- dplyr::bind_rows(lapply(names(res), function(m) {
    mem <- res[[m]]$mem; ok <- length(mem) == igraph::vcount(g)
    tibble::tibble(method = m,
      Q   = if (ok) igraph::modularity(g, mem) else NA_real_,
      NMI = if (ok) nmi(mem, truth) else NA_real_,
      ARI = if (ok) adjusted_rand(mem, truth) else NA_real_,
      n_comms = if (ok) length(unique(mem)) else NA_integer_,
      secs = res[[m]]$secs)
  }))
  print(as.data.frame(results), row.names = FALSE, digits = 3)

  meta <- make_meta(
    title = "Real coauthorship network with ground truth: Coauthor-Physics (recovery + runtime)",
    description = paste("Recovery (NMI/ARI vs 5 research-field labels), modularity Q and wall-clock",
                        "time for LCDA-ECG (eigenvector) vs Louvain/Leiden/ECG on the largest",
                        "connected component of Coauthor-Physics (Shchur et al. 2018). A",
                        "leader-meaningful real network (authors / co-authorship / fields)."),
    seed = 20260529, n = igraph::vcount(g), m_edges = igraph::ecount(g),
    fields = length(unique(truth)), b_ecg = B_ECG, b_gr = B_GR, ecg_ens = ECG_ENS,
    network_source = "Coauthor-Physics, gnn-benchmark ms_academic_phy.npz (Shchur et al. 2018)",
    limitations = c(
      "Field labels are not pure modular communities, so NMI/ARI are moderate for ALL methods (real-data external validity, not an easy benchmark).",
      "Community recovery only; leader quality against an external impact signal is validated separately (OpenAlex coauthorship study).",
      "Eigenvector centrality only (betweenness/closeness are O(nm)); single graph, single metaheuristic seed.",
      "Times are single-threaded R/C++; the B GRASP iterations are parallelisable in principle (not measured)."))
  save_dataset("realnet_coauthor", list(eval = results), meta)
  cli::cli_alert_success("Coauthor-Physics evaluation done.")
}
