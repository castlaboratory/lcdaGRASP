# data-raw/45_gnn_baseline.R
#
# Fase 1.6 (Wald): an HONEST GNN comparison. The paper argued GNN methods are
# inapplicable to attribute-free graphs; rather than deflect, we run a canonical
# Graph Auto-Encoder (Kipf & Welling 2016) on the same topology with synthesised
# STRUCTURAL features (degree, clustering, k-core, eigenvector centrality), then
# KMeans on the embedding. The GAE is even given the true number of communities
# k (generous). We compare recovery (NMI) against Louvain/Leiden/LCDA-ECG on the
# canonical LFR sweep.
#
# -> inst/extdata/gnn_baseline.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })
options(lcda_ecg_define_only = TRUE)
source(file.path(PKG_DIR, "data-raw", "97_lcda_ecg.R"))
options(lcda_ecg_define_only = FALSE)

gae_membership <- function(g, k, epochs = 200, seed = 1) {
  py <- .venv_python(); if (is.na(py)) return(NULL)
  helper <- file.path(PKG_DIR, "data-raw", "gae_baseline.py")
  ef <- tempfile(fileext = ".edges"); of <- tempfile(fileext = ".mem")
  on.exit(unlink(c(ef, of)), add = TRUE)
  utils::write.table(igraph::as_edgelist(g, names = FALSE), ef, row.names = FALSE, col.names = FALSE)
  st <- suppressWarnings(system2(py, c(helper, "--edges", ef, "--n", igraph::vcount(g),
    "--k", k, "--out", of, "--epochs", epochs, "--seed", seed), stdout = FALSE, stderr = FALSE))
  if (st != 0 || !file.exists(of)) return(NULL)
  scan(of, what = integer(), quiet = TRUE)
}

set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("GNN_MUS", "0.1,0.2,0.3,0.4,0.5"), ",")[[1]])
NGRAPH <- as.integer(Sys.getenv("GNN_NGRAPH", "5"))
N      <- as.integer(Sys.getenv("GNN_N", "500"))

rows <- list()
for (mu in MUS) for (r in seq_len(NGRAPH)) {
  gen <- lfr_generate(N, mu, seed = 60000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  g <- gen$graph; truth <- gen$truth; k <- length(unique(truth))
  one <- function(method, mem) tibble::tibble(mu = mu, graph = r, method = method,
            NMI = if (length(mem) > 1) nmi(mem, truth) else NA_real_)
  rows[[length(rows)+1]] <- one("GAE (true k)", gae_membership(g, k, seed = 1))
  rows[[length(rows)+1]] <- one("GAE (auto-k)", gae_membership(g, 0L, seed = 1))
  rows[[length(rows)+1]] <- one("Louvain", as.integer(igraph::membership(igraph::cluster_louvain(g))))
  rows[[length(rows)+1]] <- one("Leiden", as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity"))))
  rows[[length(rows)+1]] <- one("LCDA-ECG", lcda_ecg(g, B = 32, seed = 1)$membership)
  if (r == 1) cli::cli_alert("mu={mu} done")
}
res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(mu, method) |> dplyr::summarise(NMI = mean(NMI), .groups = "drop")
cat("\n=== NMI vs planted truth (GNN vs classical) ===\n")
print(as.data.frame(summ |> tidyr::pivot_wider(names_from = method, values_from = NMI)), row.names = FALSE, digits = 3)

meta <- make_meta(
  title = "Honest GNN baseline: Graph Auto-Encoder with structural features vs classical methods (LFR)",
  description = paste("Canonical GAE (Kipf & Welling 2016) on attribute-free LFR using synthesised",
                      "structural node features + KMeans (given the true k), vs Louvain/Leiden/LCDA-ECG.",
                      "Tests the paper's claim that GNNs are ill-suited to attribute-free community detection."),
  seed = 20260529, mus = MUS, ngraph = NGRAPH, n = N,
  limitations = c("Two GAE rows: 'true k' (KMeans given the planted community count, generous) and",
                  "'auto-k' (k chosen by silhouette over a candidate range, the realistic no-ground-truth setting).",
                  "Structural features only (no real node attributes, which is the regime GNNs are built for).",
                  "A 2-layer GAE; deeper/contrastive GNNs may differ, but the attribute-free regime is the point."))
save_dataset("gnn_baseline", list(recovery = res, summary = summ), meta)
cli::cli_alert_success("GNN baseline done.")
