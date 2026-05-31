# data-raw/97_lcda_ecg.R
#
# Fase 1 / algorithm improvement: LCDA-ECG -- ENSEMBLE CONSENSUS over the GRASP
# pool. The GRASP/Reactive loop already produces many diverse partitions but
# keeps only the best-by-(Q,H); the rest of the diversification is discarded.
# ECG (Poulin & Theberge 2019) beats Louvain/Leiden on LFR recovery by turning
# an ENSEMBLE of partitions into co-association edge weights and re-clustering.
# We apply that idea to the LCDA pool, then assign one leader per consensus
# community -- closing the recovery gap WHILE keeping the joint leader output.
#
# Pipeline:
#   1. run B LCDA constructions (randomised alpha, construct+repair+local_search),
#      collect the B membership vectors;
#   2. edge co-association w(u,v) = fraction of the B partitions with u,v together;
#   3. ECG floor: 2-core edges get w_min + (1-w_min)*coassoc, others w_min;
#   4. consensus = modularity clustering on the reweighted graph;
#   5. leader per consensus community = highest eigenvector-centrality member.
#
# Tested on canonical LFR: does LCDA-ECG reach Leiden/ECG recovery (NMI)?
#
# -> inst/extdata/lcda_ecg.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

# Uses the canonical exported lcdaGRASP::lcda_ecg() (no local re-definition).

# ---- evaluate on canonical LFR ---------------------------------------------
set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("ECG_MUS", "0.2,0.3,0.4,0.5"), ",")[[1]])
NGRAPH <- as.integer(Sys.getenv("ECG_NGRAPH", "5"))
N      <- as.integer(Sys.getenv("ECG_N", "500"))
BENS   <- as.integer(Sys.getenv("ECG_B", "64"))

rows <- list()
for (mu in MUS) for (r in seq_len(NGRAPH)) {
  gen <- lfr_generate(N, mu, seed = 11000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  g <- gen$graph; truth <- gen$truth
  m_gr   <- lcda_gr(g, variant = 1, B = 150, centrality = "eigen", similarity = "hpi", seed = 1)$best$membership
  m_ecg  <- lcda_ecg(g, B = BENS, seed = 1)$membership
  m_lei  <- as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity")))
  m_pecg <- ecg_membership(g)
  one <- function(method, mem) tibble::tibble(mu = mu, graph = r, method = method,
            NMI = nmi(mem, truth), ARI = adjusted_rand(mem, truth),
            Q = igraph::modularity(g, mem), k = length(unique(mem)))
  rows[[length(rows)+1]] <- one("LCDA-GR", m_gr)
  rows[[length(rows)+1]] <- one("LCDA-ECG", m_ecg)
  rows[[length(rows)+1]] <- one("Leiden", m_lei)
  if (!is.null(m_pecg)) rows[[length(rows)+1]] <- one("ECG (igraph)", m_pecg)
  cli::cli_alert("mu={mu} g={r}: LCDA-GR={round(nmi(m_gr,truth),3)} LCDA-ECG={round(nmi(m_ecg,truth),3)} Leiden={round(nmi(m_lei,truth),3)}")
}
res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(mu, method) |>
  dplyr::summarise(NMI = mean(NMI), ARI = mean(ARI), Q = mean(Q), k = mean(k), .groups = "drop")
cat("\n=== LCDA-ECG vs baselines (NMI vs planted LFR truth) ===\n")
print(as.data.frame(summ |> dplyr::select(mu, method, NMI) |> tidyr::pivot_wider(names_from = method, values_from = NMI)), row.names = FALSE, digits = 3)

meta <- make_meta(
  title = "LCDA-ECG: ensemble consensus over the GRASP pool (recovery on canonical LFR)",
  description = paste("Applies ECG-style co-association consensus to the LCDA GRASP pool, then assigns",
                      "one leader per consensus community. Tests whether consensus closes the LFR",
                      "recovery gap to Leiden/ECG while keeping the joint leader output."),
  seed = 20260529, mus = MUS, ngraph = NGRAPH, n = N, b_ensemble = BENS,
  limitations = c("Consensus members use construct+repair+full local search (the natural GRASP pool).",
                  "Leaders assigned post-hoc by eigenvector centrality within each consensus community.",
                  "n=500; pure-R IC not involved here (recovery only)."))
save_dataset("lcda_ecg", list(per_graph = res, summary = summ), meta)
cli::cli_alert_success("LCDA-ECG evaluation done.")
