# data-raw/55_weighted.R
#
# Fase 1.5 demonstration: weighted-graph support. We take LFR networks and add
# INFORMATIVE edge weights (intra-community edges heavier than inter-community),
# then check that using the weights IMPROVES recovery -- i.e. the weighted
# modularity machinery (now in the kernels) is not only correct (test T9) but
# useful. Compares LCDA-GR and Louvain run weighted vs. ignoring the weights.
#
# -> inst/extdata/weighted_demo.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("W_MUS", "0.3,0.4,0.5"), ",")[[1]])
NGRAPH <- as.integer(Sys.getenv("W_NGRAPH", "5"))
N      <- as.integer(Sys.getenv("W_N", "500"))

add_informative_weights <- function(g, truth) {
  el <- igraph::as_edgelist(g, names = FALSE)
  intra <- truth[el[, 1]] == truth[el[, 2]]
  w <- numeric(nrow(el))
  w[intra]  <- stats::runif(sum(intra),  1.5, 4.0)   # heavier inside communities
  w[!intra] <- stats::runif(sum(!intra), 0.1, 1.0)   # lighter across communities
  igraph::E(g)$weight <- w
  g
}

rows <- list()
for (mu in MUS) for (r in seq_len(NGRAPH)) {
  gen <- lfr_generate(N, mu, seed = 70000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  truth <- gen$truth
  gw <- add_informative_weights(gen$graph, truth)          # weighted graph
  gu <- igraph::delete_edge_attr(gw, "weight")             # same topology, no weights
  one <- function(method, weighted, mem) tibble::tibble(
    mu = mu, graph = r, method = method, weighted = weighted, NMI = nmi(mem, truth))
  rows[[length(rows)+1]] <- one("LCDA-GR", TRUE,  lcda_gr(gw, B = 60, seed = 1)$best$membership)
  rows[[length(rows)+1]] <- one("LCDA-GR", FALSE, lcda_gr(gu, B = 60, seed = 1)$best$membership)
  rows[[length(rows)+1]] <- one("Louvain", TRUE,  as.integer(igraph::membership(igraph::cluster_louvain(gw, weights = igraph::E(gw)$weight))))
  rows[[length(rows)+1]] <- one("Louvain", FALSE, as.integer(igraph::membership(igraph::cluster_louvain(gu))))
  if (r == 1) cli::cli_alert("mu={mu} done")
}
res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(mu, method, weighted) |> dplyr::summarise(NMI = mean(NMI), .groups = "drop")
cat("\n=== recovery: weighted vs unweighted (informative weights) ===\n")
print(as.data.frame(summ |> tidyr::pivot_wider(names_from = weighted, values_from = NMI,
                                               names_prefix = "w_")), row.names = FALSE, digits = 3)

meta <- make_meta(
  title = "Weighted-graph support: informative weights improve recovery",
  description = paste("LFR with intra-community edges weighted heavier than inter-community ones.",
                      "Compares LCDA-GR and Louvain run with vs without the weights. Demonstrates",
                      "that the weighted modularity kernels (validated against igraph in test T9)",
                      "are useful: using informative weights improves NMI recovery."),
  network_source = paste(SRC_LFR, "-- edge weights are synthetic (U(1.5,4) intra-community, U(0.1,1) inter-community), not measured."),
  seed = 20260529, mus = MUS, ngraph = NGRAPH, n = N,
  limitations = c("Weights are synthetic and informative by construction (best case for weighting).",
                  "Construction centrality/similarity remain structural; the objective is weighted."))
save_dataset("weighted_demo", list(recovery = res, summary = summ), meta)
cli::cli_alert_success("Weighted demo done.")
