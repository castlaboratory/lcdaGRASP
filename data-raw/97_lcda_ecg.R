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

# ---- LCDA-ECG: consensus over the GRASP pool -------------------------------
lcda_ecg <- function(g, B = 64, w_min = 0.05,
                     alpha_c_range = c(0.1, 0.9), alpha_s_range = c(0.1, 0.5),
                     centrality = "eigen", similarity = "hpi", seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  csr <- lcdaGRASP:::.as_csr(g)
  n   <- igraph::vcount(g)
  el  <- igraph::as_edgelist(g, names = FALSE)            # m x 2, 1-based
  coassoc <- numeric(nrow(el))
  lead_count <- integer(n)                                # leader-designation frequency
  for (b in seq_len(B)) {
    ac <- stats::runif(1, alpha_c_range[1], alpha_c_range[2])
    as_ <- stats::runif(1, alpha_s_range[1], alpha_s_range[2])
    sol <- lcda_construct(csr, ac, as_, variant = 1, centrality = centrality, similarity = similarity)
    sol <- lcda_repair(csr, sol)
    sol <- lcda_local_search(csr, sol)
    mem <- sol$membership
    coassoc <- coassoc + (mem[el[, 1]] == mem[el[, 2]])
    lead_count[sol$leaders] <- lead_count[sol$leaders] + 1L   # LCDA's own leader designations
  }
  coassoc <- coassoc / B
  core <- igraph::coreness(g)
  in2 <- core[el[, 1]] >= 2 & core[el[, 2]] >= 2
  w <- ifelse(in2, w_min + (1 - w_min) * coassoc, w_min)  # ECG-style edge weights
  cons <- igraph::cluster_louvain(g, weights = w)
  mem_c <- as.integer(igraph::membership(cons))
  ev <- igraph::eigen_centrality(g)$vector
  comms <- sort(unique(mem_c))
  # (a) post-hoc centrality leader (= the two-stage rule)
  leaders_cent <- vapply(comms, function(c) { i <- which(mem_c == c); i[which.max(ev[i])] }, integer(1))
  # (b) CONSENSUS leader = most-frequently-designated leader in the pool within
  #     the community (tie-break by eigenvector centrality). This uses LCDA's own
  #     ensemble of leader designations -- a signal ECG+centrality cannot access.
  leaders_freq <- vapply(comms, function(c) {
    i <- which(mem_c == c)
    score <- lead_count[i] + ev[i] / (max(ev) + 1)        # ev as a small tie-breaker
    i[which.max(score)]
  }, integer(1))
  # (#2) node-level confidence: mean co-association of a node's incident edges
  #      (how robustly the node sits inside its community across the ensemble)
  conf <- numeric(n)
  for (e in seq_len(nrow(el))) { conf[el[e, 1]] <- conf[el[e, 1]] + coassoc[e]; conf[el[e, 2]] <- conf[el[e, 2]] + coassoc[e] }
  deg <- igraph::degree(g); conf <- ifelse(deg > 0, conf / deg, 0)
  list(membership = mem_c, leaders = leaders_freq, leaders_cent = leaders_cent,
       lead_count = lead_count, confidence = conf, Q = igraph::modularity(g, mem_c))
}

# ---- evaluate on canonical LFR ---------------------------------------------
# Guard: when this file is sourced only to reuse lcda_ecg() (e.g. by 98_), skip
# the evaluation/save below.
if (!isTRUE(getOption("lcda_ecg_define_only"))) {
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
}  # end guard
