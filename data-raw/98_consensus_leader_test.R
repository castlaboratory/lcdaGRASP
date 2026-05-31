# data-raw/98_consensus_leader_test.R
#
# Does the CONSENSUS-DERIVED leader (most-frequently-designated across the GRASP
# pool) make a better community representative than the post-hoc top-centrality
# leader? The cleanest test: compare the two leader rules ON THE SAME LCDA-ECG
# consensus communities -- this isolates leader selection (no granularity or
# partition-quality confound). Leiden+central and ECG+central included for
# context. Metrics (NCE-free): coverage (<=2 hops) and within-community IC
# influence. Granularity-matched LFR (low mixing).
#
# -> inst/extdata/consensus_leader_test.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })
# lcda_ecg() comes from the installed lcdaGRASP package (loaded via 00_helpers).

set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("CLT_MUS", "0.1,0.2,0.3"), ",")[[1]])
NGRAPH <- as.integer(Sys.getenv("CLT_NGRAPH", "5"))
N      <- as.integer(Sys.getenv("CLT_N", "500"))
MC     <- as.integer(Sys.getenv("CLT_MC", "400"))
PIC    <- as.numeric(Sys.getenv("CLT_PIC", "0.10"))
BENS   <- as.integer(Sys.getenv("CLT_B", "64"))

central_leaders <- function(g, mem) { ev <- igraph::eigen_centrality(g)$vector
  vapply(sort(unique(mem)), function(c) { i <- which(mem == c); i[which.max(ev[i])] }, integer(1)) }
leader_metrics <- function(g, mem, leaders, nb, mc, p) {
  comms <- sort(unique(mem)); nn <- igraph::vcount(g)
  ws <- vapply(seq_along(comms), function(k) {
    members <- which(mem == comms[k]); L <- leaders[k]
    if (length(members) <= 1) return(c(1, 1, length(members)))
    d <- as.numeric(igraph::distances(g, v = L, to = members))
    inmask <- logical(nn); inmask[members] <- TRUE
    act <- vapply(seq_len(mc), function(i) { active <- logical(nn); active[L] <- TRUE; fr <- L
      while (length(fr)) { nw <- integer(0)
        for (u in fr) { cand <- nb[[u]]; cand <- cand[!active[cand]]
          if (length(cand)) { hit <- cand[stats::runif(length(cand)) < p]; if (length(hit)) { active[hit] <- TRUE; nw <- c(nw, hit) } } }
        fr <- unique(nw) }
      sum(active & inmask) }, numeric(1))
    c(mean(is.finite(d) & d <= 2), mean(act) / length(members), length(members))
  }, numeric(3))
  c(coverage = weighted.mean(ws[1, ], ws[3, ]), within_inf = weighted.mean(ws[2, ], ws[3, ]))
}

rows <- list()
for (mu in MUS) for (r in seq_len(NGRAPH)) {
  gen <- lfr_generate(N, mu, seed = 13000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  g <- gen$graph; nb <- nbr_list(g)
  ec <- lcda_ecg(g, B = BENS, seed = 1)
  mem_l <- as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity")))
  mem_e <- ecg_membership(g)
  add <- function(method, mem, lead) { x <- leader_metrics(g, mem, lead, nb, MC, PIC)
    tibble::tibble(mu = mu, graph = r, method = method, coverage = x["coverage"], within_inf = x["within_inf"]) }
  rows[[length(rows)+1]] <- add("LCDA-ECG consensus-leader", ec$membership, ec$leaders)       # freq
  rows[[length(rows)+1]] <- add("LCDA-ECG central-leader",   ec$membership, ec$leaders_central)  # same comms!
  rows[[length(rows)+1]] <- add("Leiden + central", mem_l, central_leaders(g, mem_l))
  if (!is.null(mem_e)) rows[[length(rows)+1]] <- add("ECG + central", mem_e, central_leaders(g, mem_e))
  cli::cli_alert("mu={mu} g={r} done")
}
res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(mu, method) |>
  dplyr::summarise(coverage = mean(coverage), within_inf = mean(within_inf), .groups = "drop")
cat("\n=== leader quality on consensus communities (and baselines) ===\n")
print(as.data.frame(summ), row.names = FALSE, digits = 3)
# paired, SAME communities: consensus-leader vs central-leader
w <- res |> dplyr::filter(grepl("LCDA-ECG", method)) |>
  dplyr::select(mu, graph, method, within_inf) |>
  tidyr::pivot_wider(names_from = method, values_from = within_inf)
cat(sprintf("\nConsensus vs central leader (SAME communities), within-influence: mean diff %+.3f, consensus wins %d/%d\n",
  mean(w$`LCDA-ECG consensus-leader` - w$`LCDA-ECG central-leader`),
  sum(w$`LCDA-ECG consensus-leader` > w$`LCDA-ECG central-leader`), nrow(w)))

meta <- make_meta(
  title = "Consensus-derived leader vs central leader on LCDA-ECG consensus communities",
  description = paste("Isolates the leader-selection rule: consensus (most-frequently-designated in",
                      "the GRASP pool) vs post-hoc top-eigenvector, on the SAME LCDA-ECG communities,",
                      "plus Leiden/ECG+central for context. Granularity-matched LFR."),
  seed = 20260529, mus = MUS, ngraph = NGRAPH, n = N, mc = MC, b_ensemble = BENS,
  limitations = c("Low-mixing regime; pure-R IC; n=500."))
save_dataset("consensus_leader_test", list(per_graph = res, summary = summ), meta)
cli::cli_alert_success("Consensus-leader test done.")
