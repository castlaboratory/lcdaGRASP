# data-raw/96_leader_vs_twostage_lfr.R
#
# Fase 1.9b control -- the GRANULARITY-MATCHED version of the thesis test.
# On PolBlogs the joint-vs-two-stage leader comparison was confounded: LCDA made
# 12 large communities, Leiden/ECG ~270 tiny ones, so per-community coverage is
# trivially high for the fragmented partitions. On canonical LFR the planted
# community count is fixed, and at low mixing all methods recover a SIMILAR
# number of communities (high NMI), so the comparison is fair: any difference in
# leader coverage / within-community influence reflects leader CENTRING, not
# granularity.
#
# Methods: LCDA-GR (joint) vs Leiden+central vs ECG+central. We also record the
# community count per method to demonstrate matched granularity.
#
# -> inst/extdata/leader_vs_twostage_lfr.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("LVTL_MUS", "0.1,0.2,0.3"), ",")[[1]])  # recoverable regime
NGRAPH <- as.integer(Sys.getenv("LVTL_NGRAPH", "5"))
N      <- as.integer(Sys.getenv("LVTL_N", "500"))
MC     <- as.integer(Sys.getenv("LVTL_MC", "400"))
PIC    <- as.numeric(Sys.getenv("LVTL_PIC", "0.10"))

central_leaders <- function(g, mem) {
  ev <- igraph::eigen_centrality(g)$vector
  vapply(sort(unique(mem)), function(c) { i <- which(mem == c); i[which.max(ev[i])] }, integer(1))
}
leader_metrics <- function(g, mem, leaders, nb, mc, p) {
  comms <- sort(unique(mem)); nn <- igraph::vcount(g)
  ws <- vapply(seq_along(comms), function(k) {
    members <- which(mem == comms[k]); L <- leaders[k]
    if (length(members) <= 1) return(c(1, 1, length(members)))
    d <- as.numeric(igraph::distances(g, v = L, to = members))
    cover <- mean(is.finite(d) & d <= 2)
    inmask <- logical(nn); inmask[members] <- TRUE
    act <- vapply(seq_len(mc), function(i) {
      active <- logical(nn); active[L] <- TRUE; fr <- L
      while (length(fr)) { nw <- integer(0)
        for (u in fr) { cand <- nb[[u]]; cand <- cand[!active[cand]]
          if (length(cand)) { hit <- cand[stats::runif(length(cand)) < p]; if (length(hit)) { active[hit] <- TRUE; nw <- c(nw, hit) } } }
        fr <- unique(nw) }
      sum(active & inmask) }, numeric(1))
    c(cover, mean(act) / length(members), length(members))
  }, numeric(3))
  c(coverage = weighted.mean(ws[1, ], ws[3, ]),
    within_inf = weighted.mean(ws[2, ], ws[3, ]), n_comms = ncol(ws))
}

rows <- list()
for (mu in MUS) for (r in seq_len(NGRAPH)) {
  gen <- lfr_generate(N, mu, seed = 9000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  g <- gen$graph; nb <- nbr_list(g)
  fit <- lcda_gr(g, variant = 1, B = 150, centrality = "eigen", similarity = "hpi", seed = 1)
  mem_l <- as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity")))
  mem_e <- ecg_membership(g)
  add <- function(method, m) { x <- leader_metrics(g, m$mem, m$lead, nb, MC, PIC)
    tibble::tibble(mu = mu, graph = r, method = method, coverage = x["coverage"],
                   within_inf = x["within_inf"], n_comms = x["n_comms"]) }
  rows[[length(rows)+1]] <- add("LCDA-GR (joint)", list(mem = fit$best$membership, lead = fit$best$leaders))
  rows[[length(rows)+1]] <- add("Leiden + central", list(mem = mem_l, lead = central_leaders(g, mem_l)))
  if (!is.null(mem_e)) rows[[length(rows)+1]] <- add("ECG + central", list(mem = mem_e, lead = central_leaders(g, mem_e)))
  cli::cli_alert("mu={mu} g={r}: LCDA k={length(fit$best$leaders)}, Leiden k={length(unique(mem_l))}")
}

res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(mu, method) |>
  dplyr::summarise(coverage = mean(coverage), within_inf = mean(within_inf),
                   n_comms = mean(n_comms), .groups = "drop")
cat("\n=== granularity-matched leader comparison on LFR ===\n")
print(as.data.frame(summ), row.names = FALSE, digits = 3)

meta <- make_meta(
  title = "Granularity-matched leader comparison on canonical LFR (joint vs two-stage)",
  description = paste("Leader coverage (<=2 hops) and within-community IC influence for LCDA (joint)",
                      "vs Leiden/ECG + per-community central leader, on canonical LFR at low mixing",
                      "where all methods recover a similar number of communities -- removes the",
                      "granularity confound seen on PolBlogs."),
  seed = 20260529, mus = MUS, ngraph = NGRAPH, n = N, mc = MC, p_ic = PIC,
  limitations = c("Low-mixing regime only (where recovery -- hence matched granularity -- holds).",
                  "IC simulator is pure R; two-stage leaders use eigenvector centrality."))
save_dataset("leader_vs_twostage_lfr", list(per_graph = res, summary = summ), meta)
cli::cli_alert_success("LFR granularity-matched thesis control done.")
