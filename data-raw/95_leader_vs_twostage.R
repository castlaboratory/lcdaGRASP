# data-raw/95_leader_vs_twostage.R
#
# Fase 1.9b -- THE THESIS TEST. Is the JOINTLY-detected LCDA leader a better
# representative of its community than a TWO-STAGE pipeline (detect communities
# with a strong method, then pick the most central node per community)?
#
# This is the honest competitor for "joint leader+community detection": not
# Leiden alone (Leiden gives no leaders), but Leiden/ECG + per-community
# top-centrality leader. LCDA seeds each community around its leader and grows
# it; the two-stage picks the centrality hub of an independently-detected
# community -- which may sit at a community's boundary. We test whether LCDA's
# seed-then-grow leader is better CENTRED in its community.
#
# Metrics (independent of NCE, to avoid circularity):
#   coverage  -- fraction of the community within <=2 hops of its leader;
#   within_inf-- expected fraction of the community activated by an IC cascade
#                seeded at the leader (measured inside the community).
# Methods: LCDA-GR (joint) vs Leiden+eigen-leader vs ECG+eigen-leader (two-stage).
#
# -> inst/extdata/leader_vs_twostage.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
NETS <- strsplit(Sys.getenv("LVT_NETS", "karate,dolphins,football,polbooks,polblogs"), ",")[[1]]
MC   <- as.integer(Sys.getenv("LVT_MC", "500"))
PIC  <- as.numeric(Sys.getenv("LVT_PIC", "0.10"))

# leader of each community = highest eigenvector-centrality member (two-stage rule)
central_leaders <- function(g, mem) {
  ev <- igraph::eigen_centrality(g)$vector
  vapply(sort(unique(mem)), function(c) { idx <- which(mem == c); idx[which.max(ev[idx])] },
         integer(1))
}

# per-community coverage (<=2 hops) and within-community IC influence
leader_metrics <- function(g, mem, leaders, nb, mc, p) {
  comms <- sort(unique(mem))
  rows <- lapply(seq_along(comms), function(k) {
    c <- comms[k]; members <- which(mem == c); L <- leaders[k]
    if (length(members) <= 1) return(tibble::tibble(comm = c, size = length(members), coverage = 1, within_inf = 1))
    d <- as.numeric(igraph::distances(g, v = L, to = members))
    cover <- mean(is.finite(d) & d <= 2)
    # IC cascade from L; fraction of THIS community activated (avg over mc)
    inmask <- logical(igraph::vcount(g)); inmask[members] <- TRUE
    act <- vapply(seq_len(mc), function(i) {
      active <- logical(igraph::vcount(g)); active[L] <- TRUE; frontier <- L
      while (length(frontier)) {
        nw <- integer(0)
        for (u in frontier) { cand <- nb[[u]]; cand <- cand[!active[cand]]
          if (length(cand)) { hit <- cand[stats::runif(length(cand)) < p]; if (length(hit)) { active[hit] <- TRUE; nw <- c(nw, hit) } } }
        frontier <- unique(nw)
      }
      sum(active & inmask)
    }, numeric(1))
    tibble::tibble(comm = c, size = length(members),
                   coverage = cover, within_inf = mean(act) / length(members))
  })
  dplyr::bind_rows(rows)
}

per_comm <- list(); summ <- list()
for (net in NETS) {
  g <- load_benchmark(net); nb <- nbr_list(g)
  cli::cli_h2("{net} (n={igraph::vcount(g)})")
  # joint
  fit <- lcda_gr(g, variant = 1, B = 150, centrality = "eigen", similarity = "hpi", seed = 1)
  m_lcda <- leader_metrics(g, fit$best$membership, fit$best$leaders, nb, MC, PIC) |> dplyr::mutate(method = "LCDA-GR (joint)")
  # two-stage: Leiden + central leader
  cl <- igraph::cluster_leiden(g, objective_function = "modularity"); mem_l <- as.integer(igraph::membership(cl))
  m_lei <- leader_metrics(g, mem_l, central_leaders(g, mem_l), nb, MC, PIC) |> dplyr::mutate(method = "Leiden + central")
  # two-stage: ECG + central leader (best recovery method)
  mem_e <- ecg_membership(g)
  rows_e <- if (!is.null(mem_e)) leader_metrics(g, mem_e, central_leaders(g, mem_e), nb, MC, PIC) |> dplyr::mutate(method = "ECG + central") else NULL
  pc <- dplyr::bind_rows(m_lcda, m_lei, rows_e) |> dplyr::mutate(network = net)
  per_comm[[net]] <- pc
  summ[[net]] <- pc |> dplyr::group_by(network, method) |>
    dplyr::summarise(n_comms = dplyr::n(),
                     coverage = weighted.mean(coverage, size),
                     within_inf = weighted.mean(within_inf, size), .groups = "drop")
  print(as.data.frame(summ[[net]]), row.names = FALSE, digits = 3)
}

per_comm_tbl <- dplyr::bind_rows(per_comm); summ_tbl <- dplyr::bind_rows(summ)
meta <- make_meta(
  title = "Leader quality: joint LCDA vs two-stage (Leiden/ECG + central node)",
  description = paste("Per-community leader coverage (<=2 hops) and within-community IC influence,",
                      "comparing LCDA's jointly-detected leader against the two-stage pipeline",
                      "(strong community detector + per-community top-eigenvector leader)."),
  seed = 20260529, mc = MC, p_ic = PIC,
  limitations = c("Leaders for two-stage use eigenvector centrality (the paper's recommended centrality).",
                  "Communities differ across methods, so this is an end-to-end pipeline comparison, not a fixed-community ablation.",
                  "IC simulator is pure R."))
save_dataset("leader_vs_twostage", list(per_comm = per_comm_tbl, summary = summ_tbl), meta)
cli::cli_alert_success("Thesis test done.")
