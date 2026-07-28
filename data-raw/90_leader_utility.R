# data-raw/90_leader_utility.R
#
# Fase 1.9 -- DOWNSTREAM UTILITY OF THE DETECTED LEADERS (anchor of the
# leader+community reframing, Fase 5 of PLANO_PUBLICACAO.md).
#
# WHY. The paper currently evaluates leaders only through the NCE score, which
# is circular: NCE is itself part of the contribution. To claim the leaders
# are *useful* we must show they help a downstream task that does not depend on
# NCE. The canonical task for "influential nodes" is information spreading.
#
# HYPOTHESIS (honest, Box). LCDA returns exactly one leader per community, so
# the leader set is spread across the network by construction. Standard
# seed-selection heuristics (top-k by a centrality) tend to pick clustered
# hubs that are redundant for spreading. We therefore expect the LCDA leader
# set to be a STRONG, CHEAP seed set -- competitive with or better than
# centrality heuristics of the same size -- *not* to beat greedy influence
# maximisation (IM), which is the gold-standard upper reference. The selling
# point is "a diversified, structurally-grounded seed set obtained jointly
# with the community structure, for free".
#
# DESIGN.
#   Networks      : the 5 benchmarks (load_benchmark()).
#   Seed size     : g = number of communities/leaders found by LCDA-GR (so the
#                   comparison is like-for-like; every strategy gets g seeds).
#   Strategies    : lcda      -- the LCDA-GR leaders (variant 1, HPI+eigen)
#                   degree    -- top-g by degree
#                   eigen     -- top-g by eigenvector centrality
#                   between   -- top-g by betweenness
#                   kcore     -- top-g by coreness (k-core)
#                   leaderrank-- top-g by LeaderRank (Lu et al. 2011)
#                   random    -- g random nodes (averaged; null baseline)
#                   greedy_im -- greedy IM under IC (upper reference; small nets)
#   Spread models : Independent Cascade (IC, prob p) and Linear Threshold (LT).
#                   Two models = a von Neumann robustness check.
#   Probabilities : p in {0.01, 0.05, 0.10} and the degree-normalised
#                   "weighted cascade" p_uv = 1/deg(v) variant.
#   Metric        : expected spread sigma(S) = E[ |activated| / n ], estimated
#                   by Monte Carlo; report mean + bootstrap/MC CI (Tukey).
#   Statistics    : compare the MC spread distributions across strategies with
#                   the package's permutation_test()/kruskal_then_permute()
#                   (reuses the paper's own machinery).
#   Identifiability (secondary, ties to NCE_local): fraction of each community
#                   reachable from its leader within 1-2 hops -- how well the
#                   leader "covers" its community.
#
# OUTPUT -> inst/extdata/leader_utility.rds  (per-(network,model,p,strategy)
#           spread distribution + summary + identifiability).
#
# This script is the DESIGN + a correct reference implementation. The IC/LT
# simulators are pure R here for clarity; a Rcpp kernel is the productionisation
# (Fase 2). Scale is controlled by env vars so it can be smoke-tested cheaply.

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
MC      <- as.integer(Sys.getenv("LU_MC", "2000"))     # Monte Carlo runs per seed set
NRAND   <- as.integer(Sys.getenv("LU_NRAND", "30"))    # random-seed redraws
PROBS   <- as.numeric(strsplit(Sys.getenv("LU_PROBS", "0.05,0.10"), ",")[[1]])
NETS    <- strsplit(Sys.getenv("LU_NETS", "karate,dolphins,football,polbooks"), ",")[[1]]
DO_GREEDY <- as.logical(Sys.getenv("LU_GREEDY", "TRUE"))

# nbr_list() and ic_once() now live in 00_helpers.R (shared with 30_lfr_*; the
# single place where a Rcpp kernel will replace the pure-R simulator, item 2.9).

# ---- Linear Threshold: one realisation (uniform random thresholds) ---------
lt_once <- function(nb, deg, seeds) {
  n <- length(nb)
  thr <- stats::runif(n)                  # node threshold ~ U(0,1)
  active <- logical(n); active[seeds] <- TRUE
  repeat {
    changed <- FALSE
    inactive <- which(!active)
    for (v in inactive) {
      cand <- nb[[v]]
      if (!length(cand) || deg[v] == 0) next
      share <- sum(active[cand]) / deg[v]   # equal edge weights 1/deg(v)
      if (share >= thr[v]) { active[v] <- TRUE; changed <- TRUE }
    }
    if (!changed) break
  }
  sum(active)
}

spread_dist <- function(nb, deg, seeds, model, p, mc) {
  vapply(seq_len(mc), function(i)
    if (model == "IC") ic_once(nb, seeds, p) else lt_once(nb, deg, seeds),
    numeric(1)) / length(nb)
}

# ---- LeaderRank (Lu, Zhang, Yeung, Zhou 2011) ------------------------------
# Add a ground node connected both ways to every node; PageRank on the
# augmented graph; rank original nodes by score. (igraph PageRank as the solver.)
leaderrank <- function(g) {
  n <- igraph::vcount(g)
  el <- igraph::as_edgelist(g, names = FALSE)
  gnd <- n + 1L
  aug_el <- rbind(el, el[, 2:1],                 # undirected -> both directions
                  cbind(gnd, seq_len(n)), cbind(seq_len(n), gnd))
  ag <- igraph::graph_from_edgelist(aug_el, directed = TRUE)
  pr <- igraph::page_rank(ag, directed = TRUE)$vector
  pr[seq_len(n)]
}

# ---- greedy influence maximisation under IC (upper reference) --------------
greedy_im <- function(nb, deg, p, g_size, mc) {
  n <- length(nb); chosen <- integer(0); pool <- seq_len(n)
  for (s in seq_len(g_size)) {
    gains <- vapply(setdiff(pool, chosen), function(v)
      mean(spread_dist(nb, deg, c(chosen, v), "IC", p, max(50L, mc %/% 8L))),
      numeric(1))
    chosen <- c(chosen, setdiff(pool, chosen)[which.max(gains)])
  }
  chosen
}

# ---- per-network experiment ------------------------------------------------
run_net <- function(net) {
  g <- load_benchmark(net); n <- igraph::vcount(g)
  nb <- nbr_list(g); deg <- igraph::degree(g)
  cli::cli_h2("{net} (n={n}, m={igraph::ecount(g)})")

  fit <- lcda_gr(g, variant = 1, B = 150, centrality = "eigen",
                 similarity = "hpi", seed = 1)
  leaders <- fit$best$leaders; g_size <- length(leaders)
  membership <- fit$best$membership
  cli::cli_alert_info("LCDA-GR: {g_size} leaders/communities, Q={round(fit$best$Q,4)}")

  cent_seeds <- list(
    lcda       = leaders,
    degree     = order(deg, decreasing = TRUE)[seq_len(g_size)],
    eigen      = order(igraph::eigen_centrality(g)$vector, decreasing = TRUE)[seq_len(g_size)],
    between    = order(igraph::betweenness(g), decreasing = TRUE)[seq_len(g_size)],
    kcore      = order(igraph::coreness(g), decreasing = TRUE)[seq_len(g_size)],
    leaderrank = order(leaderrank(g), decreasing = TRUE)[seq_len(g_size)]
  )

  rows <- list()
  for (model in c("IC", "LT")) {
    ps <- if (model == "IC") PROBS else NA_real_   # LT has no single p
    for (p in ps) {
      # deterministic strategies
      for (nm in names(cent_seeds)) {
        d <- spread_dist(nb, deg, cent_seeds[[nm]], model, p, MC)
        rows[[length(rows) + 1]] <- tibble::tibble(
          network = net, model = model, p = p, strategy = nm,
          g_size = g_size, mean = mean(d), sd = sd(d),
          lo = quantile(d, .025), hi = quantile(d, .975))
      }
      # random baseline: NRAND redraws, pooled
      dr <- unlist(lapply(seq_len(NRAND), function(j) {
        s <- sample.int(n, g_size); spread_dist(nb, deg, s, model, p, max(50L, MC %/% NRAND)) }))
      rows[[length(rows) + 1]] <- tibble::tibble(
        network = net, model = model, p = p, strategy = "random",
        g_size = g_size, mean = mean(dr), sd = sd(dr),
        lo = quantile(dr, .025), hi = quantile(dr, .975))
      # greedy IM upper reference (IC only, small nets)
      if (DO_GREEDY && model == "IC" && n <= 200) {
        s <- greedy_im(nb, deg, p, g_size, MC)
        d <- spread_dist(nb, deg, s, model, p, MC)
        rows[[length(rows) + 1]] <- tibble::tibble(
          network = net, model = model, p = p, strategy = "greedy_im",
          g_size = g_size, mean = mean(d), sd = sd(d),
          lo = quantile(d, .025), hi = quantile(d, .975))
      }
    }
  }

  # identifiability: fraction of each community within <=2 hops of its leader
  ident <- vapply(seq_along(leaders), function(i) {
    L <- leaders[i]; members <- which(membership == membership[L])
    if (length(members) <= 1) return(1)
    d <- igraph::distances(g, v = L, to = members)
    mean(is.finite(d) & d <= 2)
  }, numeric(1))

  list(spread = dplyr::bind_rows(rows),
       identifiability = tibble::tibble(network = net, leader = leaders,
                                        comm_size = as.integer(table(membership)[as.character(membership[leaders])]),
                                        frac_within_2hops = ident))
}

res <- lapply(NETS, run_net)
spread_tbl <- dplyr::bind_rows(lapply(res, `[[`, "spread"))
ident_tbl  <- dplyr::bind_rows(lapply(res, `[[`, "identifiability"))

meta <- make_meta(
  title = "Leader downstream utility: IC/LT spreading from detected leaders vs seed-selection baselines",
  description = paste("Expected spread (fraction activated) under Independent Cascade and Linear",
                      "Threshold, seeding from LCDA-GR leaders vs degree/eigen/betweenness/k-core/",
                      "LeaderRank/random, plus greedy IM as an upper reference on small nets.",
                      "Tests whether the community-grounded leaders are a useful seed set (Fase 1.9)."),
  network_source = src_newman(NETS),
  seed = 20260529, nrep = MC,
  limitations = c(
    "IC/LT are pure-R reference implementations; a Rcpp kernel is the productionisation.",
    "Greedy IM uses a reduced MC budget and is restricted to n<=200 (cost).",
    "LeaderRank uses igraph PageRank on the ground-node-augmented graph as the solver.",
    "Seed size fixed to the number of LCDA communities for a like-for-like comparison.",
    "PolBlogs omitted by default (LU_NETS) for runtime; enable explicitly."))

save_dataset("leader_utility", list(spread = spread_tbl, identifiability = ident_tbl), meta)
cli::cli_alert_success("Fase 1.9 leader-utility experiment done.")
