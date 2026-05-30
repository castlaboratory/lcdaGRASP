# data-raw/80_nce.R
#
# The paper's NCE uses p_1(v) = k_v / n (a GLOBAL property). A hub spanning all
# communities can maximise it yet be a poor community leader. Compare three
# scores and measure how often the chosen leader CHANGES under the proposed
# community-conditioned variant.
#
#   I_global(v) = H_b(k_v / n)                            (paper)
#   I_local(v)  = H_b(|N(v) cap C(v)| / (|C(v)| - 1))     (proposed)
#   I_mix(v)    = 0.5 I_local + 0.5 I_global              (compromise)
#
# -> inst/extdata/nce_alternatives.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

best_leader_under <- function(members, scorer)
  members[which.max(vapply(members, scorer, numeric(1)))]

make_scorers <- function(g, membership) {
  n <- igraph::vcount(g); k <- igraph::degree(g)
  I_global <- function(v) { p <- k[v] / n
    if (p > 0 && p < 1) -(p*log2(p) + (1-p)*log2(1-p)) else (1 + log2(n)) / n }
  I_local <- function(v) {
    members <- which(membership == membership[v])
    if (length(members) <= 1) return(0)
    kin <- sum(as.integer(igraph::neighbors(g, v)) %in% members)
    p <- kin / (length(members) - 1)
    if (p > 0 && p < 1) -(p*log2(p) + (1-p)*log2(1-p)) else 0 }
  I_mix <- function(v, lambda = 0.5) lambda * I_local(v) + (1 - lambda) * I_global(v)
  list(global = I_global, local = I_local, mix = I_mix)
}

graphs <- list(
  Karate   = igraph::make_graph("Zachary"),
  Dolphins = load_benchmark("dolphins"),
  Polbooks = load_benchmark("polbooks"),
  SBM_5    = igraph::sample_sbm(300,
               pref.matrix = { m <- matrix(0.02, 5, 5); diag(m) <- 0.12; m },
               block.sizes = rep(60, 5)))

records <- list()
for (gname in names(graphs)) {
  g <- graphs[[gname]]
  for (rep in 1:10) {
    res <- lcda_grasp(g, alpha_c = 0.1, alpha_s = 0.3, B = 50,
                      centrality = "eigen", similarity = "hpi", seed = rep)
    mem <- res$best$membership; leaders_paper <- res$best$leaders
    sc <- make_scorers(g, mem)
    # Community ids are not necessarily 1..g after local search merges; iterate
    # the ids actually present (the order matches lcda_repair's leader order).
    comm_ids <- sort(unique(mem))
    leaders_local <- vapply(comm_ids,
      function(c) best_leader_under(which(mem == c), sc$local), integer(1))
    leaders_mix <- vapply(comm_ids,
      function(c) best_leader_under(which(mem == c), function(v) sc$mix(v, 0.5)), integer(1))
    records[[paste0(gname, rep)]] <- tibble::tibble(
      graph = gname, rep = rep, n_communities = length(leaders_paper),
      H_paper_global = mean(vapply(leaders_paper, sc$global, numeric(1))),
      H_paper_local  = mean(vapply(leaders_paper, sc$local,  numeric(1))),
      n_leaders_changed_local = sum(leaders_paper != leaders_local),
      pct_leaders_changed = 100 * sum(leaders_paper != leaders_local) / length(leaders_paper),
      n_leaders_changed_mix = sum(leaders_paper != leaders_mix))
  }
}

meta <- make_meta(
  title = "NCE leader-score alternatives (global vs community-conditioned)",
  description = "How often the designated leader changes when switching from the paper's global NCE to the proposed community-conditioned (local) NCE. A non-trivial change fraction means the NCE choice materially decides who is called a leader.",
  seed = 1, reps = 10,
  limitations = c("Leaders are re-picked as the in-community argmax of each score, holding the partition fixed.",
                  "Mix uses lambda=0.5; the trade-off is not optimised."))
save_dataset("nce_alternatives", dplyr::bind_rows(records), meta)
