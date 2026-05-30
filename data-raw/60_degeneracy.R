# data-raw/60_degeneracy.R
#
# Good, de Montjoye & Clauset (2010): many networks have an exponentially large
# set of near-best partitions with Q within epsilon of the maximum. Audit how
# many DISTINCT near-best partitions LCDA-GRASP finds, and their pairwise NMI
# (low NMI = high degeneracy = the reported single partition is arbitrary).
#
# -> inst/extdata/degeneracy.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

graphs <- list(
  karate = igraph::make_graph("Zachary"),
  dolphins = load_benchmark("dolphins"),
  polbooks = load_benchmark("polbooks"),
  sbm_polblogs_like = igraph::sample_sbm(800,
     pref.matrix = { m <- matrix(0.005, 2, 2); diag(m) <- 0.06; m },
     block.sizes = c(400, 400)))

B <- 60L
eps_grid <- c(0.001, 0.005, 0.01, 0.02)

records <- list()
for (gname in names(graphs)) {
  g <- graphs[[gname]]
  cli::cli_alert_info("degeneracy: {gname}")
  csr <- lcdaGRASP:::.as_csr(g)
  partitions <- vector("list", B); Qs <- numeric(B)
  for (b in seq_len(B)) {
    set.seed(1000 + b)
    sol <- lcda_construct(csr, 0.1, 0.3, variant = 1, centrality = "eigen", similarity = "hpi")
    sol <- lcda_repair(csr, sol)
    sol <- lcda_local_search(csr, sol)
    partitions[[b]] <- sol$membership; Qs[b] <- sol$Q
  }
  Qstar <- max(Qs)
  for (eps in eps_grid) {
    keep <- which(Qs >= Qstar - eps)
    n_distinct <- length(unique(vapply(partitions[keep], paste, character(1), collapse = ",")))
    mean_nmi <- if (length(keep) >= 2) {
      pr <- utils::combn(keep, 2)
      mean(vapply(seq_len(ncol(pr)), function(k) nmi(partitions[[pr[1,k]]], partitions[[pr[2,k]]]), numeric(1)))
    } else NA_real_
    records[[paste0(gname, eps)]] <- tibble::tibble(
      graph = gname, B = B, epsilon = eps, Q_star = Qstar,
      n_runs_in_band = length(keep), n_distinct_partitions = n_distinct,
      mean_pairwise_NMI = mean_nmi)
  }
}

meta <- make_meta(
  title = "Modularity degeneracy audit (Good et al. 2010)",
  description = "Number of distinct near-best partitions (Q within epsilon of the max over B constructions) and their mean pairwise NMI, per network. Mean NMI well below 1 at small epsilon = the reported partition is one of many equivalent ones.",
  seed = 1000, B = B,
  limitations = c("Real networks for karate/dolphins/polbooks; an SBM proxy stands in for PolBlogs scale.",
                  "Degeneracy is assessed only among the partitions the metaheuristic actually visited."))
save_dataset("degeneracy", dplyr::bind_rows(records), meta)
