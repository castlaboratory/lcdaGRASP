# data-raw/37_vnmi_nprime.R
#
# von Neumann (b): sensitivity of the VNMI candidate-subset size n' on solution
# quality and cost. The paper's default n'=300 was a cost-based choice; this
# characterises how the final modularity Q and the local-search time depend on
# n', across DENSITY regimes (the critique notes a dense graph may need a larger
# n', a sparse one a smaller). Paired design: per (graph, rep) a single
# construction+repair is reused, and only n' varies in the local search, so the
# curve isolates the n' effect.
#
# Graphs (all n>300 so VNMI is active): sparse LFR (avg deg 12), dense LFR
# (avg deg 50), and the real Political Blogs network.
#
# -> inst/extdata/vnmi_nprime.rds  (list: sweep, meta)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
NPRIMES <- as.integer(strsplit(Sys.getenv("VN_NPRIMES", "100,200,300,500,1000"), ",")[[1]])
REPS    <- as.integer(Sys.getenv("VN_REPS", "5"))
N_LFR   <- as.integer(Sys.getenv("VN_N_LFR", "1000"))
AC <- 0.1; AS <- 0.3                      # paper defaults

# largest connected component (closeness/eigenvector want a connected graph)
lcc <- function(g) {
  cm <- igraph::components(g)
  igraph::induced_subgraph(g, which(cm$membership == which.max(cm$csize)))
}

# one graph for a given regime + rep; returns igraph or NULL
make_graph_regime <- function(regime, rep) {
  if (regime == "PolBlogs (real, deg~26)") return(lcc(load_benchmark("polblogs")))
  gen <- if (grepl("dense", regime))
    lfr_generate(N_LFR, mu = 0.3, seed = 7000 + rep, avg_degree = 60,
                 max_degree = 300, min_comm = 80, max_comm = 300)
  else
    lfr_generate(N_LFR, mu = 0.3, seed = 7000 + rep, avg_degree = 12,
                 min_comm = 30, max_comm = 100)
  if (is.null(gen)) NULL else lcc(gen$graph)
}

regimes <- c("LFR sparse (deg 12)", "LFR dense (deg 50)", "PolBlogs (real, deg~26)")
rows <- list()
for (regime in regimes) {
  pol_g <- NULL
  for (rep in seq_len(REPS)) {
    g <- if (regime == "PolBlogs (real, deg~26)") {
      if (is.null(pol_g)) pol_g <- make_graph_regime(regime, rep); pol_g  # same graph, vary construction seed
    } else make_graph_regime(regime, rep)
    if (is.null(g) || igraph::vcount(g) <= 300) next
    csr <- as_csr(g)
    set.seed(1000 * rep + 17)
    base <- lcda_repair(lcda_construct(csr, AC, AS, centrality = "eigen", similarity = "hpi"),
                        csr = csr, centrality = "eigen")
    for (np in NPRIMES) {
      if (np > igraph::vcount(g)) next
      t0 <- Sys.time()
      sol <- lcda_local_search(csr, base, centrality = "eigen", vnmi_n_prime = np)
      secs <- as.numeric(Sys.time() - t0, units = "secs")
      rows[[length(rows) + 1]] <- tibble::tibble(
        regime = regime, n = igraph::vcount(g),
        avg_deg = round(2 * igraph::ecount(g) / igraph::vcount(g), 1),
        rep = rep, n_prime = np,
        Q = igraph::modularity(g, sol$membership), secs = secs)
    }
    if (rep %% 2 == 1) cli::cli_alert_info("{regime}: rep {rep}/{REPS} (n={vcount(g)})")
  }
}
sweep <- dplyr::bind_rows(rows)
print(sweep |> group_by(regime, n_prime) |>
        summarise(Q = mean(Q), secs = mean(secs), .groups = "drop") |> as.data.frame(),
      digits = 3, row.names = FALSE)

meta <- make_meta(
  title = "VNMI candidate-subset size n': sensitivity of Q and cost across density regimes",
  description = paste("Paired sweep of the VNMI subset size n' (fixed construction+repair per",
                      "(graph, rep); only n' varies in local search). Q and local-search time vs",
                      "n' for sparse LFR (deg 12), dense LFR (deg 50) and Political Blogs.",
                      "Characterises the cost-based default n'=300: how far Q has converged and",
                      "how time grows with n'."),
  seed = 20260529, n_primes = NPRIMES, reps = REPS, n_lfr = N_LFR, alpha_c = AC, alpha_s = AS,
  limitations = c(
    "Single mu=0.3 for the LFR regimes; alpha_c/alpha_s at the paper defaults.",
    "Paired on one construction per (graph, rep): isolates the local-search n' effect, not the full GRASP pool.",
    "PolBlogs is one real graph; variation there is over construction seeds, not graphs."))
save_dataset("vnmi_nprime", list(sweep = sweep), meta)
cli::cli_alert_success("VNMI n' sweep done.")
