# data-raw/100_stability_pool.R
#
# Fase 1 / improvement #2: STABILITY-GUIDED pool size. The Reactive update
# concentrates probability on the best parameter pairs (Proposition 6); as a
# consequence the GRASP pool's co-association matrix should CONVERGE as the pool
# grows, giving a principled stopping rule instead of an arbitrary B. We measure,
# at increasing pool sizes B:
#   - coassoc_drift: mean |C_B - C_{B_prev}| (how much co-association still moves);
#   - consensus_stability: NMI(consensus_B, consensus_{B_prev}) (partition settled?);
#   - recovery: NMI(consensus_B, planted truth) (does quality plateau?).
# A stability rule "stop when drift < eps (or stability > 0.99)" then selects B.
#
# -> inst/extdata/stability_pool.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
MUS  <- as.numeric(strsplit(Sys.getenv("STB_MUS", "0.3"), ",")[[1]])
NGR  <- as.integer(Sys.getenv("STB_NGRAPH", "5"))
BMAX <- as.integer(Sys.getenv("STB_BMAX", "128"))
CHK  <- as.integer(strsplit(Sys.getenv("STB_CHK", "2,4,8,16,32,64,128"), ",")[[1]])
W_MIN <- 0.05

# collect BMAX pool memberships once, then evaluate consensus at each checkpoint
pool_memberships <- function(g, Bmax) {
  csr <- lcdaGRASP:::.as_csr(g)
  vapply(seq_len(Bmax), function(b) {
    sol <- lcda_construct(csr, stats::runif(1, 0.1, 0.9), stats::runif(1, 0.1, 0.5),
                          variant = 1, centrality = "eigen", similarity = "hpi")
    sol <- lcda_repair(csr, sol); sol <- lcda_local_search(csr, sol)
    sol$membership
  }, integer(csr$n))                                   # n x Bmax matrix
}
consensus_from <- function(g, el, coassoc) {
  core <- igraph::coreness(g); in2 <- core[el[, 1]] >= 2 & core[el[, 2]] >= 2
  w <- ifelse(in2, W_MIN + (1 - W_MIN) * coassoc, W_MIN)
  as.integer(igraph::membership(igraph::cluster_louvain(g, weights = w)))
}

rows <- list()
for (mu in MUS) for (r in seq_len(NGR)) {
  gen <- lfr_generate(500, mu, seed = 17000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  g <- gen$graph; truth <- gen$truth
  el <- igraph::as_edgelist(g, names = FALSE)
  mems <- pool_memberships(g, BMAX)
  prev_co <- NULL; prev_cons <- NULL
  for (B in CHK) {
    sub <- mems[, seq_len(B), drop = FALSE]
    coassoc <- rowMeans(matrix(sub[el[, 1], ] == sub[el[, 2], ], nrow = nrow(el)))
    cons <- consensus_from(g, el, coassoc)
    drift <- if (is.null(prev_co)) NA_real_ else mean(abs(coassoc - prev_co))
    stab  <- if (is.null(prev_cons)) NA_real_ else nmi(cons, prev_cons)
    rows[[length(rows)+1]] <- tibble::tibble(
      mu = mu, graph = r, B = B, coassoc_drift = drift,
      consensus_stability = stab, recovery = nmi(cons, truth), k = length(unique(cons)))
    prev_co <- coassoc; prev_cons <- cons
  }
  cli::cli_alert("mu={mu} g={r} done")
}
res <- dplyr::bind_rows(rows)
summ <- res |> dplyr::group_by(B) |>
  dplyr::summarise(coassoc_drift = mean(coassoc_drift, na.rm = TRUE),
                   consensus_stability = mean(consensus_stability, na.rm = TRUE),
                   recovery = mean(recovery), .groups = "drop")
cat("\n=== pool convergence vs B (mean over graphs) ===\n")
print(as.data.frame(summ), row.names = FALSE, digits = 4)

meta <- make_meta(
  title = "Stability-guided pool size: co-association convergence of the GRASP ensemble",
  description = paste("As the GRASP pool grows, the co-association matrix and the consensus partition",
                      "converge (drift -> 0, stability -> 1) and recovery plateaus -- a principled",
                      "stopping rule, consistent with the Reactive concentration of Proposition 6."),
  seed = 20260529, mus = MUS, ngraph = NGR, bmax = BMAX, checkpoints = CHK,
  limitations = c("LFR n=500; single re-clusterer (Louvain); convergence shown empirically, not bounded."))
save_dataset("stability_pool", res, meta)
cli::cli_alert_success("Stability-guided pool experiment done.")
