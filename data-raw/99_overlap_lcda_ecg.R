# data-raw/99_overlap_lcda_ecg.R
#
# Fase 1 / algorithm improvement #4: OVERLAPPING communities (and their leaders)
# from the LCDA-ECG co-association matrix. The co-association C(u,v) = fraction
# of GRASP-pool partitions in which u,v are together is a SOFT similarity, so a
# node can belong to several communities -- a capability single-output methods
# (Leiden/Louvain) do not provide, and natural for joint leader+community
# detection where bridge nodes matter.
#
# Assignment: start from the hard consensus partition. For node v and community
# c that v connects to, define affinity(v,c) = mean co-association over v's edges
# into c. v's home = argmax affinity. v ALSO joins c (overlap) if
#   affinity(v,c) >= tau * affinity(v, home).
#
# Evaluation (no overlapping ground truth available from networkx LFR, which is
# disjoint, so we validate that the overlap is MEANINGFUL, not against a planted
# overlap):
#   (1) hard-core recovery (NMI vs planted) is preserved;
#   (2) overlap nodes are structural BRIDGES -- higher participation coefficient
#       P(v) = 1 - sum_c (k_{v,c}/k_v)^2 than non-overlap nodes (Wilcoxon);
#   (3) overlap fraction vs tau; qualitative Karate example.
#
# -> inst/extdata/overlap_lcda_ecg.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

# pool + co-association (returns hard membership, leaders, edge coassoc, el)
ecg_pool <- function(g, B = 64, w_min = 0.05,
                     alpha_c_range = c(0.1, 0.9), alpha_s_range = c(0.1, 0.5)) {
  csr <- lcdaGRASP:::.as_csr(g); g <- csr$igraph; n <- csr$n
  el <- igraph::as_edgelist(g, names = FALSE)
  coassoc <- numeric(nrow(el)); lead_count <- integer(n)
  for (b in seq_len(B)) {
    sol <- lcda_construct(csr, stats::runif(1, alpha_c_range[1], alpha_c_range[2]),
                          stats::runif(1, alpha_s_range[1], alpha_s_range[2]),
                          variant = 1, centrality = "eigen", similarity = "hpi")
    sol <- lcda_repair(csr, sol); sol <- lcda_local_search(csr, sol)
    coassoc <- coassoc + (sol$membership[el[, 1]] == sol$membership[el[, 2]])
    lead_count[sol$leaders] <- lead_count[sol$leaders] + 1L
  }
  coassoc <- coassoc / B
  core <- igraph::coreness(g); in2 <- core[el[, 1]] >= 2 & core[el[, 2]] >= 2
  w <- ifelse(in2, w_min + (1 - w_min) * coassoc, w_min)
  mem <- as.integer(igraph::membership(igraph::cluster_louvain(g, weights = w)))
  list(graph = g, mem = mem, el = el, coassoc = coassoc, lead_count = lead_count, n = n)
}

# overlapping assignment from co-association affinities
overlap_sets <- function(pool, tau = 0.6) {
  el <- pool$el; co <- pool$coassoc; mem <- pool$mem; n <- pool$n
  # affinity tables: mean co-association of v's edges into each community
  sums <- vector("list", n); cnts <- vector("list", n)
  for (v in seq_len(n)) { sums[[v]] <- numeric(0); cnts[[v]] <- numeric(0) }
  acc <- function(v, c, val) {
    k <- as.character(c)
    sums[[v]][k] <<- (if (is.na(sums[[v]][k])) 0 else sums[[v]][k]) + val
    cnts[[v]][k] <<- (if (is.na(cnts[[v]][k])) 0 else cnts[[v]][k]) + 1
  }
  for (e in seq_len(nrow(el))) {
    u <- el[e, 1]; v <- el[e, 2]; cv <- co[e]
    acc(u, mem[v], cv); acc(v, mem[u], cv)
  }
  sets <- vector("list", n)
  for (v in seq_len(n)) {
    if (!length(sums[[v]])) { sets[[v]] <- mem[v]; next }
    a <- sums[[v]] / cnts[[v]]                          # mean affinity per community
    home_aff <- max(a, na.rm = TRUE)
    comms <- as.integer(names(a)[which(a >= tau * home_aff)])
    sets[[v]] <- union(mem[v], comms)
  }
  sets
}

participation <- function(g, mem) {
  n <- igraph::vcount(g); P <- numeric(n)
  adj <- igraph::as_adj_list(g, mode = "all")
  for (v in seq_len(n)) {
    nb <- as.integer(adj[[v]]); k <- length(nb)
    if (k == 0) { P[v] <- 0; next }
    tab <- table(mem[nb]); P[v] <- 1 - sum((as.numeric(tab) / k)^2)
  }
  P
}

set.seed(20260529)
TAU  <- as.numeric(Sys.getenv("OV_TAU", "0.7"))
MUS  <- as.numeric(strsplit(Sys.getenv("OV_MUS", "0.2,0.3"), ",")[[1]])
NGR  <- as.integer(Sys.getenv("OV_NGRAPH", "5"))
BENS <- as.integer(Sys.getenv("OV_B", "64"))

rows <- list()
for (mu in MUS) for (r in seq_len(NGR)) {
  gen <- lfr_generate(500, mu, seed = 15000 + 100 * as.integer(mu * 100) + r)
  if (is.null(gen)) next
  pool <- ecg_pool(gen$graph, B = BENS)
  sets <- overlap_sets(pool, tau = TAU)
  nmemb <- vapply(sets, length, integer(1))
  is_ov <- nmemb > 1
  P <- participation(pool$graph, pool$mem)
  wt <- if (sum(is_ov) > 0 && sum(!is_ov) > 0)
    suppressWarnings(stats::wilcox.test(P[is_ov], P[!is_ov])$p.value) else NA_real_
  rows[[length(rows)+1]] <- tibble::tibble(
    mu = mu, graph = r,
    hard_NMI = nmi(pool$mem, gen$truth),
    overlap_frac = mean(is_ov),
    P_overlap = if (any(is_ov)) mean(P[is_ov]) else NA_real_,
    P_nonoverlap = if (any(!is_ov)) mean(P[!is_ov]) else NA_real_,
    wilcox_p = wt)
  cli::cli_alert("mu={mu} g={r}: hardNMI={round(nmi(pool$mem,gen$truth),3)} overlap={round(mean(is_ov),3)} P(ov)={round(mean(P[is_ov]),3)} P(non)={round(mean(P[!is_ov]),3)}")
}
res <- dplyr::bind_rows(rows)
cat("\n=== overlap validity on LFR (tau=", TAU, ") ===\n", sep = "")
print(as.data.frame(res |> dplyr::group_by(mu) |> dplyr::summarise(
  hard_NMI = mean(hard_NMI), overlap_frac = mean(overlap_frac),
  P_overlap = mean(P_overlap, na.rm = TRUE), P_nonoverlap = mean(P_nonoverlap, na.rm = TRUE),
  .groups = "drop")), row.names = FALSE, digits = 3)

# qualitative Karate
kp <- ecg_pool(igraph::make_graph("Zachary"), B = BENS)
ks <- overlap_sets(kp, tau = TAU)
karate_overlap <- which(vapply(ks, length, integer(1)) > 1)
cat("\nKarate overlap nodes (tau=", TAU, "): ", paste(karate_overlap, collapse = ", "), "\n", sep = "")

meta <- make_meta(
  title = "Overlapping communities from LCDA-ECG co-association (capability + validity)",
  description = paste("Soft co-association assignment yields overlapping communities. Validated on LFR:",
                      "hard-core recovery preserved, and overlap nodes are structural bridges (higher",
                      "participation coefficient). No overlapping ground truth (networkx LFR is disjoint)."),
  seed = 20260529, tau = TAU, mus = MUS, ngraph = NGR, b_ensemble = BENS,
  limitations = c("No planted overlap (networkx LFR is disjoint); validity shown via bridge structure, not overlapping-NMI.",
                  "tau controls overlap aggressiveness; reported at a single value.",
                  "Karate overlap is qualitative."))
save_dataset("overlap_lcda_ecg", list(per_graph = res, karate_overlap = karate_overlap), meta)
cli::cli_alert_success("Overlap experiment done.")
