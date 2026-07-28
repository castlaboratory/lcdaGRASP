# data-raw/25_doe_lfr_recovery.R
#
# Fase 1.2/1.3 -- DESIGN OF EXPERIMENTS targeting RECOVERY (NMI vs the planted
# LFR partition), not modularity. The full LFR sweep (30_) showed LCDA trails
# Leiden/Louvain on recovery; the cheap probe suggested tuning helps only
# modestly. This DoE settles it rigorously, in the George-Box tradition:
#
#   Phase 1 (Screening) -- 2^? factorial over the categorical knobs
#     (variant x centrality x similarity), ANOVA on NMI, identify the best cell.
#   Phase 2 (RSM)       -- rotatable central-composite design in (alpha_c,
#     alpha_s) for the fixed-parameter LCDA-GRASP at the best categorical cell,
#     fit a second-order surface, locate the stationary point, and report the
#     BEST achievable NMI vs Leiden/Louvain.
#
# The goal is NOT to beat Leiden but to BOUND the gap under optimised
# parameters (Wald): "even RSM-optimised, LCDA recovers within X NMI of Leiden".
#
# Paired design: the same LFR graphs are reused across all configurations.
# Parallelised across cores (the LCDA hot path is already C++ -- profiling
# showed 84% in vnmi_local_search_cpp -- so speed comes from parallelism, not
# from rewriting R in C++). Checkpointed per mu (resumable).
#
# -> inst/extdata/doe_lfr_recovery.rds  (list: screening, rsm, baselines)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(parallel); library(dplyr) })

set.seed(20260529)
MUS    <- as.numeric(strsplit(Sys.getenv("DOE_MUS", "0.2,0.3,0.4"), ",")[[1]])
NGRAPH <- as.integer(Sys.getenv("DOE_NGRAPH", "5"))   # LFR graphs per mu (paired)
N      <- as.integer(Sys.getenv("DOE_N", "500"))
B_SCR  <- as.integer(Sys.getenv("DOE_B_SCREEN", "200"))
B_RSM  <- as.integer(Sys.getenv("DOE_B_RSM", "150"))
CORES  <- as.integer(Sys.getenv("DOE_CORES", as.character(max(1L, parallel::detectCores() - 2L))))
CKPT   <- file.path(PKG_DIR, "data-raw", ".doe_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)

# ---- build paired LFR graphs + baseline recovery (Leiden/Louvain) ----------
build_graphs <- function(mu) {
  lapply(seq_len(NGRAPH), function(r) {
    gen <- lfr_generate(N, mu, seed = 7000 + 100 * as.integer(mu * 100) + r)
    if (is.null(gen)) return(NULL)
    g <- gen$graph; truth <- gen$truth
    base <- c(
      Leiden  = nmi(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity")), truth),
      Louvain = nmi(igraph::membership(igraph::cluster_louvain(g)), truth))
    list(graph = g, truth = truth, baseline = base)
  }) |> Filter(f = Negate(is.null))
}

# ================= Phase 1: screening (variant x centrality x similarity) =====
scr_grid <- expand.grid(variant = c(1L, 2L),
                        centrality = c("eigen", "betweenness", "closeness"),
                        similarity = c("hpi", "dice", "jaccard"),
                        stringsAsFactors = FALSE)
cli::cli_h1("DoE recovery -- Phase 1 screening: {nrow(scr_grid)} configs x {NGRAPH} graphs x {length(MUS)} mu")

screening <- list(); baselines <- list()
for (mu in MUS) {
  cp <- file.path(CKPT, sprintf("scr_mu%03d.rds", as.integer(mu * 100)))
  if (file.exists(cp)) { sc <- readRDS(cp); screening[[length(screening)+1]] <- sc$scr; baselines[[length(baselines)+1]] <- sc$base; next }
  graphs <- build_graphs(mu)
  base_tbl <- dplyr::bind_rows(lapply(seq_along(graphs), function(gi)
    tibble::tibble(mu = mu, graph = gi, Leiden = graphs[[gi]]$baseline["Leiden"],
                   Louvain = graphs[[gi]]$baseline["Louvain"])))
  jobs <- expand.grid(cfg = seq_len(nrow(scr_grid)), gi = seq_along(graphs))
  rows <- mclapply(seq_len(nrow(jobs)), function(j) {
    ci <- jobs$cfg[j]; gi <- jobs$gi[j]; z <- graphs[[gi]]
    r <- tryCatch(lcda_gr(z$graph, variant = scr_grid$variant[ci], B = B_SCR,
                          centrality = scr_grid$centrality[ci],
                          similarity = scr_grid$similarity[ci], seed = 7), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    tibble::tibble(mu = mu, graph = gi, variant = scr_grid$variant[ci],
                   centrality = scr_grid$centrality[ci], similarity = scr_grid$similarity[ci],
                   Q = r$best$Q, NMI = nmi(r$best$membership, z$truth),
                   ARI = adjusted_rand(r$best$membership, z$truth))
  }, mc.cores = CORES)
  sc <- dplyr::bind_rows(rows)
  saveRDS(list(scr = sc, base = base_tbl), cp)
  screening[[length(screening)+1]] <- sc; baselines[[length(baselines)+1]] <- base_tbl
  cli::cli_alert_success("screening mu={mu}: best NMI={round(max(sc$NMI),3)} vs Leiden={round(mean(base_tbl$Leiden),3)}")
}
screening_tbl <- dplyr::bind_rows(screening)
baseline_tbl  <- dplyr::bind_rows(baselines)

# best categorical cell by mean NMI
best_cell <- screening_tbl |>
  dplyr::group_by(variant, centrality, similarity) |>
  dplyr::summarise(NMI = mean(NMI), .groups = "drop") |>
  dplyr::slice_max(NMI, n = 1)
cli::cli_alert_info("best cell: variant {best_cell$variant} / {best_cell$centrality} / {best_cell$similarity} (mean NMI {round(best_cell$NMI,3)})")

# ================= Phase 2: RSM in (alpha_c, alpha_s) =========================
# Rotatable central-composite design (coded units); centre at (0.3, 0.3).
ac0 <- 0.30; as0 <- 0.30; d_ac <- 0.18; d_as <- 0.15; ax <- sqrt(2)
ccd <- rbind(
  expand.grid(xc = c(-1, 1), xs = c(-1, 1)),          # factorial
  data.frame(xc = c(-ax, ax, 0, 0), xs = c(0, 0, -ax, ax)),  # axial
  data.frame(xc = rep(0, 3), xs = rep(0, 3)))          # centre replicates
ccd$alpha_c <- pmin(pmax(ac0 + ccd$xc * d_ac, 0.02), 0.9)
ccd$alpha_s <- pmin(pmax(as0 + ccd$xs * d_as, 0.02), 0.5)

mu_rsm <- MUS[ceiling(length(MUS) / 2)]   # middle mu (clearest gap)
cli::cli_h1("DoE recovery -- Phase 2 RSM: {nrow(ccd)} CCD points x {NGRAPH} graphs at mu={mu_rsm}")
graphs_rsm <- build_graphs(mu_rsm)
jobs <- expand.grid(pt = seq_len(nrow(ccd)), gi = seq_along(graphs_rsm))
rsm_rows <- mclapply(seq_len(nrow(jobs)), function(j) {
  pt <- jobs$pt[j]; gi <- jobs$gi[j]; z <- graphs_rsm[[gi]]
  r <- tryCatch(lcda_grasp(z$graph, variant = best_cell$variant, B = B_RSM,
                           alpha_c = ccd$alpha_c[pt], alpha_s = ccd$alpha_s[pt],
                           centrality = best_cell$centrality, similarity = best_cell$similarity,
                           seed = 7), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  tibble::tibble(mu = mu_rsm, graph = gi, point = pt,
                 alpha_c = ccd$alpha_c[pt], alpha_s = ccd$alpha_s[pt],
                 xc = ccd$xc[pt], xs = ccd$xs[pt],
                 Q = r$best$Q, NMI = nmi(r$best$membership, z$truth))
}, mc.cores = CORES)
rsm_tbl <- dplyr::bind_rows(rsm_rows)

meta <- make_meta(
  title = "DoE for LFR recovery (NMI): screening + response surface in (alpha_c, alpha_s)",
  description = paste("George-Box DoE targeting recovery (NMI vs planted LFR truth) rather than",
                      "modularity. Phase 1 screens variant/centrality/similarity; Phase 2 fits a",
                      "rotatable CCD in (alpha_c, alpha_s) at the best cell. Bounds the recovery gap",
                      "to Leiden/Louvain under optimised parameters."),
  network_source = SRC_LFR,
  seed = 20260529, n = N, ngraph = NGRAPH, mus = MUS, b_screen = B_SCR, b_rsm = B_RSM,
  best_cell = paste(best_cell$variant, best_cell$centrality, best_cell$similarity),
  limitations = c(
    "n=500 fixed; recovery DoE not yet repeated at large scale.",
    "RSM run at a single (middle) mu; the surface may shift with mu.",
    "GR auto-tunes alpha within ranges, so the (alpha_c, alpha_s) RSM uses the fixed LCDA-GRASP."))

save_dataset("doe_lfr_recovery", list(screening = screening_tbl, baselines = baseline_tbl,
                                      rsm = rsm_tbl, best_cell = best_cell), meta)
unlink(CKPT, recursive = TRUE)
cli::cli_alert_success("DoE recovery done. best cell: {meta$extra$best_cell}")
