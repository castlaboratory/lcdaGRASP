# data-raw/30_lfr_robustness.R
#
# Fase 1.2 (Wald + leader+community thesis): sweep the LFR mixing parameter mu
# on the CANONICAL Lancichinetti-Fortunato-Radicchi (2008) benchmark (generated
# via networkx in .venv-lfr; see lfr_generate.py), replacing the previous
# degree-corrected SBM proxy.
#
# Two questions, one sweep:
#   (A) Recovery robustness: Q, ARI and NMI vs the planted partition for
#       LCDA-GRASP/GR and the igraph baselines, across community-strength
#       regimes mu in [0.1, 0.8]. ARI/NMI reveal whether the recovered
#       partition is CORRECT, which Q alone cannot show.
#   (B) Leader-utility hypothesis (anchors Fase 5): the smoke study (1.9)
#       suggested the LCDA leaders win as spreaders precisely where community
#       structure is strong (Football). Here we test that systematically: does
#       the IC-spread advantage of the LCDA leaders over a top-degree seed set
#       of the same size GROW as mu decreases (stronger communities)?
#
# -> inst/extdata/lfr_robustness.rds  (list: quality, leader_utility)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(purrr); library(dplyr) })

set.seed(20260529)
MUS  <- as.numeric(strsplit(Sys.getenv("LFR_MUS", "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8"), ",")[[1]])
REPS <- as.integer(Sys.getenv("LFR_REPS", "10"))   # 30 for the publication run
N    <- as.integer(Sys.getenv("LFR_N", "500"))
MC   <- as.integer(Sys.getenv("LFR_MC", "500"))    # IC Monte Carlo per seed set
PIC  <- as.numeric(Sys.getenv("LFR_PIC", "0.05"))  # IC propagation probability
B_GRASP <- 50L; B_GR <- 150L

methods <- c("LCDA-GRASP 1", "LCDA-GR 1", "Louvain", "Leiden", "FastGreedy", "Walktrap")
grid <- tidyr::expand_grid(mu = MUS, rep = seq_len(REPS))
cli::cli_h1("Canonical LFR sweep: {nrow(grid)} (mu, rep) cells, n={N}")

# ---- per-cell checkpointing (resumable across restarts) --------------------
# Each completed (mu, rep) cell is written to .lfr_ckpt/; a restart skips cells
# already present, so an interrupted run loses at most one cell. The directory
# is removed only after the final dataset is saved successfully.
CKPT_DIR <- file.path(PKG_DIR, "data-raw", ".lfr_ckpt")
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
ckpt_path <- function(mu, rep) file.path(CKPT_DIR, sprintf("cell_mu%03d_rep%03d.rds", as.integer(mu * 100), rep))

quality <- list(); utility <- list()
done <- 0L
for (i in seq_len(nrow(grid))) {
  mu <- grid$mu[i]; rep <- grid$rep[i]
  cp <- ckpt_path(mu, rep)
  if (file.exists(cp)) {                         # resume: reuse cached cell
    cell <- readRDS(cp)
    if (!is.null(cell$quality)) quality[[length(quality) + 1]] <- cell$quality
    if (!is.null(cell$utility)) utility[[length(utility) + 1]] <- cell$utility
    done <- done + 1L
    next
  }
  gen <- lfr_generate(N, mu, seed = 1000 * rep + as.integer(mu * 100))
  if (is.null(gen)) { cli::cli_alert_warning("LFR gen failed (mu={mu}, rep={rep})"); next }
  g <- gen$graph; truth <- gen$truth
  if (i %% 10 == 1) cli::cli_alert_info("cell {i}/{nrow(grid)} (mu={mu}, rep={rep}, n={igraph::vcount(g)}, m={igraph::ecount(g)})")
  cell_q <- list(); cell_u <- NULL

  # ---- (A) recovery quality ----
  lcda_gr_fit <- NULL
  for (mtd in methods) {
    t0 <- Sys.time()
    out <- tryCatch(switch(mtd,
      "LCDA-GRASP 1" = { r <- lcda_grasp(g, B = B_GRASP, seed = rep); list(mem = r$best$membership, Q = r$best$Q) },
      "LCDA-GR 1"    = { r <- lcda_gr(g, variant = 1, B = B_GR, seed = rep); lcda_gr_fit <<- r; list(mem = r$best$membership, Q = r$best$Q) },
      "Louvain"      = { r <- igraph::cluster_louvain(g);     list(mem = igraph::membership(r), Q = igraph::modularity(r)) },
      "Leiden"       = { r <- igraph::cluster_leiden(g, objective_function = "modularity"); list(mem = igraph::membership(r), Q = igraph::modularity(g, igraph::membership(r))) },
      "FastGreedy"   = { r <- igraph::cluster_fast_greedy(igraph::simplify(g)); list(mem = igraph::membership(r), Q = igraph::modularity(r)) },
      "Walktrap"     = { r <- igraph::cluster_walktrap(g);    list(mem = igraph::membership(r), Q = igraph::modularity(r)) }
    ), error = function(e) NULL)
    if (is.null(out)) next
    cell_q[[length(cell_q) + 1]] <- tibble::tibble(
      mu = mu, rep = rep, method = mtd, Q = out$Q,
      ARI = adjusted_rand(out$mem, truth), NMI = nmi(out$mem, truth),
      n_comms = length(unique(out$mem)), secs = as.numeric(Sys.time() - t0, units = "secs"))
  }

  # ---- (B) leader utility vs a top-degree seed set of equal size ----
  if (!is.null(lcda_gr_fit)) {
    nb <- nbr_list(g); deg <- igraph::degree(g)
    leaders <- lcda_gr_fit$best$leaders; gs <- length(leaders)
    seeds_deg <- order(deg, decreasing = TRUE)[seq_len(gs)]
    sp_lcda <- ic_spread(nb, leaders, PIC, MC)
    sp_deg  <- ic_spread(nb, seeds_deg, PIC, MC)
    sp_rand <- mean(vapply(seq_len(10), function(j) ic_spread(nb, sample.int(igraph::vcount(g), gs), PIC, max(50L, MC %/% 10L)), numeric(1)))
    cell_u <- tibble::tibble(
      mu = mu, rep = rep, g_size = gs,
      spread_lcda = sp_lcda, spread_degree = sp_deg, spread_random = sp_rand,
      adv_vs_degree = sp_lcda - sp_deg, adv_vs_random = sp_lcda - sp_rand)
  }

  # checkpoint this cell, then accumulate
  cell_q <- if (length(cell_q)) dplyr::bind_rows(cell_q) else NULL
  saveRDS(list(quality = cell_q, utility = cell_u), ckpt_path(mu, rep))
  if (!is.null(cell_q)) quality[[length(quality) + 1]] <- cell_q
  if (!is.null(cell_u)) utility[[length(utility) + 1]] <- cell_u
  done <- done + 1L
}
cli::cli_alert_info("{done}/{nrow(grid)} cells available (checkpointed in {.path {CKPT_DIR}})")

quality_tbl <- dplyr::bind_rows(quality)
utility_tbl <- dplyr::bind_rows(utility)

meta <- make_meta(
  title = "Canonical LFR sweep (mu in [0.1, 0.8]): recovery quality + leader-utility vs community strength",
  description = paste("Q/ARI/NMI vs the planted LFR partition for LCDA-GRASP/GR and igraph baselines,",
                      "plus the IC-spread advantage of the LCDA leaders over an equal-size top-degree",
                      "seed set as a function of mu. Tests whether the leader-utility advantage grows",
                      "as communities strengthen (mu decreases)."),
  seed = 20260529, reps = REPS, n = N, mc = MC, p_ic = PIC,
  generator = "canonical LFR (Lancichinetti-Fortunato-Radicchi 2008) via networkx 3.6.1 (.venv-lfr); tau1=2.5, tau2=1.5, avg_degree=12, comm in [20,60]",
  limitations = c(
    "n=500 fixed; size effects not explored here (see large-scale item 1.4).",
    "IC simulator is pure R (item 2.9 = Rcpp kernel); MC and reps kept modest for runtime.",
    "Leader-utility uses top-degree as the centrality comparator (strongest single baseline in 1.9)."))

save_dataset("lfr_robustness", list(quality = quality_tbl, leader_utility = utility_tbl), meta)
unlink(CKPT_DIR, recursive = TRUE)   # checkpoints no longer needed after a clean save
cli::cli_alert_success("Canonical LFR sweep done.")
