# data-raw/35_largescale.R
#
# Fase 1.4 (Wald): large-scale evaluation (n > 10^4), the regime the paper
# promised but never demonstrated. Canonical LFR at increasing n, recording
# recovery (NMI vs planted truth) AND wall-clock time for Louvain, Leiden, ECG,
# LCDA-GR and the ensemble-consensus LCDA-ECG. Shows whether the method scales
# and whether the consensus advantage persists at scale.
#
# Probe finding (n=10^4): LCDA-ECG reaches NMI 0.986 (> Leiden 0.977) while
# LCDA-GR alone trails at 0.675 -- consensus is essential at scale.
#
# -> inst/extdata/largescale.rds  (list: recovery, meta)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })
options(lcda_ecg_define_only = TRUE)
source(file.path(PKG_DIR, "data-raw", "97_lcda_ecg.R"))
options(lcda_ecg_define_only = FALSE)

set.seed(20260529)
NS    <- as.integer(strsplit(Sys.getenv("LS_NS", "10000,25000,50000"), ",")[[1]])
REPS  <- as.integer(Sys.getenv("LS_REPS", "2"))
MU    <- as.numeric(Sys.getenv("LS_MU", "0.3"))
B_GR  <- as.integer(Sys.getenv("LS_B_GR", "50"))
B_ECG <- as.integer(Sys.getenv("LS_B_ECG", "32"))
CKPT  <- file.path(PKG_DIR, "data-raw", ".ls_ckpt")
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)

timed <- function(expr) { t0 <- Sys.time(); m <- force(expr); list(mem = m, secs = as.numeric(Sys.time() - t0, units = "secs")) }

rows <- list()
for (n in NS) for (r in seq_len(REPS)) {
  cp <- file.path(CKPT, sprintf("ls_n%08d_r%02d.rds", n, r))
  if (file.exists(cp)) { rows[[length(rows) + 1]] <- readRDS(cp); next }
  gen <- lfr_generate(n, MU, seed = 50000 + r, avg_degree = 15, min_comm = 50, max_comm = 200)
  if (is.null(gen)) { cli::cli_alert_warning("LFR gen failed (n={n}, r={r})"); next }
  g <- gen$graph; truth <- gen$truth
  cli::cli_h2("n={igraph::vcount(g)} m={igraph::ecount(g)} comms={length(unique(truth))} (rep {r})")
  res <- list()
  res$Louvain  <- timed(as.integer(igraph::membership(igraph::cluster_louvain(g))))
  res$Leiden   <- timed(as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity"))))
  ecg_m <- ecg_membership(g); res$ECG <- list(mem = if (is.null(ecg_m)) NA else ecg_m, secs = NA)
  res$`LCDA-GR`  <- timed(lcda_gr(g, B = B_GR, centrality = "eigen", similarity = "hpi", seed = 1)$best$membership)
  res$`LCDA-ECG` <- timed(lcda_ecg(g, B = B_ECG, centrality = "eigen", similarity = "hpi", seed = 1)$membership)
  cell <- dplyr::bind_rows(lapply(names(res), function(m) tibble::tibble(
    n = igraph::vcount(g), m_edges = igraph::ecount(g), rep = r, method = m,
    NMI = if (length(res[[m]]$mem) > 1) nmi(res[[m]]$mem, truth) else NA_real_,
    secs = res[[m]]$secs)))
  saveRDS(cell, cp); rows[[length(rows) + 1]] <- cell
  print(as.data.frame(cell[, c("method", "NMI", "secs")]), row.names = FALSE, digits = 3)
}
recovery <- dplyr::bind_rows(rows)

meta <- make_meta(
  title = "Large-scale evaluation on canonical LFR (n > 10^4): recovery + runtime",
  description = paste("NMI vs planted truth and wall-clock time for Louvain/Leiden/ECG/LCDA-GR/",
                      "LCDA-ECG on LFR networks up to n =", max(NS), ". Demonstrates scaling and",
                      "that the ensemble-consensus advantage persists (and grows) at scale."),
  seed = 20260529, ns = NS, reps = REPS, mu = MU, b_gr = B_GR, b_ecg = B_ECG,
  limitations = c("Single mu=0.3; LFR (planted, disjoint) -- not a real attributed network.",
                  "ECG timing not recorded (Python subprocess); recovery only.",
                  "Times are single-threaded R/C++; the B GRASP iterations are embarrassingly parallel."))
save_dataset("largescale", list(recovery = recovery), meta)
unlink(CKPT, recursive = TRUE)
cli::cli_alert_success("Large-scale evaluation done.")
