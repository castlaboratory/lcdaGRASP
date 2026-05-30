# data-raw/20_doe_experiment.R
#
# Design of Experiments (DoE) a la George Box for tuning LCDA-GRASP.
# Sequential strategy:
#   Phase 1  SCREENING - full factorial over the categorical/discrete factors
#            (variant, centrality, similarity) crossed with a 2-level coding of
#            (alpha_c, alpha_s). Response = best modularity Q. Replicated for
#            pure error. Run on several benchmark networks (coverage - Wald).
#   Phase 2  RESPONSE SURFACE - a rotatable Central Composite Design (CCD) in
#            (alpha_c, alpha_s) at the winning categorical settings, to fit a
#            second-order model and locate the optimum (canonical analysis is
#            done in the vignette with lm/rsm).
#
# -> inst/extdata/doe_screening.rds
# -> inst/extdata/doe_rsm.rds
# Only generates the runs; the ANOVA/Pareto/RSM analysis lives in the
# doe-parameter-tuning vignette.

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages(library(purrr))

set.seed(424242)

run_grasp <- function(g, alpha_c, alpha_s, variant, centrality, similarity, B, seed) {
  t0 <- Sys.time()
  r <- lcda_grasp(g, alpha_c = alpha_c, alpha_s = alpha_s, variant = variant,
                  B = B, centrality = centrality, similarity = similarity, seed = seed)
  tibble::tibble(Q = r$best$Q, H = r$best$H, g_comms = length(r$best$leaders),
                 secs = as.numeric(Sys.time() - t0, units = "secs"))
}

# ---- Phase 1: screening ----------------------------------------------------
SCREEN_NETS <- c("karate", "dolphins", "football", "polbooks")
SCREEN_REPS <- 10L
B_SCREEN    <- 50L
# 2-level coding of the continuous RCL parameters.
ac_lo <- 0.10; ac_hi <- 0.50
as_lo <- 0.20; as_hi <- 0.40

screen_grid <- tidyr::expand_grid(
  network    = SCREEN_NETS,
  ac_code    = c(-1, 1),
  as_code    = c(-1, 1),
  variant    = c(1L, 2L),
  centrality = c("eigen", "closeness", "betweenness"),
  similarity = c("hpi", "dice", "jaccard"),
  rep        = seq_len(SCREEN_REPS)
) |>
  dplyr::mutate(alpha_c = ifelse(ac_code < 0, ac_lo, ac_hi),
                alpha_s = ifelse(as_code < 0, as_lo, as_hi))

cli::cli_h1("DoE Phase 1: screening ({nrow(screen_grid)} runs)")
nets <- purrr::set_names(lapply(SCREEN_NETS, load_benchmark), SCREEN_NETS)

screening <- screen_grid |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  purrr::pmap_dfr(function(network, ac_code, as_code, variant, centrality,
                          similarity, rep, alpha_c, alpha_s, .row) {
    if (.row %% 200 == 0) cli::cli_alert_info("screening {(.row)}/{nrow(screen_grid)}")
    out <- run_grasp(nets[[network]], alpha_c, alpha_s, variant, centrality,
                     similarity, B_SCREEN, seed = 7000 + .row)
    dplyr::bind_cols(
      tibble::tibble(network, ac_code, as_code, alpha_c, alpha_s, variant,
                     centrality, similarity, rep), out)
  })

meta_screen <- make_meta(
  title = "DoE screening: full factorial (variant x centrality x similarity x 2-level alpha)",
  description = "LCDA-GRASP best Q over a 2x2x2x3x3 factorial, replicated, on 4 benchmark networks. Drives ANOVA / effect (Pareto) analysis to identify which factors matter.",
  seed = 424242, reps = SCREEN_REPS, B = B_SCREEN,
  coding = list(alpha_c = c(`-1` = ac_lo, `+1` = ac_hi),
                alpha_s = c(`-1` = as_lo, `+1` = as_hi)),
  limitations = c(
    "PolBlogs excluded from screening for cost; confirmed separately in RSM phase.",
    "Q is the best over B=50 GRASP iterations, not a single construction.",
    "Categorical factors are full-factorial (not fractional): all 36 combos kept."))
save_dataset("doe_screening", screening, meta_screen)

# ---- Phase 2: response surface (rotatable CCD in alpha_c x alpha_s) --------
# Coded design: 4 factorial (+/-1), 4 axial (+/-alpha, 0)/(0,+/-alpha),
# 5 center points. alpha = sqrt(2) for rotatability.
ccd_alpha <- sqrt(2)
ccd_coded <- rbind(
  expand.grid(ac = c(-1, 1), as = c(-1, 1)),      # factorial
  data.frame(ac = c(-ccd_alpha, ccd_alpha, 0, 0), # axial
             as = c(0, 0, -ccd_alpha, ccd_alpha)),
  data.frame(ac = rep(0, 5), as = rep(0, 5))      # center (pure error)
)
# Map coded -> actual. Centers chosen around the interesting region; half-ranges
# keep the axial points inside (0, 1).
ac_c0 <- 0.30; ac_hr <- 0.18
as_c0 <- 0.30; as_hr <- 0.13
to_actual <- function(code, c0, hr) pmin(0.95, pmax(0.02, c0 + code * hr))

RSM_REPS <- 12L
B_RSM    <- 50L
# Winning categorical settings are confirmed by the screening analysis; we use
# the paper's construction (variant 1, eigen, hpi) as the RSM base and also
# record variant 2 for comparison.
rsm_grid <- tidyr::expand_grid(
  network    = BENCHMARKS$network,
  point      = seq_len(nrow(ccd_coded)),
  variant    = c(1L, 2L),
  rep        = seq_len(RSM_REPS)
) |>
  dplyr::mutate(ac_code = ccd_coded$ac[point], as_code = ccd_coded$as[point],
                alpha_c = to_actual(ac_code, ac_c0, ac_hr),
                alpha_s = to_actual(as_code, as_c0, as_hr))

cli::cli_h1("DoE Phase 2: RSM/CCD ({nrow(rsm_grid)} runs)")
all_nets <- purrr::set_names(lapply(BENCHMARKS$network, load_benchmark), BENCHMARKS$network)

rsm <- rsm_grid |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  purrr::pmap_dfr(function(network, point, variant, rep, ac_code, as_code,
                           alpha_c, alpha_s, .row) {
    if (.row %% 100 == 0) cli::cli_alert_info("rsm {(.row)}/{nrow(rsm_grid)}")
    out <- run_grasp(all_nets[[network]], alpha_c, alpha_s, variant,
                     "eigen", "hpi", B_RSM, seed = 90000 + .row)
    dplyr::bind_cols(
      tibble::tibble(network, point, variant, rep, ac_code, as_code,
                     alpha_c, alpha_s), out)
  })

meta_rsm <- make_meta(
  title = "DoE response surface: rotatable CCD in (alpha_c, alpha_s)",
  description = "LCDA-GRASP best Q over a 13-point central composite design (4 factorial + 4 axial + 5 center) in coded (alpha_c, alpha_s), replicated, on all 5 networks, for variants 1 and 2. Drives the second-order model + canonical analysis in the vignette.",
  seed = 424242, reps = RSM_REPS, B = B_RSM, ccd_alpha = ccd_alpha,
  coding = list(alpha_c = c(center = ac_c0, half_range = ac_hr),
                alpha_s = c(center = as_c0, half_range = as_hr)),
  limitations = c(
    "Categorical factors fixed at the paper construction (eigen/hpi); variant kept as a 2-level factor.",
    "Axial points clamped to (0.02, 0.95) to keep alpha in a valid RCL range.",
    "On small networks the Q surface is near-flat at the top (links to modularity degeneracy)."))
save_dataset("doe_rsm", rsm, meta_rsm)
