# data-raw/998_fix_network_source.R
#
# ONE-OFF metadata repair (0.3.2). Until this fix, make_meta() hard-coded
#
#   network_source = "Newman netdata collection (...) + igraph::make_graph('Zachary')"
#
# and collected everything else into `extra`, so EVERY shipped dataset claimed
# the five classical benchmarks as its origin -- including the synthetic LFR
# sweeps, the Shchur et al. (2018) real networks and the OpenAlex extraction.
# The three scripts that did supply a correct source passed it through `...`,
# where it was swallowed by `extra` and never reached `meta$network_source`.
#
# This script rewrites ONLY `meta$network_source` (and drops the now-redundant
# `meta$extra$network_source` copy) in the shipped .rds files. It touches no
# numbers: `$results` is compared byte-for-byte before and after and the run
# aborts if anything else moved. Correcting a FALSE provenance field in place is
# legitimate; restamping a TRUE one (e.g. `pkg_version`, cf. PR #49) is not.
#
# Idempotent. Run from the package root:
#   Rscript data-raw/998_fix_network_source.R
# then regenerate the manifests:
#   Rscript data-raw/999_manifest.R

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}

# The true source of every shipped dataset, transcribed from the generator that
# writes it (script name in the comment). Mirrors the `network_source =`
# argument now passed by each of those scripts.
SOURCES <- list(
  # --- classical benchmark networks -----------------------------------------
  repro_benchmarks   = src_newman(BENCHMARKS$network),                      # 10_
  repro_summary      = src_newman(BENCHMARKS$network),                      # 10_
  blogs_timing       = src_newman("polblogs"),                              # 11_
  blogs_table9       = paste(src_newman("polblogs"),                        # 12_
                             "Collapsed to an undirected simple graph whose edge weights are the directed multiplicities."),
  doe_screening      = src_newman(c("karate", "dolphins", "football", "polbooks")),  # 20_
  doe_rsm            = src_newman(BENCHMARKS$network),                      # 20_
  eda_replicates     = src_newman(c("karate", "dolphins", "football", "polbooks")),  # 70_
  leader_utility     = src_newman(c("karate", "dolphins", "football", "polbooks")),  # 90_
  leader_vs_twostage = src_newman(BENCHMARKS$network),                      # 95_

  # --- classical benchmarks mixed with synthetic SBM ------------------------
  degeneracy         = paste(src_newman(c("karate", "dolphins", "polbooks")),        # 60_
                             SRC_SBM, "-- SBM(800, 2 blocks, p_in=0.06, p_out=0.005)"),
  nce_alternatives   = paste(src_newman(c("karate", "dolphins", "polbooks")),        # 80_
                             SRC_SBM, "-- SBM(300, 5 blocks, p_in=0.12, p_out=0.02)"),
  prop6_summary      = paste(SRC_ZACHARY, "and", SRC_SBM,                   # 40_
                             "-- SBM(400, 4 blocks, p_in=0.15, p_out=0.02)"),
  prop6_trajectories = paste(SRC_ZACHARY, "and", SRC_SBM,                   # 40_
                             "-- SBM(400, 4 blocks, p_in=0.15, p_out=0.02)"),

  # --- purely synthetic SBM -------------------------------------------------
  pool_sensitivity   = paste(SRC_SBM, "-- SBM(300, 5 blocks, p_in=0.12, p_out=0.02)"),  # 50_

  # --- synthetic LFR (networkx, .venv-lfr) ----------------------------------
  lfr_robustness         = SRC_LFR,                                         # 30_ + 31_
  doe_lfr_recovery       = SRC_LFR,                                         # 25_
  largescale             = SRC_LFR,                                         # 35_
  gnn_baseline           = SRC_LFR,                                         # 45_
  leader_vs_twostage_lfr = SRC_LFR,                                         # 96_
  lcda_ecg               = SRC_LFR,                                         # 97_
  consensus_leader_test  = SRC_LFR,                                         # 98_
  stability_pool         = SRC_LFR,                                         # 100_
  weighted_demo      = paste(SRC_LFR,                                       # 55_
                             "-- edge weights are synthetic (U(1.5,4) intra-community, U(0.1,1) inter-community), not measured."),
  overlap_lcda_ecg   = paste(SRC_LFR, "-- plus a qualitative run on", SRC_ZACHARY),  # 99_
  vnmi_nprime        = paste0(SRC_LFR,                                      # 37_
                              " -- sparse (avg degree 12) and dense (avg degree 60) regimes.",
                              " Real-network regime: largest connected component of polblogs. ",
                              src_newman("polblogs")),

  # --- external real networks / APIs ----------------------------------------
  realnet_coauthor = "Coauthor-Physics, gnn-benchmark ms_academic_phy.npz (Shchur et al. 2018)",  # 38_
  realnet_amazon   = "Amazon-Computers, gnn-benchmark npz (Shchur, Mumme, Bojchevski & Gunnemann 2018)",  # 36_
  openalex_leaders = "OpenAlex (api.openalex.org) works in topic T10064; co-authorship edges, citation impact per author"  # 39_
)

files <- sort(list.files(EXTDATA_DIR, pattern = "[.]rds$"), method = "radix")
names_on_disk <- sub("[.]rds$", "", files)
missing <- setdiff(names_on_disk, names(SOURCES))
extra_n <- setdiff(names(SOURCES), names_on_disk)
if (length(missing)) stop("no source recorded for: ", paste(missing, collapse = ", "))
if (length(extra_n)) stop("source recorded for absent dataset: ", paste(extra_n, collapse = ", "))

changed <- 0L
for (f in files) {
  nm   <- sub("[.]rds$", "", f)
  path <- file.path(EXTDATA_DIR, f)
  obj  <- readRDS(path)
  before_results <- obj$results
  old <- obj$meta$network_source

  obj$meta$network_source <- SOURCES[[nm]]
  obj$meta$extra$network_source <- NULL   # de-duplicate the shadowed copy

  # Non-negotiable invariant: nothing but $meta may move.
  stopifnot(identical(before_results, obj$results))
  stopifnot(identical(names(obj), c("results", "meta")))

  if (identical(old, obj$meta$network_source)) {
    cli::cli_alert_info("{nm}: already correct")
    next
  }
  saveRDS(obj, path, compress = "xz")
  changed <- changed + 1L
  cli::cli_alert_success("{nm}: network_source corrected")
}

cli::cli_alert_info("{changed} of {length(files)} datasets rewritten. Now run data-raw/999_manifest.R.")
