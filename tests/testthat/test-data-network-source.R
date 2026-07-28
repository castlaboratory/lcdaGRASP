# Provenance guard: meta$network_source must name the network a dataset was
# actually computed on.
#
# History: until 0.3.2, data-raw/00_helpers.R::make_meta() hard-coded
# network_source to the Newman/Zachary benchmark string and swallowed every
# other argument into `extra`. Every shipped dataset therefore claimed the five
# classical benchmarks as its origin -- including the synthetic LFR sweeps, the
# Shchur et al. (2018) real networks and the OpenAlex extraction -- and the
# three generators that DID pass a correct source had it silently redirected to
# `extra`. These tests assert the invariants that make that failure impossible
# to reintroduce, without pinning any exact wording.

extdata_dir <- function() system.file("extdata", package = "lcdaGRASP")

all_meta <- function() {
  fs <- list.files(extdata_dir(), pattern = "\\.rds$", full.names = TRUE)
  stats::setNames(lapply(fs, function(f) readRDS(f)$meta),
                  sub("\\.rds$", "", basename(fs)))
}

skip_if_no_extdata <- function() {
  skip_if(!nzchar(extdata_dir()) ||
            length(list.files(extdata_dir(), pattern = "\\.rds$")) == 0,
          "no extdata datasets installed")
}

# Datasets whose generator never touches an empirical benchmark network:
# purely synthetic (LFR via networkx, or igraph::sample_sbm), or an external
# corpus of its own. Derived from the data-raw/ script that writes each file.
SYNTHETIC_OR_EXTERNAL <- c(
  # synthetic LFR
  "lfr_robustness", "doe_lfr_recovery", "largescale", "gnn_baseline",
  "leader_vs_twostage_lfr", "lcda_ecg", "consensus_leader_test",
  "stability_pool", "weighted_demo",
  # purely synthetic SBM
  "pool_sensitivity",
  # external real networks / APIs
  "realnet_coauthor", "realnet_amazon", "openalex_leaders"
)

# Datasets computed on one or more of the five classical benchmarks.
BENCHMARK_BASED <- c(
  "repro_benchmarks", "repro_summary", "blogs_timing", "blogs_table9",
  "doe_screening", "doe_rsm", "eda_replicates", "leader_utility",
  "leader_vs_twostage", "degeneracy", "nce_alternatives",
  "prop6_summary", "prop6_trajectories", "overlap_lcda_ecg", "vnmi_nprime"
)

test_that("every shipped dataset declares a non-empty network_source", {
  skip_if_no_extdata()
  metas <- all_meta()
  bad <- names(metas)[!vapply(metas, function(m) {
    is.character(m$network_source) && length(m$network_source) == 1L &&
      nzchar(trimws(m$network_source))
  }, logical(1))]
  expect_true(length(bad) == 0,
              info = paste0("missing/empty meta$network_source: ",
                            paste(bad, collapse = ", ")))
})

test_that("network_source is not one constant pasted over every dataset", {
  skip_if_no_extdata()
  metas <- all_meta()
  skip_if(length(metas) < 3, "too few datasets to test for variation")
  srcs <- vapply(metas, function(m) as.character(m$network_source)[1], character(1))
  # The collection spans several distinct network families (classical
  # benchmarks, synthetic LFR, synthetic SBM, external corpora), so a single
  # value covering everything is the exact bug this guards against.
  expect_gt(length(unique(srcs)), 1L)
  expect_lt(max(table(srcs)) / length(srcs), 1)
})

test_that("synthetic and external datasets do not claim the benchmark collection", {
  skip_if_no_extdata()
  metas <- all_meta()
  present <- intersect(SYNTHETIC_OR_EXTERNAL, names(metas))
  skip_if(length(present) == 0, "no synthetic/external datasets installed")
  # These were never computed on Newman's GML networks or on the karate club,
  # so citing either is a false provenance claim regardless of phrasing.
  offending <- present[vapply(present, function(nm) {
    grepl("newman|umich\\.edu|zachary|karate", metas[[nm]]$network_source,
          ignore.case = TRUE)
  }, logical(1))]
  expect_true(length(offending) == 0,
              info = paste0("synthetic/external datasets citing the classical ",
                            "benchmarks: ", paste(offending, collapse = ", ")))
})

test_that("benchmark-based datasets do name a benchmark network", {
  skip_if_no_extdata()
  metas <- all_meta()
  present <- intersect(BENCHMARK_BASED, names(metas))
  skip_if(length(present) == 0, "no benchmark-based datasets installed")
  silent <- present[!vapply(present, function(nm) {
    grepl("karate|zachary|dolphins|football|polbooks|polblogs",
          metas[[nm]]$network_source, ignore.case = TRUE)
  }, logical(1))]
  expect_true(length(silent) == 0,
              info = paste0("benchmark-based datasets naming no benchmark: ",
                            paste(silent, collapse = ", ")))
})

test_that("no dataset hides a second network_source inside meta$extra", {
  skip_if_no_extdata()
  metas <- all_meta()
  # make_meta() takes network_source as a formal argument; a copy surviving in
  # `extra` means it was passed through `...` again and the real field may be
  # stale (this is precisely how the 0.3.1 datasets went wrong).
  shadowed <- names(metas)[vapply(metas, function(m) {
    !is.null(m$extra) && !is.null(m$extra$network_source)
  }, logical(1))]
  expect_true(length(shadowed) == 0,
              info = paste0("network_source shadowed in meta$extra: ",
                            paste(shadowed, collapse = ", ")))
})

test_that("every installed dataset is classified by this guard", {
  skip_if_no_extdata()
  metas <- all_meta()
  # If a new dataset ships without being placed in one of the two families
  # above, the checks silently stop covering it.
  unclassified <- setdiff(names(metas),
                          c(SYNTHETIC_OR_EXTERNAL, BENCHMARK_BASED))
  expect_true(length(unclassified) == 0,
              info = paste0("datasets not covered by the provenance guard: ",
                            paste(unclassified, collapse = ", ")))
})
