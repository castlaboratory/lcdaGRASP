# data-raw/run_all_data.R
#
# Regenerate every precomputed dataset under inst/extdata/. Run ONCE (it is the
# expensive step); the vignettes then render from the cached .rds without ever
# re-running a simulation. Safe to re-run: each step is wrapped so a failure in
# one does not abort the others.
#
# Usage (from repo root OR package dir):
#   Rscript lcdaGRASP/data-raw/run_all_data.R
# Override replication for a quick smoke run:
#   REPRO_NREP=3 Rscript lcdaGRASP/data-raw/run_all_data.R

here <- Filter(file.exists, c("data-raw", "lcdaGRASP/data-raw"))[1]
source(file.path(here, "00_helpers.R"))

# Each step maps to a sentinel output; if it already exists the step is skipped
# (set FORCE=1 to regenerate everything).
steps <- c("10_repro_benchmarks.R", "20_doe_experiment.R", "40_prop6.R",
           "80_nce.R", "60_degeneracy.R", "70_eda.R", "50_pool_sensitivity.R",
           "30_lfr_robustness.R")
sentinels <- c("10_repro_benchmarks.R" = "repro_summary",
               "20_doe_experiment.R"   = "doe_rsm",
               "40_prop6.R"            = "prop6_summary",
               "80_nce.R"              = "nce_alternatives",
               "60_degeneracy.R"       = "degeneracy",
               "70_eda.R"              = "eda_replicates",
               "50_pool_sensitivity.R" = "pool_sensitivity",
               "30_lfr_robustness.R"   = "lfr_robustness")
force <- nzchar(Sys.getenv("FORCE"))

overall <- Sys.time()
for (s in steps) {
  sentinel <- file.path(EXTDATA_DIR, paste0(sentinels[[s]], ".rds"))
  if (!force && file.exists(sentinel)) {
    cli::cli_alert_info("SKIP {s} ({sentinels[[s]]}.rds already present)")
    next
  }
  cli::cli_h1("RUN {s}")
  t0 <- Sys.time()
  ok <- tryCatch({ source(file.path(here, s), local = new.env()); TRUE },
                 error = function(e) { cli::cli_alert_danger("{s} FAILED: {conditionMessage(e)}"); FALSE })
  dt <- round(as.numeric(Sys.time() - t0, units = "mins"), 1)
  if (ok) cli::cli_alert_success("{s} done in {dt} min")
}
cli::cli_alert_success("All data steps finished in {round(as.numeric(Sys.time()-overall, units='mins'),1)} min")
cli::cli_alert_info("Datasets in: {.path {EXTDATA_DIR}}")
print(list.files(EXTDATA_DIR))
