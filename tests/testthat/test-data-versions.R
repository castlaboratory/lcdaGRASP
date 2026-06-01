# Reproducibility guard: every shipped dataset must carry the current package
# version in its metadata, so the article and the package never mix versions.

test_that("all shipped inst/extdata datasets carry the current package version", {
  ver <- as.character(utils::packageVersion("lcdaGRASP"))
  dir <- system.file("extdata", package = "lcdaGRASP")
  skip_if(!nzchar(dir) || length(list.files(dir, pattern = "\\.rds$")) == 0,
          "no extdata datasets installed")
  fs <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  bad <- vapply(fs, function(f) {
    m <- readRDS(f)$meta
    v <- if (is.null(m$pkg_version)) NA_character_ else m$pkg_version
    if (is.na(v) || v != ver) paste0(basename(f), " (", v, ")") else NA_character_
  }, character(1))
  bad <- bad[!is.na(bad)]
  expect_true(length(bad) == 0,
    info = paste0("RDS whose meta$pkg_version != ", ver, ": ",
                  paste(bad, collapse = ", ")))
})

test_that("lcda_provenance reports the current version and a SHA-256 checksum", {
  dir <- system.file("extdata", package = "lcdaGRASP")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "repro_summary.rds")),
          "extdata not installed")
  p <- lcda_provenance("repro_summary")
  expect_s3_class(p, "data.frame")
  expect_identical(p$pkg_version, as.character(utils::packageVersion("lcdaGRASP")))
  # checksum present (64 hex chars) and matching the shipped SHA256SUMS manifest
  skip_if(!file.exists(file.path(dir, "SHA256SUMS")), "no SHA256SUMS manifest")
  expect_match(p$sha256, "^[0-9a-f]{64}$")
})

test_that("lcda_provenance returns NULL for an unknown dataset", {
  expect_null(suppressMessages(lcda_provenance("does_not_exist_xyz")))
})
