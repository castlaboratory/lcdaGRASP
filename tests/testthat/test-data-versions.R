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
