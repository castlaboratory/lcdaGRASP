# data-raw/999_manifest.R
#
# Regenerate the two provenance manifests that ship with the datasets:
#
#   inst/extdata/SHA256SUMS   content identity  (what the bytes are)
#   inst/extdata/MANIFEST.csv release identity  (which package version made a
#                                                dataset, and which release
#                                                first shipped it)
#
# Why two versions per dataset. `generated_by` is the package version whose
# code produced the numbers; it is a historical fact and is NEVER rewritten
# when the package version is bumped. `first_release` is the package version
# in whose release the file first appeared. They differ whenever a release
# ships data it did not itself regenerate, which is the normal case: the
# expensive generators (LFR / real-network / OpenAlex, via the .venv-lfr
# Python venv) are re-run only when the algorithms change, not on every
# version bump.
#
# `first_release` is sticky: once recorded it is carried over verbatim from
# the existing MANIFEST.csv. Only datasets not yet in the manifest get a
# value, derived from git (the first tag containing the commit that added the
# file) and falling back to the current DESCRIPTION version for data that has
# not shipped in any release yet.
#
# Run from the package root (no simulation is re-run; this only reads files):
#   Rscript data-raw/999_manifest.R

.pkg_dir <- local({
  cands <- c(".", "lcdaGRASP", "..", "../lcdaGRASP")
  hit <- NULL
  for (d in cands) {
    desc <- file.path(d, "DESCRIPTION")
    if (file.exists(desc) &&
        any(grepl("^Package: lcdaGRASP", readLines(desc, warn = FALSE)))) {
      hit <- normalizePath(d); break
    }
  }
  if (is.null(hit)) stop("Could not locate the lcdaGRASP package directory.")
  hit
})

extdata  <- file.path(.pkg_dir, "inst", "extdata")
manifest <- file.path(extdata, "MANIFEST.csv")
sums     <- file.path(extdata, "SHA256SUMS")
cur_ver  <- as.character(read.dcf(file.path(.pkg_dir, "DESCRIPTION"), "Version")[1, 1])

# radix == C-locale byte order, so the manifests are locale-independent and
# match `LC_ALL=C shasum -a 256 *.rds | sort`.
files <- sort(list.files(extdata, pattern = "[.]rds$"), method = "radix")
if (!length(files)) stop("No .rds datasets found in ", extdata)

# ---- previously recorded first_release (sticky) ----------------------------

prev <- if (file.exists(manifest)) {
  utils::read.csv(manifest, stringsAsFactors = FALSE)
} else {
  data.frame(dataset = character(), first_release = character())
}

# ---- git fallback: first tag containing the commit that added the file -----

.first_release_from_git <- function(rel_path) {
  git <- function(...) {
    out <- suppressWarnings(system2("git", c("-C", shQuote(.pkg_dir), ...),
                                    stdout = TRUE, stderr = FALSE))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0) character() else out
  }
  add <- git("log", "--diff-filter=A", "--format=%H", "-1", "--", shQuote(rel_path))
  if (!length(add) || !nzchar(add[1])) return(NA_character_)
  tags <- git("tag", "--contains", add[1], "--sort=v:refname")
  tags <- grep("^v[0-9]", tags, value = TRUE)
  if (!length(tags)) return(NA_character_)
  sub("^v", "", tags[1])
}

# ---- build the manifest ----------------------------------------------------

rows <- do.call(rbind, lapply(files, function(f) {
  name <- sub("[.]rds$", "", f)
  m    <- readRDS(file.path(extdata, f))$meta
  fld  <- function(x) if (is.null(x)) NA_character_ else as.character(x)[1]

  fr <- prev$first_release[match(name, prev$dataset)]
  if (length(fr) != 1 || is.na(fr) || !nzchar(fr)) {
    fr <- .first_release_from_git(file.path("inst", "extdata", f))
    if (is.na(fr)) fr <- cur_ver          # not in any release yet -> this one
  }

  data.frame(
    dataset       = name,
    generated_on  = fld(m$generated_on),
    generated_by  = fld(m$pkg_version),
    first_release = fr,
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(rows, manifest, row.names = FALSE, quote = FALSE)

# ---- checksums (same format as `shasum -a 256 *.rds | sort`) ---------------

hashes <- unname(tools::sha256sum(file.path(extdata, files)))
writeLines(sprintf("%s  %s", hashes, files), sums)

cat(sprintf("MANIFEST.csv: %d datasets (%s)\n", nrow(rows),
            paste(sprintf("%s=%d", names(table(rows$first_release)),
                          as.integer(table(rows$first_release))), collapse = ", ")))
cat(sprintf("SHA256SUMS:   %d entries\n", length(files)))
cat("Verify with: cd inst/extdata && shasum -a 256 -c SHA256SUMS\n")
