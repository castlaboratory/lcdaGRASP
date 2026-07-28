# Release procedure

Single source of truth for the package **and** the JOSS paper (`paper/`).

1. **Version.** Bump `Version` and `Date` in `DESCRIPTION`; add a `NEWS.md` entry.
   Do this **in the release commit**, not after tagging: anything merged after a
   tag while `DESCRIPTION` still shows the tagged version makes two different
   contents installable under one version number (this is what happened between
   `v0.3.1` and `0.3.2`).
2. **Refresh the data manifests** (no simulation is re-run):
   ```bash
   Rscript data-raw/999_manifest.R    # -> inst/extdata/{MANIFEST.csv,SHA256SUMS}
   ```
   Data is **not** regenerated just because the version changed. `generated_by`
   in each `.rds` is the version that actually produced the numbers and stays
   as it is; `first_release` in `MANIFEST.csv` records the release a file first
   shipped in. Regenerate data only when the algorithms change:
   ```bash
   Rscript data-raw/run_all_data.R           # core panel (FORCE=1 to overwrite)
   # plus the individual heavy scripts (11,12,25,35,45,55,90,95,96,97,98,99,100)
   ```
   Confirm consistency either way:
   ```r
   testthat::test_file("tests/testthat/test-data-versions.R")  # must pass
   ```
   (The test fails if a dataset claims a version newer than `DESCRIPTION`, is
   missing from either manifest, or its checksum does not match.)
3. **Check.** `R CMD build .` then `R CMD check --as-cran lcdaGRASP_*.tar.gz`
   (expect only the benign "new submission" / local-HTML-Tidy notes).
4. **Source archive (never Finder/zip).** Use git so junk is excluded:
   ```bash
   git archive --format=zip --prefix=lcdaGRASP/ -o lcdaGRASP-vX.Y.Z.zip HEAD
   unzip -l lcdaGRASP-vX.Y.Z.zip | grep -E '__MACOSX|\.DS_Store|\.o$|\.so$|/\.git' \
     && echo "JUNK FOUND" || echo "clean"
   ```
5. **Checksums.** Already written by step 2. Verify:
   `cd inst/extdata && shasum -a 256 -c SHA256SUMS` (all `OK`).
6. **Tag the merge commit**, after the version bump is on `main`:
   `git tag -a vX.Y.Z -m "lcdaGRASP X.Y.Z"; git push origin vX.Y.Z`.
   Then check nothing drifted: `git diff --stat vX.Y.Z..main` must be empty
   until the next version bump lands.
7. **CRAN / JOSS.** Submit the tarball; the JOSS paper builds from `paper/`.

## End-to-end reproduction (one command)

```bash
Rscript -e 'remotes::install_local(".", dependencies = TRUE)'
Rscript data-raw/run_all_data.R   # + the individual scripts; needs data-raw/benchmarks and .venv-lfr
```
