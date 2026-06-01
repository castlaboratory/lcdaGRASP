# Release procedure

Single source of truth for the package **and** the JOSS paper (`paper/`).

1. **Version.** Bump `Version` and `Date` in `DESCRIPTION`; add a `NEWS.md` entry.
2. **Regenerate data at the new version.** Run the generators so every
   `inst/extdata/*.rds` is produced by the new version:
   ```bash
   Rscript data-raw/run_all_data.R           # core panel (FORCE=1 to overwrite)
   # plus the individual heavy scripts (11,12,25,35,45,55,90,95,96,97,98,99,100)
   ```
   Confirm consistency:
   ```r
   testthat::test_file("tests/testthat/test-data-versions.R")  # must pass
   ```
   (The test fails if any RDS `meta$pkg_version` differs from `DESCRIPTION`.)
3. **Check.** `R CMD build .` then `R CMD check --as-cran lcdaGRASP_*.tar.gz`
   (expect only the benign "new submission" / local-HTML-Tidy notes).
4. **Source archive (never Finder/zip).** Use git so junk is excluded:
   ```bash
   git archive --format=zip --prefix=lcdaGRASP/ -o lcdaGRASP-vX.Y.Z.zip HEAD
   unzip -l lcdaGRASP-vX.Y.Z.zip | grep -E '__MACOSX|\.DS_Store|\.o$|\.so$|/\.git' \
     && echo "JUNK FOUND" || echo "clean"
   ```
5. **Checksums.** `shasum -a 256 inst/extdata/*.rds | sort > inst/extdata/SHA256SUMS`.
6. **Tag.** `git tag -a vX.Y.Z -m "lcdaGRASP X.Y.Z"; git push origin vX.Y.Z`.
7. **CRAN / JOSS.** Submit the tarball; the JOSS paper builds from `paper/`.

## End-to-end reproduction (one command)

```bash
Rscript -e 'remotes::install_local(".", dependencies = TRUE)'
Rscript data-raw/run_all_data.R   # + the individual scripts; needs data-raw/benchmarks and .venv-lfr
```
