## Submission

New submission of lcdaGRASP 0.3.1.

## R CMD check results

Local check (`R CMD check --as-cran`, macOS, R 4.6.0):
**0 errors | 0 warnings | 2 notes.**

The only package-level NOTE is the expected *"New submission"* (CRAN incoming
feasibility) for a package not yet on CRAN.

Any additional NOTEs seen locally are properties of the **local toolchain**, not
of the package, and do not occur on CRAN's check machines — e.g. an outdated
HTML Tidy binary ("Skipping checking HTML validation"), or clock-skew /
`xcrun`-cache messages on some macOS setups.

The source tarball is ~400 KB (the Python virtual environment used only by the
`data-raw/` generators is `.Rbuildignore`d and not shipped).

## Test environments

* local: macOS (aarch64), R 4.6.0 — `devtools::test()`: 111 pass, 0 fail.
* GitHub Actions: Ubuntu, macOS, Windows (R release) — all green.
* win-builder (devel and release): pending.

## Downstream dependencies

None — this is a new package.
