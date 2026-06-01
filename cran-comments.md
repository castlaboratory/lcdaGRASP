## Submission

New submission of lcdaGRASP 0.3.1.

## R CMD check results

Local check (`R CMD check --as-cran`, macOS, R 4.6.0):
**0 errors | 0 warnings | 2 notes.**

The two NOTEs are expected and non-actionable:

1. *"New submission"* (CRAN incoming feasibility) — expected for a package not
   yet on CRAN.
2. *"Skipping checking HTML validation: 'tidy' doesn't look like recent enough
   HTML Tidy."* — this is a property of the **local** machine's outdated HTML
   Tidy binary, not of the package; it does not occur on CRAN's check machines.

The source tarball is ~400 KB (the Python virtual environment used only by the
`data-raw/` generators is `.Rbuildignore`d and not shipped).

## Test environments

* local: macOS (aarch64), R 4.6.0 — `devtools::test()`: 107 pass, 0 fail.
* GitHub Actions: Ubuntu, macOS, Windows (R release) — all green.
* win-builder (devel and release): pending.

## Downstream dependencies

None — this is a new package.
