# Provenance of a shipped dataset (version, date, checksum)

Reports the reproducibility metadata for one of the \[lcda_data()\]
datasets: the package version that generated it, the generation date,
and its SHA-256 checksum looked up from the \`SHA256SUMS\` manifest
shipped alongside the data (so it matches \`shasum -a 256 -c
inst/extdata/SHA256SUMS\`). Use it in analyses and vignettes to pin
exactly which version of a dataset a result came from. The version is
also enforced by the package's own tests, so a shipped dataset always
matches the installed \`packageVersion("lcdaGRASP")\`.

## Usage

``` r
lcda_provenance(name)
```

## Arguments

- name:

  dataset name (without the \`.rds\` extension), as in \[lcda_data()\].

## Value

a one-row data frame with columns \`dataset\`, \`pkg_version\`,
\`generated_on\`, and \`sha256\`; or \`NULL\` (invisibly) if the dataset
is not installed.

## See also

\[lcda_data()\]

## Examples

``` r
lcda_provenance("repro_summary")
#>         dataset pkg_version generated_on
#> 1 repro_summary       0.3.1   2026-05-31
#>                                                             sha256
#> 1 0a140370d6b7ce2272808592c0c82d3945a4ae21b328638c2cc3ba18f77cdbbe
```
