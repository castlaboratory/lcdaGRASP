# Changelog

## lcdaGRASP 0.3.0

### New features

- [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  — **ensemble-consensus community and leader detection**. Turns the
  GRASP pool (the diverse partitions produced across iterations,
  previously discarded) into edge co-association weights in the spirit
  of Ensemble Clustering for Graphs (Poulin & Théberge, 2019),
  re-clusters for a consensus partition, and designates one leader per
  community from the pool’s leader-designation frequencies. On canonical
  LFR benchmarks it recovers the planted structure on par with ECG and
  outperforms Leiden, with its advantage concentrated at high mixing,
  while retaining the joint leader output and adding a per-node
  confidence map. With `overlap = TRUE` it additionally returns
  overlapping community memberships (and the bridge nodes) derived from
  the soft co-association. Returns an `lcda_ecg_result` object with a
  [`print()`](https://rdrr.io/r/base/print.html) method.

### Notes

- Companion-paper reproduction now includes a canonical LFR sweep
  (recovery vs Leiden/Louvain/ECG), a design-of-experiments study of
  recovery (the `alpha_s` modularity-vs-recovery trade-off), and
  leader-utility analyses (information spreading, community
  representation vs a two-stage pipeline).

## lcdaGRASP 0.2.0

- tidyverse-style refactor with `cli` logging and a `verbose` argument
  throughout; statistics return tibbles.
- [`lcda_data()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_data.md)
  accessor for precomputed datasets under `inst/extdata/`.
- pkgdown site with nine vignettes.

## lcdaGRASP 0.1.0

- Initial version: LCDA-GRASP (fixed-parameter) and LCDA-GR (Reactive)
  algorithms with performance-critical C++ kernels via Rcpp; modularity,
  NCE leader score, centralities, and statistical comparison helpers.
