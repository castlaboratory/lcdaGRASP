# lcdaGRASP (development version)

## New features

* `lcda_metrics()` — **one metric surface for everything the paper reports**.
  Given a fitted result (and, optionally, a ground truth) it returns a tidy
  table with modularity `Q`, the NCE leader score in both its global and its
  community-conditioned form, recovery indices (NMI, ARI, Rand, VI,
  split-join), the number and size distribution of the communities, wall-clock
  runtime, and the search diagnostics (trace dispersion, best iteration, and
  how often the lexicographic tie-break on `H` was decisive).
  `level = "community"` and `level = "leader"` give the per-community and
  per-leader breakdowns, including each community's additive contribution to
  `Q` (they sum exactly to `Q`), its conductance and purity, and each leader's
  within-community degree rank and percentile. Baselines that carry no leaders
  (an `igraph` `communities` object, a bare membership vector) are accepted and
  scored the same way, with derived leaders clearly flagged, so a paper-style
  comparison table can be assembled from one function. (#45)

  `node_score =` reproduces the paper's external leader validation on any
  graph: pass a per-vertex outside signal and the leaders are ranked within
  their own community by it, with the mean percentile and the top-1 / top-3 hit
  rates summarised in a `leader_score` scope. The participation coefficient
  (the statistic the paper uses to characterise bridge nodes) and the ensemble
  overlap fraction are reported as well.

* **Community-and-leader maps drawn by the package** (#44):
  - `lcda_plot_communities()` runs the whole pipeline on a bare graph — detect
    communities, designate leaders, render the figure — in one call, and
    returns the fitted result so it can be piped into `lcda_metrics()`.
  - `plot()` methods for `lcda_grasp_result`, `lcda_gr_result` and
    `lcda_ecg_result` plot a fitted result directly.
  - `ggplot2::autoplot()` methods return the same figure as a `ggplot` object
    (\pkg{ggplot2} stays a *Suggested* dependency; the base-graphics path has
    no extra requirements).
  - `plot_partition()` gained `layout`, `leader_labels`, `legend`,
    `legend_max`, `vertex_size`, `leader_size`, `palette`, `shade_edges`,
    `mark_communities` and `main`. It now labels leaders, shades
    intra-community edges by community, greys the inter-community ones, draws
    a community/leader legend, and returns the layout and colours used so a
    companion figure can reuse them. The previous call signature still works.

* `lcda_grasp()`, `lcda_gr()` and `lcda_ecg()` results now carry `elapsed`
  (wall-clock seconds, the runtime column of the paper's timing tables) and
  `graph` (the simplified graph the kernels actually saw). Because the result
  is self-contained, `lcda_metrics(res)` and `plot(res)` need no second
  argument. The `print()` methods report the runtime and point at both.

## Notes

* The package version is deliberately left at 0.3.1: the shipped
  `inst/extdata` datasets record the version that generated them and
  `tests/testthat/test-data-versions.R` enforces the match, so the bump belongs
  to a release commit that also regenerates the data.

# lcdaGRASP 0.3.1

## Fixes

* `lcda_ecg()` now passes the chosen `centrality` to the repair and local-search
  steps (previously the consensus leader was always eigenvector-based).
* `seed = NULL` (and `NA`) is accepted and means "leave the RNG untouched";
  invalid seeds are rejected with a clear message.
* `data-raw/97_lcda_ecg.R` and the scripts that reused it now call the exported
  `lcda_ecg()` instead of a private re-implementation, so the cached
  `lcda_ecg.rds` reflects the canonical function (incl. `Q_consensus_weighted`).

# lcdaGRASP 0.3.0

## New features

* `lcda_ecg()` — **ensemble-consensus community and leader detection**. Turns
  the GRASP pool (the diverse partitions produced across iterations, previously
  discarded) into edge co-association weights in the spirit of Ensemble
  Clustering for Graphs (Poulin & Théberge, 2019), re-clusters for a consensus
  partition, and designates one leader per community from the pool's
  leader-designation frequencies. On canonical LFR benchmarks it recovers the
  planted structure on par with ECG and outperforms Leiden, with its advantage
  concentrated at high mixing, while retaining the joint leader output and
  adding a per-node confidence map.
  With `overlap = TRUE` it additionally returns overlapping community
  memberships (and the bridge nodes) derived from the soft co-association.
  Returns an `lcda_ecg_result` object with a `print()` method.

## Notes

* Companion-paper reproduction now includes a canonical LFR sweep (recovery
  vs Leiden/Louvain/ECG), a design-of-experiments study of recovery
  (the `alpha_s` modularity-vs-recovery trade-off), and leader-utility analyses
  (information spreading, community representation vs a two-stage pipeline).

# lcdaGRASP 0.2.0

* tidyverse-style refactor with `cli` logging and a `verbose` argument
  throughout; statistics return tibbles.
* `lcda_data()` accessor for precomputed datasets under `inst/extdata/`.
* pkgdown site with nine vignettes.

# lcdaGRASP 0.1.0

* Initial version: LCDA-GRASP (fixed-parameter) and LCDA-GR (Reactive)
  algorithms with performance-critical C++ kernels via Rcpp; modularity,
  NCE leader score, centralities, and statistical comparison helpers.
