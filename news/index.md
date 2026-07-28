# Changelog

## lcdaGRASP 0.3.2

This release adds the reporting surface
([`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md))
and the community-and-leader maps, and gives a version number to the
four datasets that were merged **after** the `v0.3.1` tag. Between
`v0.3.1` (commit `09cb457`) and this release, six commits added new
shipped data while `DESCRIPTION` still read `0.3.1`, so two different
contents were both installable as “0.3.1” and anyone installing from the
`v0.3.1` release did not get the data behind Section 5.5 of the
companion paper. **No algorithm behaviour changes here**: the search
itself is untouched, and every addition below is reporting and plotting.

### New features

- [`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md)
  — **one metric surface for everything the paper reports**. Given a
  fitted result (and, optionally, a ground truth) it returns a tidy
  table with modularity `Q`, the NCE leader score in both its global and
  its community-conditioned form, recovery indices (NMI, ARI, Rand, VI,
  split-join), the number and size distribution of the communities,
  wall-clock runtime, and the search diagnostics (trace dispersion, best
  iteration, and how often the lexicographic tie-break on `H` was
  decisive). `level = "community"` and `level = "leader"` give the
  per-community and per-leader breakdowns, including each community’s
  additive contribution to `Q` (they sum exactly to `Q`), its
  conductance and purity, and each leader’s within-community degree rank
  and percentile. Baselines that carry no leaders (an `igraph`
  `communities` object, a bare membership vector) are accepted and
  scored the same way, with derived leaders clearly flagged, so a
  paper-style comparison table can be assembled from one function.
  ([\#45](https://github.com/castlaboratory/lcdaGRASP/issues/45))

  `node_score =` reproduces the paper’s external leader validation on
  any graph: pass a per-vertex outside signal and the leaders are ranked
  within their own community by it, with the mean percentile and the
  top-1 / top-3 hit rates summarised in a `leader_score` scope. The
  participation coefficient (the statistic the paper uses to
  characterise bridge nodes) and the ensemble overlap fraction are
  reported as well.

- **Community-and-leader maps drawn by the package**
  ([\#44](https://github.com/castlaboratory/lcdaGRASP/issues/44)):

  - [`lcda_plot_communities()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_plot_communities.md)
    runs the whole pipeline on a bare graph — detect communities,
    designate leaders, render the figure — in one call, and returns the
    fitted result so it can be piped into
    [`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md).
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
    `lcda_grasp_result`, `lcda_gr_result` and `lcda_ecg_result` plot a
    fitted result directly.
  - [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
    methods return the same figure as a `ggplot` object ( stays a
    *Suggested* dependency; the base-graphics path has no extra
    requirements).
  - [`plot_partition()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot_partition.md)
    gained `layout`, `leader_labels`, `legend`, `legend_max`,
    `vertex_size`, `leader_size`, `palette`, `shade_edges`,
    `mark_communities` and `main`. It now labels leaders, shades
    intra-community edges by community, greys the inter-community ones,
    draws a community/leader legend, and returns the layout and colours
    used so a companion figure can reuse them. The previous call
    signature still works.

- [`lcda_grasp()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_grasp.md),
  [`lcda_gr()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_gr.md)
  and
  [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  results now carry `elapsed` (wall-clock seconds, the runtime column of
  the paper’s timing tables) and `graph` (the simplified graph the
  kernels actually saw). Because the result is self-contained,
  `lcda_metrics(res)` and `plot(res)` need no second argument. The
  [`print()`](https://rdrr.io/r/base/print.html) methods report the
  runtime and point at both. `elapsed` times the search only, starting
  after the
  [`as_csr()`](https://castlaboratory.github.io/lcdaGRASP/reference/as_csr.md)
  conversion, so it is 0.3-2% smaller than timing the whole call from
  outside.

Two reporting choices worth knowing, both made so that no number is
published where none was measured:

- The lexicographic metrics (`lex_decisive_n`, `lex_decisive_pct`) are
  emitted **only for
  [`lcda_gr()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_gr.md)**,
  the one algorithm that instruments the tie-break. For
  [`lcda_grasp()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_grasp.md)
  and
  [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  they are absent rather than reported as a measured-looking `0`.
- `internal_density` is always **structural** (a proportion of the
  possible pairs, hence always in `[0, 1]`) even on a weighted graph,
  where `internal_edges` and `boundary_edges` are weight sums. `Q` and
  `conductance` remain weight-aware.

### Notes

- **Fitted results now hold a reference to their graph.** This is what
  makes `lcda_metrics(res)` and `plot(res)` work without a second
  argument, but it matters if you
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) a pool of results:
  the graph dominates the serialised size (measured on `n = 20000`: 1.60
  MB of a 1.69 MB result). In memory it is free (R shares the object);
  on disk it is not. Drop it with `res$graph <- NULL` before saving a
  large pool, and pass the graph explicitly to
  [`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md)
  / [`plot()`](https://rdrr.io/r/graphics/plot.default.html) afterwards.
- `res$graph` is the **simplified** graph (multi-edges collapsed, as
  [`as_csr()`](https://castlaboratory.github.io/lcdaGRASP/reference/as_csr.md)
  does), not necessarily the object you passed in. Metrics are computed
  on it deliberately, so that `Q` matches the value the objective
  actually optimised.

### New datasets (added after the `v0.3.1` tag)

- `realnet_amazon` — large real network with ground truth:
  Amazon-Computers co-purchase graph (recovery and runtime).
- `vnmi_nprime` — VNMI candidate-subset size `n'` sensitivity sweep:
  modularity and local-search time vs `n'` across density regimes.
- `realnet_coauthor` — Coauthor-Physics (largest connected component, n
  = 34,493, 5 fields), the paper’s headline real network.
- `openalex_leaders` — leader validation against an external citation
  signal on an OpenAlex co-authorship graph.

The other 24 datasets carry the **same numbers** as the ones shipped in
`v0.3.1`: their `results` objects are byte-for-byte identical (verified
by serialising `$results` for all 28 datasets before and after the
metadata repair below and comparing SHA-256: 28/28 unchanged), as is
every `meta` field except `network_source`. They are **not**
bit-identical as files, and their SHA-256 checksums in
`inst/extdata/SHA256SUMS` have changed accordingly: the one provenance
field that was false has been corrected in place (see *Fixes* below). No
simulation was re-run for this release.

### Provenance model

Dataset provenance now records three distinct facts instead of
overloading one version string:

- `generated_by` — the package version whose code produced the numbers.
  This is a historical fact and is never rewritten on a version bump.
  All 28 datasets currently report `0.3.1`, which is the truth: none of
  them were regenerated for this release.

- `first_release` — the release in which a file first shipped, recorded
  in the new `inst/extdata/MANIFEST.csv` (24 datasets at `0.3.1`, 4 at
  `0.3.2`).

- `shipped_in` — the version of the installed copy, i.e.
  `packageVersion("lcdaGRASP")`.

- [`lcda_provenance()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_provenance.md)
  returns these three columns plus `generated_on` and `sha256`. The old
  `pkg_version` column is renamed `generated_by`, because it never meant
  “the version you installed”.

- `data-raw/999_manifest.R` regenerates `MANIFEST.csv` and `SHA256SUMS`
  without re-running any simulation; `first_release` is sticky once
  recorded.

- `data-raw/00_helpers.R` now stamps the git commit (`meta$git_commit`)
  into newly generated datasets, so future data pins an exact source
  state rather than a possibly-unbumped version string.

- The data-version test no longer demands that every dataset carry the
  current package version (an invariant that can only be met by
  rewriting metadata, or by re-running hours of simulation on every
  bump). It now enforces the honest invariants: no dataset generated by
  a *future* version, full coverage by both manifests, metadata
  agreement, and matching SHA-256 checksums.

- `MANIFEST.csv` is itself listed in `SHA256SUMS`, so the
  release-identity manifest is covered by the same integrity check as
  the data it describes; and `data-raw/999_manifest.R` now warns instead
  of silently falling back to the current version when it cannot reach
  git to derive `first_release`.

### Fixes

- **`meta$network_source` was a false constant in all 28 shipped
  datasets.** `data-raw/00_helpers.R::make_meta()` hard-coded it to
  `"Newman netdata collection (websites.umich.edu/~mejn) + igraph::make_graph('Zachary')"`
  and collected every other argument into `extra`. Every dataset
  therefore claimed the five classical benchmarks as its origin —
  including the synthetic LFR sweeps generated with `networkx`
  (`lfr_robustness`, `largescale`, `lcda_ecg`, …), the purely synthetic
  SBM study (`pool_sensitivity`), the Shchur et al. (2018) real networks
  (`realnet_coauthor`, `realnet_amazon`) and the OpenAlex extraction
  (`openalex_leaders`). The three generators that *did* pass a correct
  `network_source` passed it through `...`, where it was swallowed by
  `extra` and never reached the field readers see.

  Fixed at the source: `make_meta()` now takes `network_source` as a
  required formal argument with **no default**, aborts if it is missing,
  empty, or passed through `...`, and every `data-raw/*.R` generator
  supplies the network it actually ran on. The 28 shipped `.rds` files
  were repaired in place by the one-off
  `data-raw/998_fix_network_source.R`, which rewrites *only*
  `meta$network_source` (and drops the shadowed
  `meta$extra$network_source` copy) and aborts if anything else moves.
  Correcting a field that was **false** is not the same as restamping
  one that was true: `generated_by` / `generated_on` are untouched,
  exactly as in the 0.3.2 provenance model above.

  `tests/testthat/test-data-network-source.R` locks the invariants in:
  every dataset declares a non-empty source, the collection is not one
  constant, no synthetic or external dataset cites the benchmark
  collection, every benchmark-based dataset names a benchmark, no
  `network_source` survives in `meta$extra`, and every installed dataset
  is covered by the guard.

- Vignette `robustness-and-limits`: the two uncertainty ribbons (NMI vs
  mixing, and leader spreading advantage) were silently missing. Both
  [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
  calls overwrote a column with its own mean before taking
  [`sd()`](https://rdrr.io/r/stats/sd.html) of it, so the dispersion
  evaluated to `NA` and
  [`geom_ribbon()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html)
  dropped every row. The dispersion is now computed before the mean, and
  the bands render again.

## lcdaGRASP 0.3.1

### Fixes

- [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  now passes the chosen `centrality` to the repair and local-search
  steps (previously the consensus leader was always eigenvector-based).
- `seed = NULL` (and `NA`) is accepted and means “leave the RNG
  untouched”; invalid seeds are rejected with a clear message.
- `data-raw/97_lcda_ecg.R` and the scripts that reused it now call the
  exported
  [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  instead of a private re-implementation, so the cached `lcda_ecg.rds`
  reflects the canonical function (incl. `Q_consensus_weighted`).

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
