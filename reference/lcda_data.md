# Precomputed simulation datasets

The package ships the results of its publication-grade simulation panel
as compressed \`.rds\` files under \`inst/extdata/\`. The vignettes read
these so they never re-run a simulation at build time. Regenerate them
with the scripts in \`data-raw/\` (see \`data-raw/run_all_data.R\`).

## Usage

``` r
lcda_data(name = NULL)
```

## Arguments

- name:

  dataset name (without the \`.rds\` extension). If \`NULL\` (default),
  returns a tibble listing the available datasets.

## Value

If \`name\` is \`NULL\`, a tibble of available datasets. Otherwise the
\`list(results, meta)\` bundle for that dataset, or \`NULL\` with a
message if it has not been generated yet.

## Details

Available datasets (name -\> contents):

- \`repro_benchmarks\`, \`repro_summary\`:

  Benchmark reproduction on the five networks vs Louvain/Leiden and
  literature best-known Q.

- \`doe_screening\`, \`doe_rsm\`:

  Design-of-experiments runs: factorial screening and a central
  composite design in \`(alpha_c, alpha_s)\`.

- \`lfr_robustness\`:

  Q and ARI across the LFR-like mixing sweep.

- \`prop6_summary\`, \`prop6_trajectories\`:

  Reactive-update concentration.

- \`pool_sensitivity\`:

  \`(m, y, B)\` grid for LCDA-GR.

- \`degeneracy\`:

  Near-best partition counts and pairwise NMI.

- \`eda_replicates\`:

  Replicated \`(Q, H)\` for exploratory analysis.

- \`nce_alternatives\`:

  Global vs community-conditioned NCE leaders.

- \`lcda_ecg\`, \`overlap_lcda_ecg\`, \`stability_pool\`,
  \`consensus_leader_test\`:

  LCDA-ECG ensemble consensus: LFR recovery, overlapping/bridge nodes,
  pool-stability stopping rule, and the consensus-vs-central leader
  comparison.

- \`largescale\`:

  Recovery and runtime up to n = 5e4 (synthetic LFR).

- \`realnet_amazon\`:

  Large real network with ground truth: modularity Q and recovery
  (NMI/ARI vs product categories) on Amazon-Computers.

- \`gnn_baseline\`:

  Graph-auto-encoder baseline (true-k and auto-k) vs classical methods
  on LFR.

- \`weighted_demo\`:

  Weighted-graph demonstration (weights aid recovery).

- \`blogs_table9\`, \`blogs_timing\`:

  Weighted Political Blogs reproduction (mean/max Q, timings).

- \`doe_lfr_recovery\`:

  DoE screening/RSM scored by LFR recovery.

- \`leader_utility\`, \`leader_vs_twostage\`,
  \`leader_vs_twostage_lfr\`:

  Downstream leader utility (IC/LT spread, coverage) vs two-stage
  detect-then-centrality pipelines.

## Examples

``` r
lcda_data()                       # list what is available
#> # A tibble: 25 × 1
#>    dataset              
#>    <chr>                
#>  1 blogs_table9         
#>  2 blogs_timing         
#>  3 consensus_leader_test
#>  4 degeneracy           
#>  5 doe_lfr_recovery     
#>  6 doe_rsm              
#>  7 doe_screening        
#>  8 eda_replicates       
#>  9 gnn_baseline         
#> 10 largescale           
#> # ℹ 15 more rows
if (FALSE) { # \dontrun{
d <- lcda_data("repro_summary")
d$meta$limitations
head(d$results)
} # }
```
