# LCDA-ECG: ensemble-consensus community and leader detection.

Builds a pool of `B` randomised LCDA constructions, turns it into edge
co-association weights (ECG-style), re-clusters the reweighted graph for
the consensus partition, and designates one leader per community from
the pool's leader-designation frequencies. Recovers planted structure on
par with ECG and outperforms Leiden (advantage concentrated at high
mixing), while retaining the joint leader output and a node-confidence
map.

## Usage

``` r
lcda_ecg(
  g,
  B = 64,
  w_min = 0.05,
  alpha_c_range = c(0.1, 0.9),
  alpha_s_range = c(0.1, 0.5),
  variant = 1,
  centrality = "eigen",
  similarity = "hpi",
  overlap = FALSE,
  tau = 0.7,
  verbose = FALSE,
  seed = NA_integer_
)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- B:

  pool size (number of GRASP constructions to ensemble).

- w_min:

  ECG floor weight for 2-core edges; off-2-core edges get exactly
  `w_min`. Default 0.05, as in Poulin & Theberge (2019).

- alpha_c_range, alpha_s_range:

  bounds of the uniform RCL parameters sampled per pool member
  (diversification source).

- variant:

  construction variant, 1 or 2.

- centrality, similarity:

  metric names passed to the construction.

- overlap:

  logical; if `TRUE`, also return overlapping community memberships
  derived from the co-association (soft) similarity.

- tau:

  overlap threshold in (0,1\]: a node joins community `c` when its mean
  co-association to `c` reaches `tau` times its home-community affinity.
  Only used when `overlap = TRUE`.

- verbose:

  logical; show a cli progress bar and a final summary.

- seed:

  integer RNG seed, or `NA` to leave the RNG untouched.

## Value

an object of class `lcda_ecg_result`: the consensus `membership`
(1-based), `leaders` (consensus-derived, 1-based), `leaders_central`
(top-eigenvector per community, for comparison), a per-node `confidence`
vector, the leader-designation counts `lead_count`, the input-graph
modularity `Q` (weight-aware; comparable to `lcda_grasp`/`lcda_gr`), and
`Q_consensus_weighted` (modularity under the ECG co-association weights,
i.e. the objective the consensus optimised). When `overlap = TRUE` it
additionally carries `overlap_membership` (a length-n list of the
community ids each node belongs to) and `is_overlap` (logical, the
bridge nodes). It also carries the wall-clock `elapsed` time in seconds
and the (simplified) input `graph`, so that
[`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md)
and
[`plot.lcda_ecg_result()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot.lcda.md)
can be called on the result alone.

## Details

Weighted graphs: a numeric `weight` edge attribute is honoured by the
modularity objective and the local search inside each pool construction,
but the similarity, centrality and NCE leader score remain *structural*
(unweighted). The consensus re-clustering uses the ECG co-association
weights, not the input weights.

## See also

[`lcda_metrics()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_metrics.md),
[`plot.lcda_ecg_result()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot.lcda.md).

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_ecg(g, B = 24, overlap = TRUE, tau = 0.6, seed = 1)
res$Q
#> [1] 0.4197896
which(res$is_overlap)   # bridge nodes
#> [1] 10
```
