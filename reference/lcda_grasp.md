# LCDA-GRASP: fixed-parameter variant.

LCDA-GRASP: fixed-parameter variant.

## Usage

``` r
lcda_grasp(
  g,
  alpha_c = 0.1,
  alpha_s = 0.3,
  variant = 1,
  B = 50,
  centrality = "eigen",
  similarity = "hpi",
  verbose = FALSE,
  seed = NA_integer_
)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- alpha_c, alpha_s:

  RCL parameters.

- variant:

  1 or 2.

- B:

  number of GRASP iterations.

- centrality, similarity:

  metric names.

- verbose:

  logical; show a cli progress bar and a final summary.

- seed:

  integer RNG seed, or \`NA\` to leave the RNG untouched.

## Value

an object of class \`lcda_grasp_result\`: best partition, Q/H traces,
and the parameters used.

## Details

Weighted graphs: a numeric \`weight\` edge attribute is honoured by the
modularity objective and the local search, but the construction
(similarity and centrality) and the NCE leader score remain
\*structural\* (unweighted).

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_grasp(g, B = 20, seed = 1)
res$best$Q
#> [1] 0.4020381
res$best$leaders
#> [1] 34  1  6
```
