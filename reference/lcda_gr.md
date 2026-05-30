# LCDA-GR: Reactive variant with self-tuning of (alpha_c, alpha_s). Implements Algorithm 4 of the paper.

LCDA-GR: Reactive variant with self-tuning of (alpha_c, alpha_s).
Implements Algorithm 4 of the paper.

## Usage

``` r
lcda_gr(
  g,
  variant = 1,
  B = 150,
  centrality = "eigen",
  similarity = "hpi",
  m = 20,
  y = NULL,
  alpha_c_range = c(0.1, 0.9),
  alpha_s_range = c(0.1, 0.5),
  verbose = FALSE,
  seed = NA_integer_
)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- variant:

  construction variant, 1 or 2.

- B:

  number of GRASP iterations.

- centrality, similarity:

  metric names.

- m:

  pool size (default 20, per paper).

- y:

  refresh period (default 3m).

- alpha_c_range, alpha_s_range:

  bounds for the uniform initial pool.

- verbose:

  logical; show a cli progress bar and a final summary.

- seed:

  integer RNG seed, or \`NA\` to leave the RNG untouched.

## Value

an object of class \`lcda_gr_result\`: best partition, traces, the
reactive pool state, and the H-decisive iterations.

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_gr(g, B = 30, seed = 1)
res$best$Q
#> [1] 0.4197896
length(res$lex_decisive_iters)  # how often H broke a Q-tie
#> [1] 0
```
