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
  p_floor = 0.05,
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

- p_floor:

  minimum probability mass reserved across the pool at each refresh,
  spread uniformly so every pair keeps \`p_k \>= p_floor/m \> 0\`. This
  prevents a pair that happened not to be sampled in a block from being
  permanently excluded (and matches the \`p_k \>= delta \> 0\` premise
  of Proposition 6). Set to 0 to recover the raw proportional rule.

- verbose:

  logical; show a cli progress bar and a final summary.

- seed:

  integer RNG seed, or \`NA\` to leave the RNG untouched.

## Value

an object of class \`lcda_gr_result\`: best partition, traces, the
reactive pool state, the H-decisive iterations, the wall-clock
\`elapsed\` time in seconds, and the (simplified) input \`graph\`, so
that \[lcda_metrics()\] and \[plot()\] can be called on the result
alone.

## Details

Weighted graphs: a numeric \`weight\` edge attribute is honoured by the
modularity objective and the local search, but the construction
(similarity and centrality) and the NCE leader score remain
\*structural\* (unweighted).

## See also

\[lcda_metrics()\], \[plot.lcda_gr_result()\].

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_gr(g, B = 30, seed = 1)
res$best$Q
#> [1] 0.4197896
length(res$lex_decisive_iters)  # how often H broke a Q-tie
#> [1] 0
```
