# LCDA construction - variant 1 (centrality computed once) or 2 (adaptive).

LCDA construction - variant 1 (centrality computed once) or 2
(adaptive).

## Usage

``` r
lcda_construct(
  csr,
  alpha_c,
  alpha_s,
  variant = 1,
  centrality = "eigen",
  similarity = "hpi",
  verbose = FALSE,
  cent = NULL
)
```

## Arguments

- csr:

  CSR object from the internal converter.

- alpha_c:

  numeric in \[0,1\] - centrality RCL parameter.

- alpha_s:

  numeric in \[0,1\] - similarity RCL parameter.

- variant:

  1 (static centrality) or 2 (recomputed each iteration).

- centrality:

  one of "eigen", "betweenness", "closeness".

- similarity:

  one of "hpi", "dice", "jaccard".

- verbose:

  logical; emit a cli trace of the construction.

- cent:

  optional numeric vector of length \`csr\$n\`: a precomputed global
  centrality, used only by \`variant = 1\`. The global centrality is
  deterministic per graph, so the outer GRASP loops compute it once and
  pass it here instead of recomputing it at every one of the \`B\`
  iterations; results are identical either way.

## Value

list(membership, leaders, d) - membership a 1-based vector, leaders a
1-based integer vector of leader indices, d the number of communities.

## See also

\[as_csr()\] to build \`csr\`; \[lcda_repair()\] and
\[lcda_local_search()\] for the remaining pipeline stages;
\[lcda_grasp()\] for the all-in-one driver.

## Examples

``` r
g <- igraph::make_graph("Zachary")
csr <- as_csr(g)
sol <- lcda_construct(csr, alpha_c = 0.1, alpha_s = 0.3)
sol <- lcda_repair(csr, sol)
sol <- lcda_local_search(csr, sol)
c(communities = sol$d, leaders = length(sol$leaders))
#> communities     leaders 
#>           3           3 
```
