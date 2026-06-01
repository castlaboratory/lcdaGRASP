# Repair: ensure every community has exactly one leader by recomputing centrality within the community and designating the top-scoring node.

Repair: ensure every community has exactly one leader by recomputing
centrality within the community and designating the top-scoring node.

## Usage

``` r
lcda_repair(csr, construction, centrality = "eigen", verbose = FALSE)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

- construction:

  a list as returned by \[lcda_construct()\].

- centrality:

  leader-selection centrality (\`"eigen"\`, \`"betweenness"\`, or
  \`"closeness"\`); should match the one used in construction so the
  final leader is consistent with the chosen measure.

- verbose:

  logical; emit a cli trace.

## Value

the input \`construction\` with its \`leaders\` field rebuilt.

## See also

\[as_csr()\], \[lcda_construct()\], \[lcda_local_search()\].

## Examples

``` r
csr <- as_csr(igraph::make_graph("Zachary"))
sol <- lcda_construct(csr, alpha_c = 0.1, alpha_s = 0.3)
sol <- lcda_repair(csr, sol)
length(sol$leaders)            # one leader per community
#> [1] 3
```
