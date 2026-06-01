# Convert an igraph graph to the CSR representation used by the kernels

The composable building blocks \[lcda_construct()\], \[lcda_repair()\],
and \[lcda_local_search()\] operate on a compressed-sparse-row (CSR)
view of the graph rather than on the igraph object directly.
\`as_csr()\` builds that view: an undirected simple graph becomes a list
of \`indptr\`/\`indices\`/\`weights\` (0-based, as the C++ kernels
expect) together with the originating igraph. A numeric \`weight\` edge
attribute, if present, is carried through (the kernels reduce exactly to
the unweighted case when all weights are 1).

## Usage

``` r
as_csr(graph)
```

## Arguments

- graph:

  an undirected \[igraph::igraph\] object. Multi-edges are collapsed to
  a simple graph; directed graphs are not supported.

## Value

a CSR list with elements \`n\`, \`indptr\`, \`indices\`, \`weights\`,
and \`igraph\`, suitable as the \`csr\` argument of
\[lcda_construct()\], \[lcda_repair()\], and \[lcda_local_search()\].

## Examples

``` r
g <- igraph::make_graph("Zachary")
csr <- as_csr(g)
str(csr[c("n", "indptr", "indices")])
#> List of 3
#>  $ n      : num 34
#>  $ indptr : int [1:35] 0 16 25 35 41 44 48 52 56 61 ...
#>  $ indices: int [1:156] 1 2 3 4 5 6 7 8 10 11 ...
```
