# Repair: ensure every community has exactly one leader by recomputing centrality within the community and designating the top-scoring node.

Repair: ensure every community has exactly one leader by recomputing
centrality within the community and designating the top-scoring node.

## Usage

``` r
lcda_repair(csr, construction, verbose = FALSE)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

- construction:

  a list as returned by \[lcda_construct()\].

- verbose:

  logical; emit a cli trace.

## Value

the input \`construction\` with its \`leaders\` field rebuilt.
