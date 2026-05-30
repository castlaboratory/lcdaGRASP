# Betweenness centrality - delegates to igraph::betweenness for now. Tagged for a future native Rcpp implementation (Brandes 2001).

Betweenness centrality - delegates to igraph::betweenness for now.
Tagged for a future native Rcpp implementation (Brandes 2001).

## Usage

``` r
centrality_betweenness(csr)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

## Value

numeric vector of per-node centrality scores.
