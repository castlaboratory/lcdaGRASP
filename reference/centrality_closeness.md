# Closeness centrality - delegates to igraph::closeness. On a disconnected graph, plain closeness is ill-defined across components (nodes in tiny components score spuriously high, distorting leader selection); we fall back to harmonic centrality, the robust generalisation that handles unreachable pairs (Boldi & Vigna 2014). On connected graphs the two coincide in ranking, so this changes nothing for the standard benchmarks.

Closeness centrality - delegates to igraph::closeness. On a disconnected
graph, plain closeness is ill-defined across components (nodes in tiny
components score spuriously high, distorting leader selection); we fall
back to harmonic centrality, the robust generalisation that handles
unreachable pairs (Boldi & Vigna 2014). On connected graphs the two
coincide in ranking, so this changes nothing for the standard
benchmarks.

## Usage

``` r
centrality_closeness(csr)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

## Value

numeric vector of per-node centrality scores.
