# Closeness centrality - delegates to igraph::closeness. On a disconnected graph, plain closeness is ill-defined across components (nodes in tiny components score spuriously high, distorting leader selection); we fall back to harmonic centrality, the robust generalisation that handles unreachable pairs (Boldi & Vigna 2014). Note the two are NOT identical even on a connected graph – closeness inverts the mean distance, harmonic averages the inverse distances – so on disconnected inputs (e.g. PolBlogs, and the frequently-disconnected induced subgraphs of variant 2) the selected leader may differ from a per-component closeness; this is a deliberate robustness choice, recorded here so benchmark results are interpreted accordingly.

Closeness centrality - delegates to igraph::closeness. On a disconnected
graph, plain closeness is ill-defined across components (nodes in tiny
components score spuriously high, distorting leader selection); we fall
back to harmonic centrality, the robust generalisation that handles
unreachable pairs (Boldi & Vigna 2014). Note the two are NOT identical
even on a connected graph – closeness inverts the mean distance,
harmonic averages the inverse distances – so on disconnected inputs
(e.g. PolBlogs, and the frequently-disconnected induced subgraphs of
variant 2) the selected leader may differ from a per-component
closeness; this is a deliberate robustness choice, recorded here so
benchmark results are interpreted accordingly.

## Usage

``` r
centrality_closeness(csr)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

## Value

numeric vector of per-node centrality scores.
