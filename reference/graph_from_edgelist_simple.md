# Build an igraph from an integer edgelist (1-based), undirected simple.

Build an igraph from an integer edgelist (1-based), undirected simple.

## Usage

``` r
graph_from_edgelist_simple(edges, n = NULL)
```

## Arguments

- edges:

  a two-column matrix (or coercible) of 1-based vertex pairs.

- n:

  optional vertex count; padded with isolates if larger than observed.

## Value

an undirected simple igraph object.
