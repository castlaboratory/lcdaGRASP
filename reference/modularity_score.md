# Compute modularity for an arbitrary partition.

Compute modularity for an arbitrary partition.

## Usage

``` r
modularity_score(g, membership, verbose = FALSE)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- membership:

  integer vector of 1-based community ids, length \`vcount(g)\`.

- verbose:

  logical; emit a cli trace.

## Value

scalar modularity Q.

## Examples

``` r
g <- igraph::make_graph("Zachary")
cl <- igraph::cluster_louvain(g)
modularity_score(g, igraph::membership(cl))
#> [1] 0.395217
```
