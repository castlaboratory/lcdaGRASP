# Coerce a variety of inputs to igraph.

A pure constructor (no logging).

## Usage

``` r
as_graph(x)
```

## Arguments

- x:

  an igraph object, a two-column edgelist matrix, or an adjacency
  matrix.

## Value

an undirected igraph object.

## Examples

``` r
el <- matrix(c(1, 2, 2, 3, 3, 1), ncol = 2, byrow = TRUE)
as_graph(el)
#> IGRAPH 8212642 U--- 3 3 -- 
#> + edges from 8212642:
#> [1] 1--2 1--3 2--3
```
