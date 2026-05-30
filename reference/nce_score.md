# Global NCE leader score (Eq. 14).

Global NCE leader score (Eq. 14).

## Usage

``` r
nce_score(g, leaders, verbose = FALSE)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- leaders:

  integer vector of 1-based leader vertex indices.

- verbose:

  logical; emit a cli trace.

## Value

scalar NCE score H(C).

## Examples

``` r
g <- igraph::make_graph("Zachary")
nce_score(g, leaders = c(1, 34))
#> [1] 0.9987513
```
