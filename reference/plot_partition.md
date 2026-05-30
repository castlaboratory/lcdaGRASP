# Plot a leader-community partition (leaders drawn larger).

Plot a leader-community partition (leaders drawn larger).

## Usage

``` r
plot_partition(g, membership, leaders, ...)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- membership:

  integer vector of 1-based community ids.

- leaders:

  integer vector of 1-based leader vertex indices.

- ...:

  further arguments passed to \[igraph::plot.igraph()\].

## Value

invisibly, \`NULL\` (called for its plotting side effect).

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_grasp(g, B = 20, seed = 1)
plot_partition(g, res$best$membership, res$best$leaders)
```
