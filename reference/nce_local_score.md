# Community-conditioned NCE (proposed alternative).

Community-conditioned NCE (proposed alternative).

## Usage

``` r
nce_local_score(g, membership, leaders, verbose = FALSE)
```

## Arguments

- g:

  an igraph object (undirected, simple).

- membership:

  integer vector of 1-based community ids, length \`vcount(g)\`.

- leaders:

  integer vector of 1-based leader vertex indices.

- verbose:

  logical; emit a cli trace.

## Value

scalar community-conditioned NCE score.
