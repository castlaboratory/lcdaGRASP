# Two-stage non-parametric comparison: Kruskal-Wallis omnibus followed by pairwise permutation tests with Bonferroni correction.

Mirrors the procedure described in section 5.1 of the paper.

## Usage

``` r
kruskal_then_permute(
  value,
  group,
  B_perm = 10000,
  alpha = 0.05,
  verbose = FALSE
)
```

## Arguments

- value:

  numeric vector of observations.

- group:

  factor or character vector of group labels.

- B_perm:

  number of permutations per pairwise test.

- alpha:

  global significance level.

- verbose:

  logical; emit a cli trace of the decisions.

## Value

list(omnibus, pairwise, letters, decision); \`pairwise\` is a tibble and
\`letters\` is a compact-letter display compatible with the paper's
tables.
