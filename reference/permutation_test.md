# Two-sample permutation test for a difference in means. Robust to ties and to zero-variance pools (as discussed in section 5.1).

Two-sample permutation test for a difference in means. Robust to ties
and to zero-variance pools (as discussed in section 5.1).

## Usage

``` r
permutation_test(x, y, B = 10000, verbose = FALSE)
```

## Arguments

- x, y:

  numeric vectors of observations for the two groups.

- B:

  number of permutations.

- verbose:

  logical; emit a cli trace.

## Value

list with the observed mean difference and the permutation p-value.
