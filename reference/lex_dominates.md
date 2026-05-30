# Lexicographic dominance under (Q, H).

A pure predicate (no logging): returns whether (Q_new, H_new)
lexicographically dominates (Q_old, H_old) under the (modularity, NCE)
objective with tie tolerance \`tol\`.

## Usage

``` r
lex_dominates(Q_new, H_new, Q_old, H_old, tol = 1e-12)
```

## Arguments

- Q_new, H_new:

  modularity and NCE of the candidate solution.

- Q_old, H_old:

  modularity and NCE of the incumbent solution.

- tol:

  numeric tolerance for treating two Q values as tied.

## Value

logical scalar.
