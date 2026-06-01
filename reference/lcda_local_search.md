# First-improvement local search, with automatic VNMI dispatch for n \> 300.

First-improvement local search, with automatic VNMI dispatch for n \>
300.

## Usage

``` r
lcda_local_search(
  csr,
  construction,
  centrality = "eigen",
  n_threshold = 300,
  vnmi_n_prime = 300,
  vnmi_eps = 1e-04,
  verbose = FALSE
)
```

## Arguments

- csr:

  a CSR list (indptr, indices, igraph) from the internal converter.

- construction:

  a list as returned by \[lcda_construct()\]/\[lcda_repair()\].

- centrality:

  leader-selection centrality passed on to \[lcda_repair()\] so the
  re-validated leader matches the chosen measure.

- n_threshold:

  node count above which VNMI is used instead of full search.

- vnmi_n_prime:

  size of the sampled subset V' used by VNMI.

- vnmi_eps:

  minimum per-move modularity gain accepted by VNMI.

- verbose:

  logical; emit a cli trace.

## Value

the input \`construction\` with updated \`membership\`, \`Q\`, and
leaders.

## See also

\[as_csr()\], \[lcda_construct()\], \[lcda_repair()\].

## Examples

``` r
csr <- as_csr(igraph::make_graph("Zachary"))
sol <- lcda_local_search(csr, lcda_construct(csr, 0.1, 0.3))
sol$Q                          # modularity after local search
#> [1] 0.3717949
```
