# Package index

## Algorithms

The LCDA variants proposed in the paper, plus the ensemble-consensus
variant.

- [`lcda_grasp()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_grasp.md)
  : LCDA-GRASP: fixed-parameter variant.
- [`lcda_gr()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_gr.md)
  : LCDA-GR: Reactive variant with self-tuning of (alpha_c, alpha_s).
  Implements Algorithm 4 of the paper.
- [`lcda_ecg()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_ecg.md)
  : LCDA-ECG: ensemble-consensus community and leader detection.

## Construction phases

Building blocks of a single GRASP iteration (construct, repair, local
search).

- [`lcda_construct()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_construct.md)
  : LCDA construction - variant 1 (centrality computed once) or 2
  (adaptive).
- [`lcda_repair()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_repair.md)
  : Repair: ensure every community has exactly one leader by recomputing
  centrality within the community and designating the top-scoring node.
- [`lcda_local_search()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_local_search.md)
  : First-improvement local search, with automatic VNMI dispatch for n
  \> 300.

## Quality scores

Modularity, the NCE leader score (global and community-conditioned), and
the lexicographic objective.

- [`modularity_score()`](https://castlaboratory.github.io/lcdaGRASP/reference/modularity_score.md)
  : Compute modularity for an arbitrary partition.
- [`nce_score()`](https://castlaboratory.github.io/lcdaGRASP/reference/nce_score.md)
  : Global NCE leader score (Eq. 14).
- [`nce_local_score()`](https://castlaboratory.github.io/lcdaGRASP/reference/nce_local_score.md)
  : Community-conditioned NCE (proposed alternative).
- [`lex_dominates()`](https://castlaboratory.github.io/lcdaGRASP/reference/lex_dominates.md)
  : Lexicographic dominance under (Q, H).

## Centrality measures

Leader-selection centralities.

- [`centrality_eigen()`](https://castlaboratory.github.io/lcdaGRASP/reference/centrality_eigen.md)
  : Eigenvector centrality (own implementation, O(m log n)).
- [`centrality_betweenness()`](https://castlaboratory.github.io/lcdaGRASP/reference/centrality_betweenness.md)
  : Betweenness centrality - delegates to igraph::betweenness for now.
  Tagged for a future native Rcpp implementation (Brandes 2001).
- [`centrality_closeness()`](https://castlaboratory.github.io/lcdaGRASP/reference/centrality_closeness.md)
  : Closeness centrality - delegates to igraph::closeness. On a
  disconnected graph, plain closeness is ill-defined across components
  (nodes in tiny components score spuriously high, distorting leader
  selection); we fall back to harmonic centrality, the robust
  generalisation that handles unreachable pairs (Boldi & Vigna 2014).
  Note the two are NOT identical even on a connected graph – closeness
  inverts the mean distance, harmonic averages the inverse distances –
  so on disconnected inputs (e.g. PolBlogs, and the
  frequently-disconnected induced subgraphs of variant 2) the selected
  leader may differ from a per-component closeness; this is a deliberate
  robustness choice, recorded here so benchmark results are interpreted
  accordingly.

## Statistical comparison

The non-parametric procedure used in the paper’s experiments.

- [`kruskal_then_permute()`](https://castlaboratory.github.io/lcdaGRASP/reference/kruskal_then_permute.md)
  : Two-stage non-parametric comparison: Kruskal-Wallis omnibus followed
  by pairwise permutation tests with Bonferroni correction.
- [`permutation_test()`](https://castlaboratory.github.io/lcdaGRASP/reference/permutation_test.md)
  : Two-sample permutation test for a difference in means. Robust to
  ties and to zero-variance pools (as discussed in section 5.1).

## Plotting

- [`plot_partition()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot_partition.md)
  : Plot a leader-community partition (leaders drawn larger).
- [`plot_grasp_trajectory()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot_grasp_trajectory.md)
  : Plot the Q trajectory across GRASP iterations, with running maximum.
- [`plot_reactive_pk()`](https://castlaboratory.github.io/lcdaGRASP/reference/plot_reactive_pk.md)
  : Plot the evolution of selection probabilities p_k across iterations.

## Graph constructors & data

- [`as_csr()`](https://castlaboratory.github.io/lcdaGRASP/reference/as_csr.md)
  : Convert an igraph graph to the CSR representation used by the
  kernels
- [`as_graph()`](https://castlaboratory.github.io/lcdaGRASP/reference/as_graph.md)
  : Coerce a variety of inputs to igraph.
- [`graph_from_edgelist_simple()`](https://castlaboratory.github.io/lcdaGRASP/reference/graph_from_edgelist_simple.md)
  : Build an igraph from an integer edgelist (1-based), undirected
  simple.
- [`lcda_data()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_data.md)
  : Precomputed simulation datasets
- [`lcda_provenance()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_provenance.md)
  : Provenance of a shipped dataset (versions, date, checksum)
