# Paper-grade metrics for a community-and-leader solution

Computes, from a fitted result (or any partition) plus the graph it was
fitted on, every quantity the tables and figures of Ospina et al. (2026)
report: modularity \`Q\`, the NCE leader score \`H\` in both its global
and its community-conditioned form, recovery against a ground truth
(NMI, ARI, VI, split-join), the number and size distribution of the
communities, wall-clock runtime, and the search diagnostics (trace
dispersion, the best iteration, and how often the lexicographic
tie-break on \`H\` was actually decisive).

## Usage

``` r
lcda_metrics(
  object,
  graph = NULL,
  truth = NULL,
  level = c("overall", "community", "leader"),
  node_score = NULL
)
```

## Arguments

- object:

  a fitted \[lcda_grasp()\], \[lcda_gr()\] or \[lcda_ecg()\] result; an
  igraph \`communities\` object (e.g. from
  \[igraph::cluster_louvain()\]); a \`list(membership =, leaders =)\`;
  or a bare membership vector. When no leaders are available they are
  derived as the top-eigenvector node of each community, and the output
  records that.

- graph:

  the \[igraph::igraph\] the solution was computed on. Optional for
  results produced by this package, which carry their own (simplified)
  graph.

- truth:

  optional ground-truth community labels (a vector of length
  \`vcount(graph)\`, or a \`communities\` object). Enables the recovery
  metrics. The labels you supply are preserved: \`truth_dominant\` and
  \`truth_label\` report \`"X"\`, not the integer code \`as.factor()\`
  gives it.

- level:

  one of \`"overall"\`, \`"community"\`, \`"leader"\`; see \*Levels\*.

- node_score:

  optional numeric vector of length \`vcount(graph)\` holding an
  external per-node signal against which to score the leaders; see
  \*Validating leaders against an external signal\*.

## Value

A \[tibble::tibble\]. For \`level = "overall"\`, the long
\`algorithm\`/\`scope\`/\`metric\`/\`value\` form; otherwise one row per
community or per leader.

## Details

This is the bridge between the algorithms and the reporting: a user can
go from a bare \[igraph::igraph\] to a publication table in three calls
(\`lcda_gr()\` -\> \`lcda_metrics()\` -\> \`plot()\`), without touching
the precomputed datasets in \[lcda_data()\].

## Levels

- \`"overall"\`:

  One row per metric, in tidy long form with columns \`algorithm\`,
  \`scope\`, \`metric\`, \`value\`. \`scope\` groups the metrics into
  \`"partition"\`, \`"leaders"\`, \`"recovery"\` (only when \`truth\` is
  given), \`"search"\`, and \`"consensus"\` (LCDA-ECG only). The two
  lexicographic metrics in \`"search"\` appear only for \[lcda_gr()\],
  the one algorithm that instruments the tie-break; they are absent,
  rather than reported as a measured-looking zero, for \[lcda_grasp()\].
  \`elapsed_sec\` times the search only: it starts after the graph has
  been converted by \[as_csr()\], so it is slightly smaller than timing
  the whole call from the outside (measured at 0.3-2 the metric set is
  heterogeneous and grows with the input; pivot with
  \`tidyr::pivot_wider()\` for a paper-style row.

- \`"community"\`:

  One row per community: \`size\`, its \`leader\`, the internal edge
  weight and the cut (\`boundary_edges\`), internal density,
  conductance, and the community's additive contribution to \`Q\` (these
  sum exactly to \`Q\`). With \`truth\`, also the dominant ground-truth
  label and the community's purity.

- \`"leader"\`:

  One row per leader: degree (total, within- and between-community),
  eigenvector centrality, the node-level NCE term, the participation
  coefficient, and the leader's rank and percentile by within-community
  degree, which is the form the paper's external leader validation
  takes. Supply \`node_score\` to rank the leaders by an external signal
  instead (see below).

## Validating leaders against an external signal

The paper validates its leaders by asking where each one sits in the
\*within-community\* distribution of an outside quantity (citations, in
its OpenAlex study). Pass that quantity as \`node_score\`, one value per
vertex, and the leader level gains \`score\`,
\`score_rank_in_community\` and \`score_pctile_in_community\`, while the
overall level gains a \`leader_score\` scope with the mean percentile
and the top-1 / top-3 hit rates. A percentile near 0.5 means the leaders
are no better than chance.

## Weighted graphs

All edge sums are taken over the graph's \`weight\` attribute when
present, so on a weighted graph the "edge" columns (\`internal_edges\`,
\`boundary_edges\`, \`total_edge_weight\`) are weight sums and the
\`degree\` columns are strengths. \`Q\` and \`conductance\` are
weight-aware and stay in their usual ranges. Two quantities are
deliberately \*\*structural\*\* (computed from edge counts, ignoring
weights), because a weight sum in their numerator would push them out of
their defining range: \`internal_density\`, which is a proportion of the
possible pairs and therefore always in \\0, 1\\, and the NCE leader
score, which matches how the algorithms themselves treat weights.

## See also

\[lcda_grasp()\], \[lcda_gr()\], \[lcda_ecg()\],
\[plot.lcda_grasp_result()\] for the matching figure, and
\[lcda_data()\] for the precomputed panels of the paper.

## Examples

``` r
g <- igraph::make_graph("Zachary")
res <- lcda_grasp(g, B = 20, seed = 1)

m <- lcda_metrics(res)
subset(m, scope == "partition")
#> # A tibble: 11 × 4
#>    algorithm  scope     metric                 value
#>    <chr>      <chr>     <chr>                  <dbl>
#>  1 LCDA-GRASP partition n_nodes               34    
#>  2 LCDA-GRASP partition n_edges               78    
#>  3 LCDA-GRASP partition total_edge_weight     78    
#>  4 LCDA-GRASP partition n_communities          3    
#>  5 LCDA-GRASP partition Q                      0.402
#>  6 LCDA-GRASP partition size_min               5    
#>  7 LCDA-GRASP partition size_median           12    
#>  8 LCDA-GRASP partition size_mean             11.3  
#>  9 LCDA-GRASP partition size_max              17    
#> 10 LCDA-GRASP partition size_sd                6.03 
#> 11 LCDA-GRASP partition singleton_communities  0    

# per-community and per-leader breakdowns
lcda_metrics(res, level = "community")
#> # A tibble: 3 × 11
#>   algorithm  community  size leader leader_name leader_degree internal_edges
#>   <chr>          <int> <int>  <int> <chr>               <dbl>          <dbl>
#> 1 LCDA-GRASP         1    17     34 34                     17             34
#> 2 LCDA-GRASP         2    12      1 1                      16             24
#> 3 LCDA-GRASP         3     5      6 6                       4              6
#> # ℹ 4 more variables: boundary_edges <dbl>, internal_density <dbl>,
#> #   conductance <dbl>, q_contribution <dbl>
lcda_metrics(res, level = "leader")
#> # A tibble: 3 × 14
#>   algorithm  community leader leader_name community_size degree degree_within
#>   <chr>          <int>  <int> <chr>                <int>  <dbl>         <dbl>
#> 1 LCDA-GRASP         1     34 34                      17     17            14
#> 2 LCDA-GRASP         2      1 1                       12     16            10
#> 3 LCDA-GRASP         3      6 6                        5      4             3
#> # ℹ 7 more variables: degree_between <dbl>, eigen_centrality <dbl>,
#> #   nce_node <dbl>, participation <dbl>, degree_rank_in_community <int>,
#> #   degree_pctile_in_community <dbl>, source <chr>

# recovery against a ground truth, and a baseline scored the same way
truth <- c(rep(1, 17), rep(2, 17))
subset(lcda_metrics(res, truth = truth), scope == "recovery")
#> # A tibble: 6 × 4
#>   algorithm  scope    metric               value
#>   <chr>      <chr>    <chr>                <dbl>
#> 1 LCDA-GRASP recovery nmi                  0.310
#> 2 LCDA-GRASP recovery ari                  0.289
#> 3 LCDA-GRASP recovery rand                 0.647
#> 4 LCDA-GRASP recovery vi                   1.17 
#> 5 LCDA-GRASP recovery split_join          17    
#> 6 LCDA-GRASP recovery n_communities_truth  2    
subset(lcda_metrics(igraph::cluster_louvain(g), g, truth = truth),
       scope == "recovery")
#> # A tibble: 6 × 4
#>   algorithm   scope    metric               value
#>   <chr>       <chr>    <chr>                <dbl>
#> 1 multi level recovery nmi                  0.277
#> 2 multi level recovery ari                  0.139
#> 3 multi level recovery rand                 0.576
#> 4 multi level recovery vi                   1.46 
#> 5 multi level recovery split_join          25    
#> 6 multi level recovery n_communities_truth  2    

# do the leaders top their own community on an external signal?
external <- igraph::betweenness(g)
subset(lcda_metrics(res, node_score = external), scope == "leader_score")
#> # A tibble: 3 × 4
#>   algorithm  scope        metric            value
#>   <chr>      <chr>        <chr>             <dbl>
#> 1 LCDA-GRASP leader_score score_pctile_mean     1
#> 2 LCDA-GRASP leader_score score_top1_rate       1
#> 3 LCDA-GRASP leader_score score_top3_rate       1
```
