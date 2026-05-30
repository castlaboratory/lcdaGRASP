---
title: 'lcdaGRASP: Joint community and leader detection in complex networks via GRASP, Reactive GRASP, and ensemble consensus'
tags:
  - R
  - C++
  - network science
  - community detection
  - leader detection
  - modularity
  - GRASP metaheuristic
authors:
  - name: Raydonal Ospina
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Geiza Silva
    affiliation: 1
  - name: Francisco Jucelino Matos Junior
    affiliation: 1
  - name: André Leite
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Luiz Satoru Ochi
    affiliation: 2
affiliations:
  - name: Department of Statistics, Universidade Federal de Pernambuco (UFPE), Brazil
    index: 1
  - name: Institute of Computing, Universidade Federal Fluminense (UFF), Brazil
    index: 2
date: 29 May 2026
bibliography: paper.bib
---

# Summary

`lcdaGRASP` is an R package for the **joint** detection of communities **and
their leaders** in undirected networks. Most tools answer only one of these
questions: community-detection methods such as Louvain [@blondel2008louvain] and
Leiden [@traag2019leiden] partition the graph but designate no representative
node, while centrality measures rank influential nodes without delimiting the
communities they lead. `lcdaGRASP` produces both at once: a partition together
with one structurally grounded leader per community, where every assignment is
traceable to a named centrality or similarity score rather than to an opaque
embedding.

The package implements four algorithms from @ospina2026grasp — two
fixed-parameter variants (LCDA-GRASP) and two self-calibrating Reactive variants
(LCDA-GR) — and adds an ensemble-consensus variant, **LCDA-ECG**, that
aggregates the diverse pool of partitions a multi-start search produces
(otherwise discarded) into a single high-quality partition with a
consensus-derived leader, overlapping-community memberships, and a per-node
confidence map. Performance-critical kernels (incremental modularity and its
$\Delta Q$, community affinities, a Variable-Neighbourhood Multi-Improvement
local search, and power-iteration eigenvector centrality) are written in C++ via
`Rcpp` [@eddelbuettel2011rcpp]; the modularity machinery supports weighted graphs.
A companion `pkgdown` site, nine vignettes, precomputed reproduction datasets,
and a continuous-integration suite (multi-platform `R CMD check`, code coverage,
linting) accompany the package.

# Statement of need

Identifying *who leads which community* is central to applications in social,
biological, and information networks — influence analysis, targeted
intervention, and the interpretation of group structure. Two practical gaps
motivate `lcdaGRASP`:

1. **No unified, reproducible tool for joint leader-and-community detection.**
   Classical community detection (in `igraph` [@csardi2006igraph], `leidenAlg`,
   `linkcomm`) returns partitions without leaders; leader-driven methods
   [@kanawati2014licod; @ahajjam2018new; @sun2020leaderaware; @akachar2025leadcd]
   are described in papers but lack maintained, tested, open implementations.
   `lcdaGRASP` provides one, with unit-tested invariants and reproducible
   benchmarks.

2. **Determinism and interpretability.** Deterministic leader-driven heuristics
   return a single partition and cannot escape weak local optima; GNN-based
   methods require node attributes and designate leaders only implicitly through
   embedding geometry [@kipf2016variational], which is hard to justify in
   high-stakes settings. `lcdaGRASP` couples GRASP diversification with a
   Reactive parameter-tuning mechanism that needs no grid search, and its
   ensemble-consensus variant recovers planted structure on a par with the
   strongest modularity optimisers — and better at high mixing and large scale —
   while keeping every leader designation interpretable.

The package targets researchers and practitioners in network science who need
**leaders and communities together**, reproducibly, on attribute-free graphs of
up to tens of thousands of nodes. It also serves as the reference implementation
and reproducibility artifact for @ospina2026grasp.

# Functionality

- `lcda_grasp()`, `lcda_gr()` — fixed-parameter and Reactive GRASP algorithms.
- `lcda_ecg()` — ensemble consensus over the GRASP pool, returning the consensus
  partition, a consensus-derived leader per community, overlapping communities
  (`overlap = TRUE`) with bridge nodes, and a node-confidence map.
- Quality scores: modularity, the Node-Connection Entropy leader score
  (global and community-conditioned), and the lexicographic objective.
- Centralities, weighted/structural similarities, statistical comparison
  (Kruskal–Wallis + permutation), plotting, and cached reproduction datasets via
  `lcda_data()`.

# Acknowledgements

This work was developed at the CAST Laboratory, UFPE. [Funding to be completed.]

# References
