# lcdaGRASP — companion code for Ospina et al. (2026)

Companion code and critical review for the preprint:

> Ospina, R., Silva, G., Matos Junior, F. J., Leite, A., & Ochi, L. S.
> (2026). *A GRASP Framework for Community and Leader Detection in
> Complex Networks.* Preprint submitted to *Knowledge-Based Systems*.

The paper proposes four algorithms (LCDA-GRASP 1/2 and LCDA-GR 1/2) for
**joint community-and-leader detection**. This repository contains:

1.  **`lcdaGRASP/`** — a working R package with the four algorithms and
    their performance-critical kernels in C++ via Rcpp. Tidyverse style,
    `cli` logging, a `verbose` flag on every substantive function.
2.  **`lcdaGRASP/vignettes/`** — articles that render the
    publication-grade simulation panel (benchmarks, a **George-Box
    design-of-experiments** study, robustness, convergence, degeneracy,
    EDA, leader scores) from precomputed data, so they build without
    re-running any simulation.
3.  **`lcdaGRASP/data-raw/`** — the scripts that generate that data into
    `lcdaGRASP/inst/extdata/*.rds` (run once; see `run_all_data.R`).
4.  **`scripts/`** — the original standalone simulation drivers.
5.  **`artigo/`** — the LaTeX source of the manuscript, with errata
    fixes and literature updates applied this revision (see
    `artigo/CHANGES.md`).
6.  **`analise_critica.md`** — the critical review (four reasoning
    lenses: Wald, von Neumann, Tukey, Box).

Package website (pkgdown):
**<https://castlaboratory.github.io/lcdaGRASP/>**

## Install

``` r

# from the repository root
remotes::install_local("lcdaGRASP", dependencies = TRUE)
```

## Quick start

``` r

library(lcdaGRASP)
library(igraph)

g   <- make_graph("Zachary")
res <- lcda_grasp(g, alpha_c = 0.1, alpha_s = 0.3, B = 50, seed = 1, verbose = TRUE)
print(res)

res_r <- lcda_gr(g, variant = 1, B = 150, seed = 1)   # Reactive, self-tuning
print(res_r)
```

Precomputed simulation results back the articles and are listed with
[`lcda_data()`](https://castlaboratory.github.io/lcdaGRASP/reference/lcda_data.md);
load one with `lcda_data("repro_summary")`.

## Articles (by reasoning lens)

| Article | Lens | What it shows |
|----|----|----|
| Getting started | — | API tour |
| Benchmark reproduction | Wald | Paper’s numbers + erratum E3 (Louvain on Karate = 0.4198) |
| **DoE parameter tuning** | Box | Factorial screening → response surface (CCD) → canonical optimum for `(α_c, α_s)` |
| Robustness (LFR) | Wald | Q **and** ARI across the mixing sweep |
| Reactive convergence | von Neumann | Proposition 6 empirically (entropy decay) |
| Pool sensitivity | von Neumann | `(m, y, B)` grid + lex-decisive audit |
| Modularity degeneracy | Box | Near-best partition counts + NMI |
| EDA distributions | Tukey | Boxplots / ECDFs / (Q, H) scatter |
| Leader score | Box | Global vs community-conditioned NCE |

## Rebuild the data and the site

``` bash
Rscript lcdaGRASP/data-raw/run_all_data.R       # regenerate inst/extdata/*.rds (slow, once)
Rscript -e 'pkgdown::build_site("lcdaGRASP")'   # build docs/
```

## Errata flagged to the authors (see `artigo/CHANGES.md`)

- **E1** Table 9 duplicates LCDA-GRASP and LCDA-GR timings (impossible
  with different `B`).
- **E2** Lemma proof: smallest positive ΔQ is `1/(2m²)`, so `T ≤ 2m²`
  (was `2m`).
- **E3** Table 10 Louvain on Karate `0.3715`; current igraph
  Louvain/Leiden reach `0.4198` over 30 restarts — the `+12.9%` claim is
  baseline-dependent.
- **E4** Proposition 6 proof: Borel–Cantelli step made explicit (Lévy
  extension).
- **E5** Abstract framing softened.
- **E6** LeaderRank values look anomalous; verify the source.

## License

MIT. See `LICENSE`.
