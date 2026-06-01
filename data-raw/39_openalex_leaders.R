# data-raw/39_openalex_leaders.R
#
# Fase 5 anchor (LEADER validation against an EXTERNAL impact signal). The other
# real-network study (Coauthor-Physics) validates community recovery; this one
# validates the *leader* half of the contribution, which no labelled dataset
# carries directly. We build a co-authorship graph from OpenAlex on the paper's
# own field (topic T10064, "Complex Network Analysis Techniques") and attach each
# author's total citation count -- a signal the method never sees. We then ask
# whether the framework's designated community leaders are high-impact authors.
#
# For each detected (consensus) community we score four leader rules by the
# citation PERCENTILE of the chosen author within the community:
#   consensus  -- LCDA-ECG's consensus leader (most-frequently designated)
#   central    -- top-eigenvector author in the community (post-hoc baseline)
#   top-degree -- most-connected author (naive baseline)
#   random     -- expected ~0.5 (Monte Carlo)
# plus the top-1 / top-3 rate (leader is the most / among 3 most cited).
#
# -> inst/extdata/openalex_leaders.rds  (list: per_community, summary, meta)

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

set.seed(20260529)
FILTER  <- Sys.getenv("OA_FILTER", "primary_topic.id:T10064,publication_year:2015-2024")
MAXW    <- as.integer(Sys.getenv("OA_MAXW", "5000"))
MINW    <- as.integer(Sys.getenv("OA_MINW", "2"))
B_ECG   <- as.integer(Sys.getenv("OA_B_ECG", "64"))
MIN_CSZ <- as.integer(Sys.getenv("OA_MIN_CSZ", "5"))   # min community size to score

py <- .venv_python()
helper <- file.path(PKG_DIR, "data-raw", "openalex_coauthor.py")
ef <- file.path(tempdir(), "oa_edges.txt"); imf <- file.path(tempdir(), "oa_impact.txt")
cache <- file.path(tempdir(), "oa_works.json")
st <- if (is.na(py)) 1L else suppressWarnings(system2(py,
  c(helper, "--filter", FILTER, "--max-works", MAXW, "--min-works", MINW,
    "--edges", ef, "--impact", imf, "--cache", cache), stdout = FALSE, stderr = FALSE))

if (is.na(py) || st != 0 || !file.exists(ef) || !file.exists(imf)) {
  cli::cli_alert_warning("OpenAlex fetch failed; dataset not (re)generated.")
} else {
  el <- as.matrix(utils::read.table(ef))
  imp <- utils::read.table(imf); citations_all <- imp[[1]]; nworks_all <- imp[[2]]
  g <- igraph::graph_from_edgelist(el, directed = FALSE)
  if (igraph::vcount(g) < length(citations_all))
    g <- igraph::add_vertices(g, length(citations_all) - igraph::vcount(g))
  igraph::V(g)$cit <- citations_all; igraph::V(g)$nw <- nworks_all
  g <- igraph::simplify(g)
  cm <- igraph::components(g)
  g <- igraph::induced_subgraph(g, which(cm$membership == which.max(cm$csize)))
  cit <- igraph::V(g)$cit; deg <- igraph::degree(g)
  cli::cli_h1("OpenAlex T10064 coauthorship LCC: n={vcount(g)} m={ecount(g)} (authors with >= {MINW} works)")

  res <- lcda_ecg(g, B = B_ECG, centrality = "eigen", similarity = "hpi", seed = 1)
  mem <- res$membership
  # percentile of a node's citations within its own community (1 = most cited)
  pctl <- function(member_idx, node) mean(cit[member_idx] <= cit[node])

  rows <- list()
  for (c in sort(unique(mem))) {
    mi <- which(mem == c)
    if (length(mi) < MIN_CSZ) next
    cons <- res$leaders[res$membership[res$leaders] == c][1]          # consensus leader
    cent <- res$leaders_central[res$membership[res$leaders_central] == c][1]  # top-eigenvector
    topd <- mi[which.max(deg[mi])]                                    # highest degree
    rnd  <- mean(vapply(1:200, function(i) pctl(mi, sample(mi, 1)), numeric(1)))
    rows[[length(rows) + 1]] <- tibble::tibble(
      community = c, size = length(mi),
      cons_pctl = if (!is.na(cons)) pctl(mi, cons) else NA_real_,
      cent_pctl = if (!is.na(cent)) pctl(mi, cent) else NA_real_,
      topd_pctl = pctl(mi, topd), rand_pctl = rnd,
      cons_is_top1 = !is.na(cons) && cit[cons] == max(cit[mi]),
      cons_in_top3 = !is.na(cons) && cit[cons] >= sort(cit[mi], decreasing = TRUE)[min(3, length(mi))])
  }
  per_community <- dplyr::bind_rows(rows)
  summary <- tibble::tibble(
    n_communities = nrow(per_community),
    consensus = mean(per_community$cons_pctl, na.rm = TRUE),
    central   = mean(per_community$cent_pctl, na.rm = TRUE),
    top_degree= mean(per_community$topd_pctl, na.rm = TRUE),
    random    = mean(per_community$rand_pctl, na.rm = TRUE),
    top1_rate = mean(per_community$cons_is_top1, na.rm = TRUE),
    top3_rate = mean(per_community$cons_in_top3, na.rm = TRUE))
  cat("\n-- mean leader citation-percentile within community --\n")
  print(as.data.frame(summary), row.names = FALSE, digits = 3)

  meta <- make_meta(
    title = "OpenAlex co-authorship: do detected leaders match external impact (citations)?",
    description = paste("Leader validation against an external signal. Co-authorship graph from",
                        "OpenAlex topic T10064 (Complex Network Analysis Techniques); each author's",
                        "total citation count is attached but never used by the method. For each",
                        "detected community we report the citation percentile of LCDA-ECG's",
                        "consensus leader vs top-eigenvector / top-degree / random, and the",
                        "top-1/top-3 rates."),
    seed = 20260529, n = igraph::vcount(g), m_edges = igraph::ecount(g),
    oa_filter = FILTER, max_works = MAXW, min_works = MINW, b_ecg = B_ECG, min_comm_size = MIN_CSZ,
    network_source = "OpenAlex (api.openalex.org) works in topic T10064; co-authorship edges, citation impact per author",
    limitations = c(
      "Citations correlate with co-authorship centrality, so high leader impact is partly expected; the value is that citations are an EXTERNAL signal and the comparison to baselines.",
      "Communities are detected (no ground-truth field labels here); only the leader signal is externally validated.",
      "A single topic slice and a citation snapshot; results depend on the OpenAlex query date."))
  save_dataset("openalex_leaders", list(per_community = per_community, summary = summary), meta)
  cli::cli_alert_success("OpenAlex leader validation done.")
}
