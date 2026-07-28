# R/metrics.R
#
# lcda_metrics(): the single metric surface behind the tables and figures of
# Ospina et al. (2026). Everything the paper reports about a solution --
# modularity Q, the NCE leader score H (global and community-conditioned),
# recovery against a ground truth (NMI / ARI / VI / split-join), community
# counts and size distribution, runtime, search diagnostics (trace dispersion,
# lexicographic decisiveness), and the per-community / per-leader breakdowns --
# is computed here from a result object plus the graph it was fitted on.
#
# Design notes:
#   * ONE entry point with a `level` argument, rather than a scatter of small
#     accessors, so that a user can go graph -> algorithm -> metrics -> figure.
#   * Kernels are reused, never re-implemented: Q comes from `modularity_cpp`
#     (via `modularity_score()`), H from `nce_score_cpp` / `nce_local_cpp`,
#     eigenvector centrality from `eigen_centrality_cpp`, and the partition
#     comparison indices from `igraph::compare()`.
#   * Weighted graphs: every structural sum below is taken over the CSR weight
#     vector, which is all-ones for an unweighted graph, so the weighted and
#     unweighted cases share one code path (and the "edge counts" become weight
#     sums on a weighted graph -- stated in the docs, per the Box lens).

.lcda_or <- function(a, b) if (is.null(a)) b else a

# --- normalise the many accepted inputs to (graph, membership, leaders) ------

.lcda_solution <- function(object, graph = NULL) {
  extras <- list()
  algo <- "unknown"
  leaders <- NULL
  leaders_source <- "given"

  if (inherits(object, c("lcda_grasp_result", "lcda_gr_result"))) {
    algo <- if (inherits(object, "lcda_gr_result")) "LCDA-GR" else "LCDA-GRASP"
    membership <- object$best$membership
    leaders <- object$best$leaders
    graph <- .lcda_or(graph, object$graph)
    extras <- list(
      Q_reported = object$best$Q, H_reported = object$best$H,
      best_iter = object$best$iter, elapsed = object$elapsed,
      trace_Q = object$trace_Q, trace_H = object$trace_H,
      params = object$params,
      lex_decisive = length(.lcda_or(object$lex_decisive_iters, integer(0))))
  } else if (inherits(object, "lcda_ecg_result")) {
    algo <- "LCDA-ECG"
    membership <- object$membership
    leaders <- object$leaders
    graph <- .lcda_or(graph, object$graph)
    extras <- list(
      Q_reported = object$Q, elapsed = object$elapsed, params = object$params,
      Q_consensus_weighted = object$Q_consensus_weighted,
      confidence = object$confidence, lead_count = object$lead_count,
      leaders_central = object$leaders_central, is_overlap = object$is_overlap)
  } else if (inherits(object, "communities")) {
    algo <- .lcda_or(object$algorithm, "igraph")
    membership <- as.integer(igraph::membership(object))
  } else if (is.list(object) && !is.null(object$membership)) {
    algo <- .lcda_or(object$algorithm, "custom")
    membership <- as.integer(object$membership)
    leaders <- object$leaders
  } else if (is.numeric(object) || is.factor(object)) {
    algo <- "custom"
    membership <- as.integer(as.factor(object))
  } else {
    cli::cli_abort(c(
      "Don't know how to read a solution out of a {.cls {class(object)[1]}}.",
      i = "Supply an {.cls lcda_grasp_result}, {.cls lcda_gr_result}, {.cls lcda_ecg_result}, an {.cls igraph} {.cls communities} object, a {.code list(membership=, leaders=)}, or a bare membership vector."))
  }

  if (is.null(graph))
    cli::cli_abort(c("{.arg graph} is required.",
                     i = "Only results produced by this package carry their own graph."))
  .check_graph(graph)
  csr <- .as_csr(graph)
  graph <- csr$igraph            # the simplified graph the kernels actually see

  membership <- as.integer(membership)
  if (length(membership) != csr$n)
    cli::cli_abort("{.arg membership} has length {length(membership)} but the graph has {csr$n} vertices.")
  if (anyNA(membership)) cli::cli_abort("{.arg membership} must not contain {.val NA}.")
  # dense, contiguous, 1-based community ids
  membership <- match(membership, sort(unique(membership)))

  if (is.null(leaders)) {
    # No leaders supplied (e.g. a Louvain/Leiden baseline): derive the paper's
    # two-stage leader, the top-eigenvector node of each community. Recorded in
    # the output so the reader knows these were not jointly optimised.
    ev <- eigen_centrality_cpp(csr$indptr, csr$indices)
    leaders <- vapply(sort(unique(membership)), function(cc) {
      idx <- which(membership == cc); idx[which.max(ev[idx])]
    }, integer(1))
    leaders_source <- "derived (top eigenvector per community)"
  }
  leaders <- as.integer(leaders)
  if (any(leaders < 1L | leaders > csr$n))
    cli::cli_abort("{.arg leaders} must be 1-based vertex indices in [1, {csr$n}].")

  list(graph = graph, csr = csr, membership = membership, leaders = leaders,
       leaders_source = leaders_source, algorithm = algo, extras = extras)
}

# --- edge-level scatter shared by the community and leader summaries ---------

.lcda_edge_frame <- function(csr) {
  n <- csr$n
  src <- rep.int(seq_len(n), diff(csr$indptr))
  list(src = src, dst = csr$indices + 1L, w = csr$weights)
}

# group-sum with a fixed number of groups (empty groups become 0)
.group_sum <- function(x, g, ngroups) {
  out <- numeric(ngroups)
  if (!length(x)) return(out)
  rs <- rowsum(as.numeric(x), group = g, reorder = TRUE)
  out[as.integer(rownames(rs))] <- rs[, 1L]
  out
}

# Participation coefficient P_v = 1 - sum_c (k_vc / k_v)^2: 0 when all of a
# node's edges stay inside one community, near 1 when they are spread evenly.
# This is the statistic the paper uses to characterise the bridge (overlap)
# nodes of the ensemble consensus.
.participation <- function(csr, membership) {
  n <- csr$n
  d <- max(membership)
  ef <- .lcda_edge_frame(csr)
  deg <- .group_sum(ef$w, ef$src, n)
  if (!length(ef$w)) return(rep(0, n))
  # one group per (node, neighbouring community) pair
  key <- (as.numeric(ef$src) - 1) * d + membership[ef$dst]
  rs <- rowsum(ef$w, group = key, reorder = FALSE)
  kk <- as.numeric(rownames(rs))
  node <- (kk - 1) %/% d + 1
  sq <- .group_sum(as.numeric(rs)^2, node, n)
  ifelse(deg > 0, 1 - sq / deg^2, 0)
}

.lcda_truth <- function(truth, n) {
  if (is.null(truth)) return(NULL)
  if (inherits(truth, "communities")) truth <- igraph::membership(truth)
  truth <- as.integer(as.factor(truth))
  if (length(truth) != n)
    cli::cli_abort("{.arg truth} has length {length(truth)} but the graph has {n} vertices.")
  truth
}

# --- the public surface ------------------------------------------------------

#' Paper-grade metrics for a community-and-leader solution
#'
#' Computes, from a fitted result (or any partition) plus the graph it was
#' fitted on, every quantity the tables and figures of Ospina et al. (2026)
#' report: modularity `Q`, the NCE leader score `H` in both its global and its
#' community-conditioned form, recovery against a ground truth (NMI, ARI, VI,
#' split-join), the number and size distribution of the communities, wall-clock
#' runtime, and the search diagnostics (trace dispersion, the best iteration,
#' and how often the lexicographic tie-break on `H` was actually decisive).
#'
#' This is the bridge between the algorithms and the reporting: a user can go
#' from a bare [igraph::igraph] to a publication table in three calls
#' (`lcda_gr()` -> `lcda_metrics()` -> `plot()`), without touching the
#' precomputed datasets in [lcda_data()].
#'
#' @section Levels:
#' \describe{
#'   \item{`"overall"`}{One row per metric, in tidy long form with columns
#'     `algorithm`, `scope`, `metric`, `value`. `scope` groups the metrics into
#'     `"partition"`, `"leaders"`, `"recovery"` (only when `truth` is given),
#'     `"search"`, and `"consensus"` (LCDA-ECG only). Long form is used because
#'     the metric set is heterogeneous and grows with the input; pivot with
#'     `tidyr::pivot_wider()` for a paper-style row.}
#'   \item{`"community"`}{One row per community: `size`, its `leader`, the
#'     internal edge weight and the cut (`boundary_edges`), internal density,
#'     conductance, and the community's additive contribution to `Q` (these sum
#'     exactly to `Q`). With `truth`, also the dominant ground-truth label and
#'     the community's purity.}
#'   \item{`"leader"`}{One row per leader: degree (total, within- and
#'     between-community), eigenvector centrality, the node-level NCE term, the
#'     participation coefficient, and the leader's rank and percentile by
#'     within-community degree, which is the form the paper's external leader
#'     validation takes. Supply `node_score` to rank the leaders by an external
#'     signal instead (see below).}
#' }
#'
#' @section Validating leaders against an external signal:
#' The paper validates its leaders by asking where each one sits in the
#' *within-community* distribution of an outside quantity (citations, in its
#' OpenAlex study). Pass that quantity as `node_score`, one value per vertex,
#' and the leader level gains `score`, `score_rank_in_community` and
#' `score_pctile_in_community`, while the overall level gains a `leader_score`
#' scope with the mean percentile and the top-1 / top-3 hit rates. A percentile
#' near 0.5 means the leaders are no better than chance.
#'
#' @section Weighted graphs:
#' All edge sums are taken over the graph's `weight` attribute when present, so
#' on a weighted graph the "edge" columns are weight sums and `degree` columns
#' are strengths. `Q` is weight-aware; the NCE leader score is deliberately
#' structural (unweighted), matching the algorithms.
#'
#' @param object a fitted [lcda_grasp()], [lcda_gr()] or [lcda_ecg()] result; an
#'   igraph `communities` object (e.g. from [igraph::cluster_louvain()]); a
#'   `list(membership =, leaders =)`; or a bare membership vector. When no
#'   leaders are available they are derived as the top-eigenvector node of each
#'   community, and the output records that.
#' @param graph the [igraph::igraph] the solution was computed on. Optional for
#'   results produced by this package, which carry their own (simplified) graph.
#' @param truth optional ground-truth community labels (a vector of length
#'   `vcount(graph)`, or a `communities` object). Enables the recovery metrics.
#' @param level one of `"overall"`, `"community"`, `"leader"`; see *Levels*.
#' @param node_score optional numeric vector of length `vcount(graph)` holding
#'   an external per-node signal against which to score the leaders; see
#'   *Validating leaders against an external signal*.
#'
#' @return A [tibble::tibble]. For `level = "overall"`, the long
#'   `algorithm`/`scope`/`metric`/`value` form; otherwise one row per community
#'   or per leader.
#'
#' @seealso [lcda_grasp()], [lcda_gr()], [lcda_ecg()],
#'   [plot.lcda_grasp_result()] for the matching figure, and [lcda_data()] for
#'   the precomputed panels of the paper.
#'
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_grasp(g, B = 20, seed = 1)
#'
#' m <- lcda_metrics(res)
#' subset(m, scope == "partition")
#'
#' # per-community and per-leader breakdowns
#' lcda_metrics(res, level = "community")
#' lcda_metrics(res, level = "leader")
#'
#' # recovery against a ground truth, and a baseline scored the same way
#' truth <- c(rep(1, 17), rep(2, 17))
#' subset(lcda_metrics(res, truth = truth), scope == "recovery")
#' subset(lcda_metrics(igraph::cluster_louvain(g), g, truth = truth),
#'        scope == "recovery")
#'
#' # do the leaders top their own community on an external signal?
#' external <- igraph::betweenness(g)
#' subset(lcda_metrics(res, node_score = external), scope == "leader_score")
#' @export
lcda_metrics <- function(object, graph = NULL, truth = NULL,
                         level = c("overall", "community", "leader"),
                         node_score = NULL) {
  level <- match.arg(level)
  s <- .lcda_solution(object, graph)
  truth <- .lcda_truth(truth, s$csr$n)
  if (!is.null(node_score)) {
    if (!is.numeric(node_score) || length(node_score) != s$csr$n)
      cli::cli_abort("{.arg node_score} must be a numeric vector of length {s$csr$n}.")
  }
  switch(level,
    overall   = .lcda_metrics_overall(s, truth, node_score),
    community = .lcda_metrics_community(s, truth),
    leader    = .lcda_metrics_leader(s, truth, node_score))
}

# --- level = "overall" -------------------------------------------------------

.lcda_metrics_overall <- function(s, truth, node_score = NULL) {
  g <- s$graph; csr <- s$csr; mem <- s$membership; lead <- s$leaders
  n <- csr$n
  d <- max(mem)
  sizes <- tabulate(mem, nbins = d)
  two_m <- sum(csr$weights)

  Q <- modularity_score(g, mem)
  H_global <- nce_score(g, lead)
  H_local  <- nce_local_score(g, mem, lead)

  ef <- .lcda_edge_frame(csr)
  deg <- .group_sum(ef$w, ef$src, n)                       # degree / strength
  # closed-neighbourhood coverage of the leader set
  covered <- unique(c(lead, ef$dst[ef$src %in% lead]))

  add <- function(scope, ...) {
    v <- c(...)
    tibble::tibble(algorithm = s$algorithm, scope = scope,
                   metric = names(v), value = as.numeric(v))
  }

  out <- list(
    add("partition",
        n_nodes = n, n_edges = igraph::ecount(g), total_edge_weight = two_m / 2,
        n_communities = d, Q = Q,
        size_min = min(sizes), size_median = stats::median(sizes),
        size_mean = mean(sizes), size_max = max(sizes),
        size_sd = if (d > 1) stats::sd(sizes) else 0,
        singleton_communities = sum(sizes == 1L)),
    add("leaders",
        n_leaders = length(lead), H_global = H_global, H_local = H_local,
        leader_degree_mean = mean(deg[lead]), leader_degree_max = max(deg[lead]),
        leader_coverage = length(covered) / n,
        leader_participation_mean = mean(.participation(csr, mem)[lead]),
        leaders_derived = as.numeric(s$leaders_source != "given"))
  )

  if (!is.null(node_score)) {
    lm <- .lcda_metrics_leader(s, NULL, node_score)
    out[[length(out) + 1L]] <- add("leader_score",
      score_pctile_mean = mean(lm$score_pctile_in_community),
      score_top1_rate = mean(lm$score_rank_in_community == 1L),
      score_top3_rate = mean(lm$score_rank_in_community <= 3L))
  }

  if (!is.null(truth)) {
    out[[length(out) + 1L]] <- add("recovery",
      nmi = igraph::compare(truth, mem, method = "nmi"),
      ari = igraph::compare(truth, mem, method = "adjusted.rand"),
      rand = igraph::compare(truth, mem, method = "rand"),
      vi = igraph::compare(truth, mem, method = "vi"),
      split_join = igraph::compare(truth, mem, method = "split.join"),
      n_communities_truth = length(unique(truth)))
  }

  ex <- s$extras
  search <- c(
    if (!is.null(ex$params$B)) c(B = ex$params$B),
    if (!is.null(ex$best_iter)) c(best_iter = ex$best_iter),
    if (!is.null(ex$elapsed)) c(elapsed_sec = ex$elapsed),
    if (!is.null(ex$trace_Q)) c(
      trace_Q_mean = mean(ex$trace_Q), trace_Q_sd = stats::sd(ex$trace_Q),
      trace_Q_min = min(ex$trace_Q), trace_Q_max = max(ex$trace_Q)),
    if (!is.null(ex$trace_H)) c(trace_H_mean = mean(ex$trace_H)),
    if (!is.null(ex$lex_decisive) && !is.null(ex$params$B)) c(
      lex_decisive_n = ex$lex_decisive,
      lex_decisive_pct = 100 * ex$lex_decisive / ex$params$B))
  if (length(search)) out[[length(out) + 1L]] <- add("search", search)

  consensus <- c(
    if (!is.null(ex$Q_consensus_weighted)) c(Q_consensus_weighted = ex$Q_consensus_weighted),
    if (!is.null(ex$confidence)) c(
      confidence_mean = mean(ex$confidence), confidence_min = min(ex$confidence)),
    if (!is.null(ex$is_overlap)) c(overlap_nodes = sum(ex$is_overlap),
                                   overlap_frac = mean(ex$is_overlap)),
    if (!is.null(ex$leaders_central)) c(
      leaders_agree_with_central = sum(sort(lead) == sort(ex$leaders_central))))
  if (length(consensus)) out[[length(out) + 1L]] <- add("consensus", consensus)

  do.call(rbind, out)
}

# --- level = "community" -----------------------------------------------------

.lcda_metrics_community <- function(s, truth) {
  csr <- s$csr; mem <- s$membership; lead <- s$leaders
  d <- max(mem)
  sizes <- tabulate(mem, nbins = d)
  two_m <- sum(csr$weights)

  ef <- .lcda_edge_frame(csr)
  cs <- mem[ef$src]; cd <- mem[ef$dst]
  vol <- .group_sum(ef$w, cs, d)                       # sum of degrees/strengths
  keep <- cs == cd
  internal2 <- .group_sum(ef$w[keep], cs[keep], d)     # 2 * internal weight
  cut <- vol - internal2
  denom <- pmin(vol, two_m - vol)

  leader_of <- rep(NA_integer_, d)
  leader_of[mem[lead]] <- lead
  deg <- .group_sum(ef$w, ef$src, csr$n)

  out <- tibble::tibble(
    algorithm      = s$algorithm,
    community      = seq_len(d),
    size           = sizes,
    leader         = leader_of,
    leader_name    = .lcda_vnames(s$graph, leader_of),
    leader_degree  = ifelse(is.na(leader_of), NA_real_, deg[leader_of]),
    # internal_edges sums to the graph's internal weight; boundary_edges is the
    # community's cut, so it counts every boundary edge once per side and its
    # sum is twice the total boundary weight.
    internal_edges = internal2 / 2,
    boundary_edges = cut,
    internal_density = ifelse(sizes > 1, (internal2 / 2) / (sizes * (sizes - 1) / 2), NA_real_),
    conductance    = ifelse(denom > 0, cut / denom, NA_real_),
    # additive decomposition of modularity: these sum exactly to Q
    q_contribution = internal2 / two_m - (vol / two_m)^2)

  if (!is.null(s$extras$confidence))
    out$mean_confidence <- vapply(seq_len(d),
      function(k) mean(s$extras$confidence[mem == k]), numeric(1))

  if (!is.null(truth)) {
    dom <- vapply(seq_len(d), function(k) {
      tb <- table(truth[mem == k]); as.integer(names(tb)[which.max(tb)])
    }, integer(1))
    pur <- vapply(seq_len(d), function(k) {
      tb <- table(truth[mem == k]); max(tb) / sum(tb)
    }, numeric(1))
    out$truth_dominant <- dom
    out$purity <- pur
  }
  out
}

.lcda_vnames <- function(g, v) {
  if (!("name" %in% igraph::vertex_attr_names(g))) return(as.character(v))
  nm <- igraph::vertex_attr(g, "name")
  ifelse(is.na(v), NA_character_, nm[v])
}

# --- level = "leader" --------------------------------------------------------

.lcda_metrics_leader <- function(s, truth, node_score = NULL) {
  csr <- s$csr; mem <- s$membership; lead <- s$leaders
  n <- csr$n
  ef <- .lcda_edge_frame(csr)
  deg <- .group_sum(ef$w, ef$src, n)
  keep <- mem[ef$src] == mem[ef$dst]
  deg_within <- .group_sum(ef$w[keep], ef$src[keep], n)
  ev <- eigen_centrality_cpp(csr$indptr, csr$indices)

  # where does the leader sit in its own community's distribution of a score?
  # (used for within-community degree, and for any external score supplied)
  rank_pct <- function(score) {
    r <- integer(length(lead)); p <- numeric(length(lead))
    for (i in seq_along(lead)) {
      peers <- score[mem == mem[lead[i]]]
      r[i] <- sum(peers > score[lead[i]]) + 1L
      p[i] <- if (length(peers) > 1) mean(peers <= score[lead[i]]) else 1
    }
    list(rank = r, pct = p)
  }
  rp <- rank_pct(deg_within)
  rank_in <- rp$rank; pct_in <- rp$pct

  out <- tibble::tibble(
    algorithm       = s$algorithm,
    community       = mem[lead],
    leader          = lead,
    leader_name     = .lcda_vnames(s$graph, lead),
    community_size  = tabulate(mem, nbins = max(mem))[mem[lead]],
    degree          = deg[lead],
    degree_within   = deg_within[lead],
    degree_between  = deg[lead] - deg_within[lead],
    eigen_centrality = ev[lead],
    nce_node        = vapply(lead, function(v)
                        nce_node_cpp(as.integer(round(igraph::degree(s$graph, v))), n),
                        numeric(1)),
    participation   = .participation(csr, mem)[lead],
    degree_rank_in_community = rank_in,
    degree_pctile_in_community = pct_in,
    source          = s$leaders_source)

  if (!is.null(node_score)) {
    rp2 <- rank_pct(node_score)
    out$score <- node_score[lead]
    out$score_rank_in_community <- rp2$rank
    out$score_pctile_in_community <- rp2$pct
  }

  ex <- s$extras
  if (!is.null(ex$lead_count)) {
    out$pool_lead_count <- ex$lead_count[lead]
    if (!is.null(ex$params$B)) out$pool_lead_freq <- ex$lead_count[lead] / ex$params$B
  }
  if (!is.null(ex$confidence)) out$confidence <- ex$confidence[lead]
  if (!is.null(ex$leaders_central))
    out$is_also_central <- lead %in% ex$leaders_central
  if (!is.null(truth)) out$truth_label <- truth[lead]
  out[order(out$community), ]
}
