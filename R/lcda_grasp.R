# R/lcda_grasp.R
#
# Top-level orchestrators for LCDA-GRASP (fixed parameters) and LCDA-GR
# (Reactive). Heavy kernels live in src/*.cpp via Rcpp.
#
# Conventions (see CLAUDE.md): C++ kernels are 0-based; the R wrappers
# translate to and from 1-based indexing at the boundary. User-facing logging
# goes through cli; every substantive function takes a `verbose` flag.

# --- internal: convert an igraph object to the CSR representation -----

.as_csr <- function(g) {
  if (!inherits(g, "igraph")) cli::cli_abort("{.arg g} must be an {.cls igraph} object.")
  if (igraph::is_directed(g)) cli::cli_abort("Only undirected graphs are supported.")
  if (!igraph::is_simple(g)) g <- igraph::simplify(g)

  n <- igraph::vcount(g)
  # CSR via Matrix; A is symmetric for an undirected simple graph, so CSC == CSR.
  # If the graph carries a numeric "weight" edge attribute, the kernels use it
  # (weighted modularity / strength); otherwise weights are all 1 and the
  # kernels reduce exactly to the unweighted case.
  has_w <- "weight" %in% igraph::edge_attr_names(g)
  A <- igraph::as_adjacency_matrix(g, type = "both",
                                   attr = if (has_w) "weight" else NULL, sparse = TRUE)
  A <- methods::as(A, "CsparseMatrix")
  list(
    n       = n,
    indptr  = as.integer(A@p),   # length n+1, 0-based
    indices = as.integer(A@i),   # length 2m, 0-based
    weights = as.numeric(A@x),   # length 2m, edge weights (all 1 if unweighted)
    igraph  = g
  )
}

# --- centralities -----------------------------------------------------

#' Eigenvector centrality (own implementation, O(m log n)).
#' @param csr a CSR list (indptr, indices, igraph) from the internal converter.
#' @return numeric vector of per-node centrality scores.
#' @export
centrality_eigen <- function(csr) {
  eigen_centrality_cpp(csr$indptr, csr$indices)
}

#' Betweenness centrality - delegates to igraph::betweenness for now.
#' Tagged for a future native Rcpp implementation (Brandes 2001).
#' @param csr a CSR list (indptr, indices, igraph) from the internal converter.
#' @return numeric vector of per-node centrality scores.
#' @export
centrality_betweenness <- function(csr) {
  igraph::betweenness(csr$igraph, directed = FALSE, normalized = FALSE)
}

#' Closeness centrality - delegates to igraph::closeness.
#' On a disconnected graph, plain closeness is ill-defined across components
#' (nodes in tiny components score spuriously high, distorting leader selection);
#' we fall back to harmonic centrality, the robust generalisation that handles
#' unreachable pairs (Boldi & Vigna 2014). Note the two are NOT identical even
#' on a connected graph -- closeness inverts the mean distance, harmonic averages
#' the inverse distances -- so on disconnected inputs (e.g. PolBlogs, and the
#' frequently-disconnected induced subgraphs of variant 2) the selected leader
#' may differ from a per-component closeness; this is a deliberate robustness
#' choice, recorded here so benchmark results are interpreted accordingly.
#' @param csr a CSR list (indptr, indices, igraph) from the internal converter.
#' @return numeric vector of per-node centrality scores.
#' @export
centrality_closeness <- function(csr) {
  g <- csr$igraph
  if (igraph::is_connected(g)) {
    igraph::closeness(g, normalized = TRUE)
  } else {
    igraph::harmonic_centrality(g, normalized = TRUE)
  }
}

.dispatch_centrality <- function(name) {
  switch(name,
    eigen       = centrality_eigen,
    betweenness = centrality_betweenness,
    closeness   = centrality_closeness,
    cli::cli_abort("Unknown centrality: {.val {name}}.")
  )
}

.dispatch_similarity <- function(name) {
  switch(name,
    hpi     = similarity_hpi_cpp,
    dice    = similarity_dice_cpp,
    jaccard = similarity_jaccard_cpp,
    cli::cli_abort("Unknown similarity: {.val {name}}.")
  )
}

# --- construction phases (Algorithms 1 and 2 of the paper) -----------

#' LCDA construction - variant 1 (centrality computed once) or 2 (adaptive).
#'
#' @param csr CSR object from the internal converter.
#' @param alpha_c numeric in [0,1] - centrality RCL parameter.
#' @param alpha_s numeric in [0,1] - similarity RCL parameter.
#' @param variant 1 (static centrality) or 2 (recomputed each iteration).
#' @param centrality one of "eigen", "betweenness", "closeness".
#' @param similarity one of "hpi", "dice", "jaccard".
#' @param verbose logical; emit a cli trace of the construction.
#'
#' @return list(membership, leaders, d) - membership a 1-based vector, leaders
#'   a 1-based integer vector of leader indices, d the number of communities.
#' @export
lcda_construct <- function(csr,
                           alpha_c, alpha_s,
                           variant = 1,
                           centrality = "eigen",
                           similarity = "hpi",
                           verbose = FALSE) {
  if (!variant %in% c(1, 2)) cli::cli_abort("{.arg variant} must be 1 or 2, not {.val {variant}}.")
  n <- csr$n
  membership <- rep(NA_integer_, n)
  leaders <- integer(0)
  Q_remaining <- seq_len(n) - 1L     # 0-based remaining pool

  if (variant == 1) cent_global <- .dispatch_centrality(centrality)(csr)
  sim_fn <- .dispatch_similarity(similarity)

  comm_id <- 0L
  while (length(Q_remaining) > 0) {
    if (variant == 1) {
      cent <- cent_global[Q_remaining + 1L]
    } else {
      sub <- igraph::induced_subgraph(csr$igraph, Q_remaining + 1L)
      cent <- .dispatch_centrality(centrality)(.as_csr(sub))
    }

    # Closeness/betweenness on a disconnected induced subgraph (common in
    # variant 2) can return NaN/Inf; treat non-finite scores as the lowest
    # finite value so they never win the leader RCL but never empty it either.
    if (any(!is.finite(cent))) {
      fill <- if (any(is.finite(cent))) min(cent[is.finite(cent)]) else 0
      cent[!is.finite(cent)] <- fill
    }

    rcl_L_idx <- rcl_max_cpp(cent, alpha_c)            # 0-based into Q_remaining
    if (length(rcl_L_idx) == 0) cli::cli_abort("Empty leader RCL (check {.arg alpha_c}).")
    pick <- rcl_L_idx[sample.int(length(rcl_L_idx), 1)]
    L <- Q_remaining[pick + 1L]                        # 0-based node id
    leaders <- c(leaders, L + 1L)                      # store 1-based
    membership[L + 1L] <- comm_id
    Q_remaining <- setdiff(Q_remaining, L)

    if (length(Q_remaining) == 0) break

    sim_vals <- sim_fn(csr$indptr, csr$indices, L, as.integer(Q_remaining))
    rcl_S_idx <- rcl_min_cpp(sim_vals, alpha_s)        # 0-based into Q_remaining
    members <- Q_remaining[rcl_S_idx + 1L]
    membership[members + 1L] <- comm_id
    Q_remaining <- setdiff(Q_remaining, members)
    if (verbose) cli::cli_alert_info(
      "community {comm_id + 1L}: leader node {L + 1L}, +{length(members)} member{?s}")
    comm_id <- comm_id + 1L
  }

  # comm_id is not incremented for a final singleton-leader community (the loop
  # breaks before the increment); derive d from the membership instead.
  list(
    membership = membership + 1L,     # 1-based
    leaders    = leaders,
    d          = max(membership) + 1L # number of 0-based communities
  )
}

# --- repair step -----------------------------------------------------

#' Repair: ensure every community has exactly one leader by recomputing
#' centrality within the community and designating the top-scoring node.
#' @param csr a CSR list (indptr, indices, igraph) from the internal converter.
#' @param construction a list as returned by [lcda_construct()].
#' @param centrality leader-selection centrality (`"eigen"`, `"betweenness"`,
#'   or `"closeness"`); should match the one used in construction so the final
#'   leader is consistent with the chosen measure.
#' @param verbose logical; emit a cli trace.
#' @return the input `construction` with its `leaders` field rebuilt.
#' @export
lcda_repair <- function(csr, construction, centrality = "eigen", verbose = FALSE) {
  mem <- construction$membership          # 1-based community ids
  cent_fn <- .dispatch_centrality(centrality)
  # Iterate the community ids actually present (iterating 0..d-1 would mix
  # 0-based ids against the 1-based `mem` and hit empty communities).
  comms <- sort(unique(mem))
  new_leaders <- integer(0)
  for (cc in comms) {
    members_1based <- which(mem == cc)
    if (length(members_1based) == 1L) {
      new_leaders <- c(new_leaders, members_1based)
      next
    }
    sub <- igraph::induced_subgraph(csr$igraph, members_1based)
    cent_local <- cent_fn(.as_csr(sub))
    # Edgeless community => degenerate centrality; fall back to the first member.
    idx <- if (length(cent_local) == 0L || !any(is.finite(cent_local))) 1L else which.max(cent_local)
    new_leaders <- c(new_leaders, members_1based[idx])
  }
  if (verbose) cli::cli_alert_info("repaired {length(new_leaders)} leader{?s}")
  construction$leaders <- new_leaders
  construction
}

# --- local search ----------------------------------------------------

#' First-improvement local search, with automatic VNMI dispatch for n > 300.
#' @param csr a CSR list (indptr, indices, igraph) from the internal converter.
#' @param construction a list as returned by [lcda_construct()]/[lcda_repair()].
#' @param centrality leader-selection centrality passed on to [lcda_repair()]
#'   so the re-validated leader matches the chosen measure.
#' @param n_threshold node count above which VNMI is used instead of full search.
#' @param vnmi_n_prime size of the sampled subset V' used by VNMI.
#' @param vnmi_eps minimum per-move modularity gain accepted by VNMI.
#' @param verbose logical; emit a cli trace.
#' @return the input `construction` with updated `membership`, `Q`, and leaders.
#' @export
lcda_local_search <- function(csr, construction, centrality = "eigen",
                              n_threshold = 300,
                              vnmi_n_prime = 300, vnmi_eps = 1e-4,
                              verbose = FALSE) {
  mem0 <- as.integer(construction$membership - 1L)   # 0-based for C++
  d    <- as.integer(construction$d)
  use_vnmi <- csr$n > n_threshold
  if (use_vnmi) {
    res <- vnmi_local_search_cpp(csr$indptr, csr$indices, mem0, d, csr$weights,
                                 n_prime = vnmi_n_prime, eps = vnmi_eps)
  } else {
    res <- local_search_cpp(csr$indptr, csr$indices, mem0, d, csr$weights)
  }
  construction$membership <- as.integer(res$membership) + 1L
  construction$Q          <- res$modularity
  construction <- lcda_repair(csr, construction, centrality = centrality)  # leaders re-validated
  if (verbose) cli::cli_alert_info(
    "local search ({if (use_vnmi) 'VNMI' else 'first-improvement'}): Q = {round(res$modularity, 6)}")
  construction
}

# --- quality scores --------------------------------------------------

#' Compute modularity for an arbitrary partition.
#' @param g an igraph object (undirected, simple).
#' @param membership integer vector of 1-based community ids, length `vcount(g)`.
#' @param verbose logical; emit a cli trace.
#' @return scalar modularity Q.
#' @examples
#' g <- igraph::make_graph("Zachary")
#' cl <- igraph::cluster_louvain(g)
#' modularity_score(g, igraph::membership(cl))
#' @export
modularity_score <- function(g, membership, verbose = FALSE) {
  csr <- .as_csr(g)
  q <- modularity_cpp(csr$indptr, csr$indices, as.integer(membership - 1L), csr$weights)
  if (verbose) cli::cli_alert_info("Q = {round(q, 6)} over {length(unique(membership))} communities")
  q
}

#' Global NCE leader score (Eq. 14).
#' @param g an igraph object (undirected, simple).
#' @param leaders integer vector of 1-based leader vertex indices.
#' @param verbose logical; emit a cli trace.
#' @return scalar NCE score H(C).
#' @examples
#' g <- igraph::make_graph("Zachary")
#' nce_score(g, leaders = c(1, 34))
#' @export
nce_score <- function(g, leaders, verbose = FALSE) {
  deg <- igraph::degree(g, v = leaders)
  h <- nce_score_cpp(as.integer(deg), igraph::vcount(g))
  if (verbose) cli::cli_alert_info("H (global) = {round(h, 4)} over {length(leaders)} leaders")
  h
}

#' Community-conditioned NCE (proposed alternative).
#' @param g an igraph object (undirected, simple).
#' @param membership integer vector of 1-based community ids, length `vcount(g)`.
#' @param leaders integer vector of 1-based leader vertex indices.
#' @param verbose logical; emit a cli trace.
#' @return scalar community-conditioned NCE score.
#' @export
nce_local_score <- function(g, membership, leaders, verbose = FALSE) {
  csr <- .as_csr(g)
  h <- nce_local_cpp(csr$indptr, csr$indices,
                     as.integer(membership - 1L), as.integer(leaders - 1L))
  if (verbose) cli::cli_alert_info("H (community-conditioned) = {round(h, 4)}")
  h
}

#' Lexicographic dominance under (Q, H).
#'
#' A pure predicate (no logging): returns whether (Q_new, H_new)
#' lexicographically dominates (Q_old, H_old) under the (modularity, NCE)
#' objective with tie tolerance `tol`.
#' @param Q_new,H_new modularity and NCE of the candidate solution.
#' @param Q_old,H_old modularity and NCE of the incumbent solution.
#' @param tol numeric tolerance for treating two Q values as tied.
#' @return logical scalar.
#' @export
lex_dominates <- function(Q_new, H_new, Q_old, H_old, tol = 1e-12) {
  if (Q_new > Q_old + tol) return(TRUE)
  if (abs(Q_new - Q_old) <= tol && H_new > H_old + tol) return(TRUE)
  FALSE
}

# --- LCDA-GRASP outer loop ------------------------------------------

#' LCDA-GRASP: fixed-parameter variant.
#'
#' @param g an igraph object (undirected, simple).
#' @param alpha_c,alpha_s RCL parameters.
#' @param variant 1 or 2.
#' @param B number of GRASP iterations.
#' @param centrality,similarity metric names.
#' @param verbose logical; show a cli progress bar and a final summary.
#' @param seed integer RNG seed, or `NA` to leave the RNG untouched.
#' @return an object of class `lcda_grasp_result`: best partition, Q/H traces,
#'   and the parameters used.
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_grasp(g, B = 20, seed = 1)
#' res$best$Q
#' res$best$leaders
#' @export
lcda_grasp <- function(g,
                       alpha_c = 0.1, alpha_s = 0.3,
                       variant = 1, B = 50,
                       centrality = "eigen", similarity = "hpi",
                       verbose = FALSE,
                       seed = NA_integer_) {
  if (!is.na(seed)) set.seed(seed)
  csr <- .as_csr(g)
  g <- csr$igraph   # use the simplified graph so H (degree/n) matches the CSR
  best <- list(Q = -Inf, H = -Inf, membership = NULL, leaders = NULL, iter = NA)
  trace_Q <- numeric(B)
  trace_H <- numeric(B)

  if (verbose) cli::cli_progress_bar("LCDA-GRASP", total = B, clear = FALSE)
  for (it in seq_len(B)) {
    sol <- lcda_construct(csr, alpha_c, alpha_s, variant = variant,
                          centrality = centrality, similarity = similarity)
    sol <- lcda_repair(csr, sol, centrality = centrality)
    sol <- lcda_local_search(csr, sol, centrality = centrality)
    Q <- sol$Q
    H <- nce_score_cpp(as.integer(igraph::degree(g, v = sol$leaders)), igraph::vcount(g))
    trace_Q[it] <- Q
    trace_H[it] <- H
    if (lex_dominates(Q, H, best$Q, best$H)) {
      best <- list(Q = Q, H = H, membership = sol$membership,
                   leaders = sol$leaders, iter = it)
    }
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) {
    cli::cli_progress_done()
    cli::cli_alert_success(
      "best Q = {round(best$Q, 6)} (H = {round(best$H, 4)}) at iteration {best$iter}, {length(best$leaders)} communities")
  }
  structure(
    list(best = best, trace_Q = trace_Q, trace_H = trace_H,
         params = list(alpha_c = alpha_c, alpha_s = alpha_s, variant = variant,
                       B = B, centrality = centrality, similarity = similarity)),
    class = "lcda_grasp_result")
}

# --- LCDA-GR (Reactive) outer loop ----------------------------------

#' LCDA-GR: Reactive variant with self-tuning of (alpha_c, alpha_s).
#' Implements Algorithm 4 of the paper.
#'
#' @param g an igraph object (undirected, simple).
#' @param variant construction variant, 1 or 2.
#' @param B number of GRASP iterations.
#' @param centrality,similarity metric names.
#' @param m pool size (default 20, per paper).
#' @param y refresh period (default 3m).
#' @param alpha_c_range,alpha_s_range bounds for the uniform initial pool.
#' @param p_floor minimum probability mass reserved across the pool at each
#'   refresh, spread uniformly so every pair keeps `p_k >= p_floor/m > 0`. This
#'   prevents a pair that happened not to be sampled in a block from being
#'   permanently excluded (and matches the `p_k >= delta > 0` premise of
#'   Proposition 6). Set to 0 to recover the raw proportional rule.
#' @param verbose logical; show a cli progress bar and a final summary.
#' @param seed integer RNG seed, or `NA` to leave the RNG untouched.
#' @return an object of class `lcda_gr_result`: best partition, traces, the
#'   reactive pool state, and the H-decisive iterations.
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_gr(g, B = 30, seed = 1)
#' res$best$Q
#' length(res$lex_decisive_iters)  # how often H broke a Q-tie
#' @export
lcda_gr <- function(g,
                    variant = 1, B = 150,
                    centrality = "eigen", similarity = "hpi",
                    m = 20, y = NULL,
                    alpha_c_range = c(0.1, 0.9),
                    alpha_s_range = c(0.1, 0.5),
                    p_floor = 0.05,
                    verbose = FALSE,
                    seed = NA_integer_) {
  if (!is.na(seed)) set.seed(seed)
  if (is.null(y)) y <- 3L * m
  csr <- .as_csr(g)
  g <- csr$igraph   # use the simplified graph so H (degree/n) matches the CSR

  Psi <- tibble::tibble(
    alpha_c = stats::runif(m, alpha_c_range[1], alpha_c_range[2]),
    alpha_s = stats::runif(m, alpha_s_range[1], alpha_s_range[2])
  )
  p_k  <- rep(1 / m, m)
  mu_k <- rep(0.0, m)
  n_k  <- rep(0L, m)
  Q_star <- -Inf
  best <- list(Q = -Inf, H = -Inf, membership = NULL, leaders = NULL, iter = NA, k = NA)
  trace_Q <- numeric(B); trace_H <- numeric(B); trace_k <- integer(B)
  pk_history <- matrix(NA_real_, nrow = B, ncol = m)
  lex_decisive <- integer(0)   # iters where H broke a Q-tie

  if (verbose) cli::cli_progress_bar("LCDA-GR", total = B, clear = FALSE)
  for (it in seq_len(B)) {
    k <- sample.int(m, 1, prob = p_k)
    a_c <- Psi$alpha_c[k]; a_s <- Psi$alpha_s[k]

    sol <- lcda_construct(csr, a_c, a_s, variant = variant,
                          centrality = centrality, similarity = similarity)
    sol <- lcda_repair(csr, sol, centrality = centrality)
    sol <- lcda_local_search(csr, sol, centrality = centrality)
    Q <- sol$Q
    H <- nce_score_cpp(as.integer(igraph::degree(g, v = sol$leaders)), igraph::vcount(g))
    trace_Q[it] <- Q; trace_H[it] <- H; trace_k[it] <- k

    Q_tie <- abs(Q - best$Q) < 1e-12 && Q > -Inf
    if (lex_dominates(Q, H, best$Q, best$H)) {
      if (Q_tie) lex_decisive <- c(lex_decisive, it)
      best <- list(Q = Q, H = H, membership = sol$membership,
                   leaders = sol$leaders, iter = it, k = k)
    }
    if (Q > Q_star) Q_star <- Q

    n_k[k] <- n_k[k] + 1L
    mu_k[k] <- mu_k[k] + (Q - mu_k[k]) / n_k[k]   # Welford mean

    if (it %% y == 0 && Q_star > 0) {
      qk <- pmax(mu_k, 0) / Q_star
      p_k <- if (sum(qk) > 0) qk / sum(qk) else rep(1 / m, m)
      # Epsilon floor: keep every pair reachable (p_k >= p_floor/m > 0) so a
      # pair not sampled in a block is not excluded forever; also restores the
      # p_k >= delta > 0 premise of Proposition 6.
      if (p_floor > 0) p_k <- (1 - p_floor) * p_k + p_floor / m
    }
    pk_history[it, ] <- p_k
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) {
    cli::cli_progress_done()
    cli::cli_alert_success(
      "best Q = {round(best$Q, 6)} (H = {round(best$H, 4)}); H broke a Q-tie in {length(lex_decisive)}/{B} iterations")
  }

  structure(
    list(best = best, trace_Q = trace_Q, trace_H = trace_H, trace_k = trace_k,
         pk_history = pk_history, Psi = Psi, mu_k = mu_k, n_k = n_k,
         lex_decisive_iters = lex_decisive,
         params = list(B = B, m = m, y = y, variant = variant,
                       centrality = centrality, similarity = similarity)),
    class = "lcda_gr_result")
}

#' @export
print.lcda_grasp_result <- function(x, ...) {
  cli::cli_h1("LCDA-GRASP result")
  cli::cli_dl(c(
    "variant"        = "{x$params$variant}",
    "B"              = "{x$params$B}",
    "centrality / similarity" = "{x$params$centrality} / {x$params$similarity}",
    "best Q"         = "{format(x$best$Q, digits = 6)}",
    "best H"         = "{format(x$best$H, digits = 4)}",
    "best iteration" = "{x$best$iter}",
    "communities"    = "{length(x$best$leaders)}"
  ))
  invisible(x)
}

#' @export
print.lcda_gr_result <- function(x, ...) {
  pct <- 100 * length(x$lex_decisive_iters) / x$params$B
  cli::cli_h1("LCDA-GR (Reactive) result")
  cli::cli_dl(c(
    "variant"          = "{x$params$variant}",
    "B / m / y"        = "{x$params$B} / {x$params$m} / {x$params$y}",
    "best Q"           = "{format(x$best$Q, digits = 6)}",
    "best H"           = "{format(x$best$H, digits = 4)}",
    "best iter / pair" = "{x$best$iter} / {x$best$k}",
    "H-decisive updates" = "{length(x$lex_decisive_iters)} ({sprintf('%.1f%%', pct)} of B)"
  ))
  invisible(x)
}
