# R/lcda_ecg.R
#
# LCDA-ECG: ensemble-consensus community-and-leader detection. The GRASP pool
# (the diverse partitions produced across iterations) is normally discarded
# except for the best-by-(Q,H) solution. Following Ensemble Clustering for
# Graphs (Poulin & Theberge 2019), we turn that pool into edge co-association
# weights and re-cluster, which recovers the planted structure on par with ECG
# and outperforms Leiden (its advantage concentrated at high mixing). The leader
# of each consensus community
# is the node most frequently designated a leader across the pool -- an ensemble
# signal a generic detector-plus-centrality pipeline cannot access. A per-node
# confidence map (mean co-association of incident edges) comes for free.

#' LCDA-ECG: ensemble-consensus community and leader detection.
#'
#' Builds a pool of \code{B} randomised LCDA constructions, turns it into edge
#' co-association weights (ECG-style), re-clusters the reweighted graph for the
#' consensus partition, and designates one leader per community from the pool's
#' leader-designation frequencies. Recovers planted structure on par with ECG
#' and outperforms Leiden (advantage concentrated at high mixing), while
#' retaining the joint leader output and a node-confidence map.
#'
#' @details
#' Weighted graphs: a numeric \code{weight} edge attribute is honoured by the
#' modularity objective and the local search inside each pool construction, but
#' the similarity, centrality and NCE leader score remain \emph{structural}
#' (unweighted). The consensus re-clustering uses the ECG co-association weights,
#' not the input weights.
#'
#' @param g an igraph object (undirected, simple).
#' @param B pool size (number of GRASP constructions to ensemble).
#' @param w_min ECG floor weight for 2-core edges; off-2-core edges get exactly
#'   \code{w_min}. Default 0.05, as in Poulin & Theberge (2019).
#' @param alpha_c_range,alpha_s_range bounds of the uniform RCL parameters
#'   sampled per pool member (diversification source).
#' @param variant construction variant, 1 or 2.
#' @param centrality,similarity metric names passed to the construction.
#' @param overlap logical; if \code{TRUE}, also return overlapping community
#'   memberships derived from the co-association (soft) similarity.
#' @param tau overlap threshold in (0,1]: a node joins community \code{c} when
#'   its mean co-association to \code{c} reaches \code{tau} times its home-community
#'   affinity. Only used when \code{overlap = TRUE}.
#' @param verbose logical; show a cli progress bar and a final summary.
#' @param seed integer RNG seed, or \code{NA} to leave the RNG untouched.
#' @return an object of class \code{lcda_ecg_result}: the consensus
#'   \code{membership} (1-based), \code{leaders} (consensus-derived, 1-based),
#'   \code{leaders_central} (top-eigenvector per community, for comparison), a
#'   per-node \code{confidence} vector, the leader-designation counts
#'   \code{lead_count}, the input-graph modularity \code{Q} (weight-aware;
#'   comparable to \code{lcda_grasp}/\code{lcda_gr}), and
#'   \code{Q_consensus_weighted} (modularity under the ECG co-association
#'   weights, i.e. the objective the consensus optimised). When
#'   \code{overlap = TRUE} it
#'   additionally carries \code{overlap_membership} (a length-n list of the
#'   community ids each node belongs to) and \code{is_overlap} (logical, the
#'   bridge nodes).
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_ecg(g, B = 24, overlap = TRUE, tau = 0.6, seed = 1)
#' res$Q
#' which(res$is_overlap)   # bridge nodes
#' @export
lcda_ecg <- function(g, B = 64, w_min = 0.05,
                     alpha_c_range = c(0.1, 0.9), alpha_s_range = c(0.1, 0.5),
                     variant = 1, centrality = "eigen", similarity = "hpi",
                     overlap = FALSE, tau = 0.7,
                     verbose = FALSE, seed = NA_integer_) {
  if (!is.na(seed)) set.seed(seed)
  csr <- .as_csr(g)
  g   <- csr$igraph                                   # simplified, undirected
  n   <- csr$n
  el  <- igraph::as_edgelist(g, names = FALSE)        # m x 2, 1-based
  coassoc    <- numeric(nrow(el))
  lead_count <- integer(n)

  if (verbose) cli::cli_progress_bar("LCDA-ECG pool", total = B, clear = FALSE)
  for (b in seq_len(B)) {
    a_c <- stats::runif(1, alpha_c_range[1], alpha_c_range[2])
    a_s <- stats::runif(1, alpha_s_range[1], alpha_s_range[2])
    sol <- lcda_construct(csr, a_c, a_s, variant = variant,
                          centrality = centrality, similarity = similarity)
    sol <- lcda_repair(csr, sol)
    sol <- lcda_local_search(csr, sol)
    mem <- sol$membership
    coassoc <- coassoc + (mem[el[, 1]] == mem[el[, 2]])   # same-community indicator
    lead_count[sol$leaders] <- lead_count[sol$leaders] + 1L
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) cli::cli_progress_done()
  coassoc <- coassoc / B

  # ECG-style edge weights: 2-core edges blend the floor with co-association.
  core <- igraph::coreness(g)
  in2  <- core[el[, 1]] >= 2 & core[el[, 2]] >= 2
  w    <- ifelse(in2, w_min + (1 - w_min) * coassoc, w_min)

  cons       <- igraph::cluster_louvain(g, weights = w)
  membership <- as.integer(igraph::membership(cons))
  comms      <- sort(unique(membership))
  ev         <- igraph::eigen_centrality(g)$vector

  # consensus leader = most-frequently-designated leader in the pool (eigenvector
  # centrality as a small tie-breaker); plus the plain central leader for reference.
  leaders <- vapply(comms, function(cc) {
    idx <- which(membership == cc)
    idx[which.max(lead_count[idx] + ev[idx] / (max(ev) + 1))]
  }, integer(1))
  leaders_central <- vapply(comms, function(cc) {
    idx <- which(membership == cc); idx[which.max(ev[idx])]
  }, integer(1))

  # per-node confidence: mean co-association over incident edges
  conf <- numeric(n)
  agg  <- tapply(c(coassoc, coassoc), c(el[, 1], el[, 2]), sum)
  conf[as.integer(names(agg))] <- as.numeric(agg)
  deg  <- igraph::degree(g)
  conf <- ifelse(deg > 0, conf / deg, 0)

  # Q on the input graph (honours an input `weight` attribute) is the headline
  # value, directly comparable to the Q reported by lcda_grasp/lcda_gr.
  # Q_consensus_weighted is the modularity under the ECG co-association weights
  # `w` -- the objective the consensus re-clustering actually optimised.
  w_in <- if ("weight" %in% igraph::edge_attr_names(g)) igraph::E(g)$weight else NULL
  Q                    <- igraph::modularity(g, membership, weights = w_in)
  Q_consensus_weighted <- igraph::modularity(g, membership, weights = w)

  # optional overlapping communities from the soft co-association: a node joins
  # any community whose mean co-association reaches tau * its home affinity.
  overlap_membership <- NULL; is_overlap <- NULL
  if (overlap) {
    ends   <- c(el[, 1], el[, 2])
    comm_o <- membership[c(el[, 2], el[, 1])]            # community of the other endpoint
    cv     <- c(coassoc, coassoc)
    s  <- tapply(cv, list(node = ends, comm = comm_o), sum)
    ct <- tapply(cv, list(node = ends, comm = comm_o), length)
    aff <- s / ct                                        # (nodes-with-edges) x communities
    node_ids  <- as.integer(rownames(aff))
    comm_cols <- as.integer(colnames(aff))
    overlap_membership <- as.list(membership)            # default: hard membership
    for (i in seq_along(node_ids)) {
      a <- aff[i, ]; ha <- max(a, na.rm = TRUE)
      joined <- comm_cols[which(!is.na(a) & a >= tau * ha)]
      overlap_membership[[node_ids[i]]] <- sort(union(membership[node_ids[i]], joined))
    }
    is_overlap <- vapply(overlap_membership, length, integer(1)) > 1L
  }

  if (verbose) cli::cli_alert_success(
    "LCDA-ECG: Q = {round(Q, 6)}, {length(leaders)} communities, mean confidence {round(mean(conf), 3)}{if (overlap) paste0(', ', sum(is_overlap), ' overlap nodes') else ''}")

  structure(
    list(membership = membership, leaders = leaders,
         leaders_central = leaders_central, confidence = conf,
         lead_count = lead_count, Q = Q,
         Q_consensus_weighted = Q_consensus_weighted,
         overlap_membership = overlap_membership, is_overlap = is_overlap,
         params = list(B = B, w_min = w_min, variant = variant,
                       centrality = centrality, similarity = similarity,
                       overlap = overlap, tau = tau)),
    class = "lcda_ecg_result")
}

#' @export
print.lcda_ecg_result <- function(x, ...) {
  cli::cli_h1("LCDA-ECG (ensemble consensus) result")
  items <- c(
    "communities"          = "{length(x$leaders)}",
    "modularity Q (input)" = "{format(x$Q, digits = 6)}",
    "Q (consensus weights)" = "{format(x$Q_consensus_weighted, digits = 6)}",
    "pool size B"          = "{x$params$B}",
    "mean node confidence" = "{format(mean(x$confidence), digits = 3)}")
  if (isTRUE(x$params$overlap))
    items <- c(items, "overlap nodes" = "{sum(x$is_overlap)} (tau = {x$params$tau})")
  cli::cli_dl(items)
  invisible(x)
}
