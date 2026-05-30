# R/stats_and_plots.R
#
# Statistical comparison procedures used in section 5 of the paper, plus helper
# plotting functions for partitions, GRASP trajectories, and the Reactive p_k.

#' Two-stage non-parametric comparison: Kruskal-Wallis omnibus followed by
#' pairwise permutation tests with Bonferroni correction.
#'
#' Mirrors the procedure described in section 5.1 of the paper.
#'
#' @param value numeric vector of observations.
#' @param group factor or character vector of group labels.
#' @param B_perm number of permutations per pairwise test.
#' @param alpha global significance level.
#' @param verbose logical; emit a cli trace of the decisions.
#' @return list(omnibus, pairwise, letters, decision); `pairwise` is a tibble and
#'   `letters` is a compact-letter display compatible with the paper's tables.
#' @export
kruskal_then_permute <- function(value, group, B_perm = 1e4, alpha = 0.05,
                                 verbose = FALSE) {
  group <- factor(group)
  lv <- levels(group)
  k <- length(lv)
  if (k < 2) cli::cli_abort("Need at least 2 groups, got {k}.")

  omnibus <- stats::kruskal.test(value, group)
  if (verbose) cli::cli_alert_info("Kruskal-Wallis p = {signif(omnibus$p.value, 3)}")
  if (omnibus$p.value > alpha) {
    if (verbose) cli::cli_alert_info("Omnibus not rejected at alpha = {alpha}; all groups share a letter.")
    return(list(omnibus = omnibus, pairwise = NULL,
                letters = stats::setNames(rep("a", k), lv),
                decision = "omnibus_not_rejected"))
  }

  pairs <- utils::combn(k, 2)
  np <- ncol(pairs)
  p_raw <- numeric(np); diff <- numeric(np)
  for (i in seq_len(np)) {
    a <- lv[pairs[1, i]]; b <- lv[pairs[2, i]]
    pt <- permutation_test(value[group == a], value[group == b], B = B_perm)
    p_raw[i] <- pt$p.value; diff[i] <- pt$diff
  }
  p_adj <- stats::p.adjust(p_raw, method = "bonferroni")
  pw <- tibble::tibble(
    group1 = lv[pairs[1, ]], group2 = lv[pairs[2, ]],
    diff = diff, p_raw = p_raw, p_adj = p_adj, significant = p_adj < alpha)

  if (verbose) cli::cli_alert_success("{sum(pw$significant)}/{np} pairwise contrasts significant after Bonferroni.")
  list(omnibus = omnibus, pairwise = pw,
       letters = .cld_from_pairwise(lv, pw, alpha = alpha),
       decision = "omnibus_rejected")
}

#' Two-sample permutation test for a difference in means.
#' Robust to ties and to zero-variance pools (as discussed in section 5.1).
#' @param x,y numeric vectors of observations for the two groups.
#' @param B number of permutations.
#' @param verbose logical; emit a cli trace.
#' @return list with the observed mean difference and the permutation p-value.
#' @export
permutation_test <- function(x, y, B = 1e4, verbose = FALSE) {
  obs <- mean(x) - mean(y)
  pool <- c(x, y); nx <- length(x); n <- length(pool)
  if (n < 2) return(list(diff = obs, p.value = NA_real_))
  if (stats::var(pool) == 0) {
    if (verbose) cli::cli_alert_info("Zero-variance pool; p = 1 by construction.")
    return(list(diff = obs, p.value = 1.0))
  }
  cnt <- 0L
  for (b in seq_len(B)) {
    idx <- sample.int(n, nx)
    if (abs(mean(pool[idx]) - mean(pool[-idx])) >= abs(obs)) cnt <- cnt + 1L
  }
  if (verbose) cli::cli_alert_info("observed diff = {round(obs, 4)}, permutation p = {round((cnt + 1) / (B + 1), 4)}")
  list(diff = obs, p.value = (cnt + 1) / (B + 1))
}

# Compact letter display: groups sharing a letter are not significantly
# different at level alpha (greedy assignment).
.cld_from_pairwise <- function(lv, pw, alpha = 0.05) {
  k <- length(lv)
  not_diff <- matrix(TRUE, k, k); rownames(not_diff) <- colnames(not_diff) <- lv
  for (i in seq_len(nrow(pw))) {
    if (pw$significant[i]) {
      not_diff[pw$group1[i], pw$group2[i]] <- FALSE
      not_diff[pw$group2[i], pw$group1[i]] <- FALSE
    }
  }
  letters_seq <- c(letters, LETTERS)
  group_letter <- list(); next_idx <- 1L
  for (i in seq_len(k)) {
    candidates <- character(0)
    for (j in seq_len(i - 1)) if (not_diff[i, j]) candidates <- c(candidates, group_letter[[lv[j]]])
    if (length(candidates) == 0) {
      group_letter[[lv[i]]] <- letters_seq[next_idx]; next_idx <- next_idx + 1L
    } else {
      group_letter[[lv[i]]] <- unique(candidates)
    }
  }
  vapply(lv, function(g) paste(sort(unique(group_letter[[g]])), collapse = ""), character(1))
}

# --- plotting helpers ------------------------------------------------

#' Plot a leader-community partition (leaders drawn larger).
#' @param g an igraph object (undirected, simple).
#' @param membership integer vector of 1-based community ids.
#' @param leaders integer vector of 1-based leader vertex indices.
#' @param ... further arguments passed to [igraph::plot.igraph()].
#' @return invisibly, `NULL` (called for its plotting side effect).
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_grasp(g, B = 20, seed = 1)
#' plot_partition(g, res$best$membership, res$best$leaders)
#' @export
plot_partition <- function(g, membership, leaders, ...) {
  cols <- grDevices::hcl.colors(max(membership), palette = "Dynamic")
  vsize <- rep(4, igraph::vcount(g)); vsize[leaders] <- 9
  igraph::plot.igraph(g, vertex.color = cols[membership],
                      vertex.size = vsize, vertex.label = NA, ...)
}

#' Plot the Q trajectory across GRASP iterations, with running maximum.
#' @param res an `lcda_grasp_result` or `lcda_gr_result` object.
#' @param ... further arguments passed to [plot()].
#' @return invisibly, a data frame of per-iteration Q and running maximum.
#' @export
plot_grasp_trajectory <- function(res, ...) {
  if (!inherits(res, c("lcda_grasp_result", "lcda_gr_result")))
    cli::cli_abort("{.arg res} must be an {.cls lcda_grasp_result} or {.cls lcda_gr_result}.")
  Q <- res$trace_Q; Qmax <- cummax(Q)
  plot(seq_along(Q), Q, type = "p", pch = 19, col = "grey60",
       xlab = "iteration", ylab = "modularity",
       main = "GRASP trajectory: per-iter Q and running best", ...)
  lines(seq_along(Q), Qmax, lwd = 2)
  invisible(data.frame(iter = seq_along(Q), Q = Q, Q_running_max = Qmax))
}

#' Plot the evolution of selection probabilities p_k across iterations.
#' @param res an `lcda_gr_result` object.
#' @param ... further arguments passed to [matplot()].
#' @return invisibly, the p_k history matrix.
#' @export
plot_reactive_pk <- function(res, ...) {
  if (!inherits(res, "lcda_gr_result")) cli::cli_abort("{.arg res} must be an {.cls lcda_gr_result}.")
  pk <- res$pk_history; B <- nrow(pk); m <- ncol(pk)
  matplot(seq_len(B), pk, type = "l", lty = 1,
          xlab = "iteration", ylab = "p_k",
          main = "Reactive update: probability mass over pool",
          col = grDevices::hcl.colors(m, palette = "Zissou 1"), ...)
  invisible(pk)
}

# --- convenience graph constructors ---------------------------------

#' Build an igraph from an integer edgelist (1-based), undirected simple.
#' @param edges a two-column matrix (or coercible) of 1-based vertex pairs.
#' @param n optional vertex count; padded with isolates if larger than observed.
#' @return an undirected simple igraph object.
#' @export
graph_from_edgelist_simple <- function(edges, n = NULL) {
  if (is.null(n)) n <- max(edges)
  g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  if (igraph::vcount(g) < n) g <- igraph::add_vertices(g, n - igraph::vcount(g))
  g
}

#' Coerce a variety of inputs to igraph.
#'
#' A pure constructor (no logging).
#' @param x an igraph object, a two-column edgelist matrix, or an adjacency matrix.
#' @return an undirected igraph object.
#' @examples
#' el <- matrix(c(1, 2, 2, 3, 3, 1), ncol = 2, byrow = TRUE)
#' as_graph(el)
#' @export
as_graph <- function(x) {
  if (inherits(x, "igraph")) return(x)
  if (is.matrix(x) && ncol(x) == 2) return(graph_from_edgelist_simple(x))
  if (methods::is(x, "Matrix") || is.matrix(x)) {
    return(igraph::graph_from_adjacency_matrix(x, mode = "undirected", diag = FALSE, weighted = NULL))
  }
  cli::cli_abort("Don't know how to coerce an object of class {.cls {class(x)[1]}}.")
}
