# R/plot_leaders.R
#
# The community-and-leader map: the figure the paper uses to show a partition
# together with the node each community elected. Three entry points, one
# renderer:
#
#   lcda_plot_communities(g)            graph -> detect -> designate -> plot
#   plot(<result>)                      plot an already fitted result
#   ggplot2::autoplot(<result>)         the same figure as a ggplot object
#
# Base graphics is the default so the figure never depends on a Suggested
# package; the ggplot2 path is guarded by requireNamespace().

# --- shared layout / styling -------------------------------------------------

.lcda_layout <- function(g, layout = NULL) {
  if (is.null(layout)) return(igraph::layout_with_fr(g))
  if (is.function(layout)) return(layout(g))
  layout <- as.matrix(layout)
  if (nrow(layout) != igraph::vcount(g) || ncol(layout) < 2)
    cli::cli_abort("{.arg layout} must be a {igraph::vcount(g)} x 2 matrix.")
  layout
}

.lcda_palette <- function(d, palette = "Dynamic") {
  grDevices::hcl.colors(max(d, 2L), palette = palette)[seq_len(d)]
}

.lcda_leader_labels <- function(g, leaders) {
  if ("name" %in% igraph::vertex_attr_names(g))
    as.character(igraph::vertex_attr(g, "name"))[leaders]
  else as.character(leaders)
}

.lcda_check_partition <- function(g, membership, leaders) {
  if (length(membership) != igraph::vcount(g))
    cli::cli_abort("{.arg membership} must have one entry per vertex.")
  if (anyNA(membership)) cli::cli_abort("{.arg membership} must not contain {.val NA}.")
  if (length(leaders) && any(leaders < 1L | leaders > igraph::vcount(g)))
    cli::cli_abort("{.arg leaders} must be 1-based vertex indices.")
  invisible(TRUE)
}

# --- base-graphics renderer --------------------------------------------------

#' Plot a community partition with its leaders highlighted
#'
#' Draws the graph with vertices coloured by community and the leader of each
#' community drawn larger, outlined, and labelled. Intra-community edges take
#' the community colour and inter-community edges are greyed, so the block
#' structure and the elected leaders are both readable at a glance. This is the
#' figure style used for the community-and-leader maps in Ospina et al. (2026).
#'
#' @param g an [igraph::igraph] object (undirected, simple).
#' @param membership integer vector of 1-based community ids, one per vertex.
#' @param leaders integer vector of 1-based leader vertex indices.
#' @param layout optional layout: a two-column matrix with one row per vertex,
#'   or a layout function such as [igraph::layout_with_kk]. Defaults to
#'   Fruchterman-Reingold. Pass an explicit layout for a reproducible figure.
#' @param leader_labels logical; label the leader vertices (with the vertex
#'   `name` attribute when present, otherwise the vertex index).
#' @param legend logical; draw a community/leader legend. At most
#'   `legend_max` communities are listed, with a "+k more" entry.
#' @param legend_max maximum number of communities listed in the legend.
#' @param vertex_size,leader_size plotting sizes for ordinary vertices and for
#'   leaders.
#' @param palette an [grDevices::hcl.colors()] qualitative palette name.
#' @param shade_edges logical; colour intra-community edges by community and
#'   grey out the inter-community ones.
#' @param mark_communities logical; additionally draw shaded hulls around the
#'   communities (via `mark.groups`).
#' @param main plot title; `NULL` for a sensible default.
#' @param ... further arguments passed to [igraph::plot.igraph()].
#' @return invisibly, a list with the `layout` used and the `colors` per
#'   community (so a caller can reuse them for a companion figure).
#' @seealso [lcda_plot_communities()] for the end-to-end pipeline and
#'   [plot.lcda_grasp_result()] to plot a fitted result directly.
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_grasp(g, B = 20, seed = 1)
#' set.seed(1)
#' plot_partition(g, res$best$membership, res$best$leaders)
#' @export
plot_partition <- function(g, membership, leaders,
                           layout = NULL,
                           leader_labels = TRUE,
                           legend = TRUE, legend_max = 12L,
                           vertex_size = 4, leader_size = 9,
                           palette = "Dynamic",
                           shade_edges = TRUE,
                           mark_communities = FALSE,
                           main = NULL, ...) {
  .lcda_check_partition(g, membership, leaders)
  membership <- as.integer(membership)
  d <- max(membership)
  cols <- .lcda_palette(d, palette)
  n <- igraph::vcount(g)

  is_leader <- rep(FALSE, n); is_leader[leaders] <- TRUE
  vsize <- rep(vertex_size, n); vsize[leaders] <- leader_size
  vfill <- grDevices::adjustcolor(cols[membership], alpha.f = 0.85)
  vframe <- ifelse(is_leader, "black", grDevices::adjustcolor("grey30", alpha.f = 0.5))

  vlabel <- rep(NA_character_, n)
  if (isTRUE(leader_labels) && length(leaders))
    vlabel[leaders] <- .lcda_leader_labels(g, leaders)

  ecol <- "grey85"
  if (isTRUE(shade_edges) && igraph::ecount(g) > 0) {
    el <- igraph::as_edgelist(g, names = FALSE)
    intra <- membership[el[, 1]] == membership[el[, 2]]
    ecol <- ifelse(intra,
                   grDevices::adjustcolor(cols[membership[el[, 1]]], alpha.f = 0.55),
                   grDevices::adjustcolor("grey60", alpha.f = 0.35))
  }

  lay <- .lcda_layout(g, layout)
  groups <- if (isTRUE(mark_communities)) split(seq_len(n), membership) else NULL
  if (is.null(main))
    main <- sprintf("%d communities, %d leaders", d, length(leaders))

  igraph::plot.igraph(
    g, layout = lay,
    vertex.color = vfill, vertex.size = vsize,
    vertex.frame.color = vframe,
    vertex.label = vlabel, vertex.label.cex = 0.7,
    vertex.label.color = "black", vertex.label.family = "sans",
    vertex.label.dist = 0, vertex.label.font = 2,
    edge.color = ecol, edge.width = 0.8,
    mark.groups = groups,
    mark.col = if (is.null(groups)) NULL else grDevices::adjustcolor(cols, alpha.f = 0.10),
    mark.border = if (is.null(groups)) NULL else grDevices::adjustcolor(cols, alpha.f = 0.35),
    main = main, ...)

  if (isTRUE(legend) && d >= 1) {
    lab_of <- rep(NA_character_, d)
    if (length(leaders)) lab_of[membership[leaders]] <- .lcda_leader_labels(g, leaders)
    k <- min(d, legend_max)
    txt <- sprintf("C%d (leader %s)", seq_len(k), lab_of[seq_len(k)])
    fill <- cols[seq_len(k)]
    if (d > k) {
      txt <- c(txt, sprintf("+%d more", d - k))
      fill <- c(fill, NA)
    }
    graphics::legend("bottomleft", legend = txt, fill = fill, border = NA,
                     bty = "n", cex = 0.7)
  }

  invisible(list(layout = lay, colors = cols))
}

# --- end-to-end pipeline: graph in, figure out -------------------------------

#' Detect communities and leaders on a graph, then plot them
#'
#' The one-call path from a bare graph to the paper's community-and-leader
#' figure: run one of the LCDA algorithms, take the partition and the leader it
#' designates for each community, and render them with [plot_partition()].
#'
#' @param graph an [igraph::igraph] object (undirected, simple).
#' @param method which algorithm to run: `"gr"` ([lcda_gr()], the reactive
#'   variant and the default), `"grasp"` ([lcda_grasp()]), or `"ecg"`
#'   ([lcda_ecg()], the ensemble consensus, which also yields node confidences).
#' @param args a named list of extra arguments for the chosen algorithm, e.g.
#'   `list(B = 100, variant = 2)`.
#' @param seed integer RNG seed for the run, or `NA` to leave the RNG untouched.
#'   Also seeds the layout, so the figure is reproducible.
#' @param plot logical; set to `FALSE` to fit and return without drawing.
#' @param verbose logical; passed to the algorithm.
#' @param ... further arguments passed to [plot_partition()] (`layout`,
#'   `legend`, `mark_communities`, ...).
#' @return invisibly, the fitted result object (an `lcda_gr_result`,
#'   `lcda_grasp_result`, or `lcda_ecg_result`), with the layout used attached
#'   as the attribute `"lcda_layout"`. Feed it straight to [lcda_metrics()].
#' @seealso [lcda_metrics()], [plot_partition()].
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_plot_communities(g, method = "grasp", args = list(B = 20), seed = 1)
#' lcda_metrics(res, level = "leader")
#' @export
lcda_plot_communities <- function(graph,
                                  method = c("gr", "grasp", "ecg"),
                                  args = list(),
                                  seed = NA_integer_,
                                  plot = TRUE,
                                  verbose = FALSE,
                                  ...) {
  method <- match.arg(method)
  .check_graph(graph)
  if (!is.list(args)) cli::cli_abort("{.arg args} must be a named list.")
  seed <- .check_seed(seed)

  fn <- switch(method, gr = lcda_gr, grasp = lcda_grasp, ecg = lcda_ecg)
  base_args <- list(verbose = verbose, seed = seed)
  base_args <- base_args[setdiff(names(base_args), names(args))]  # args wins
  res <- do.call(fn, c(list(graph), base_args, args))

  sol <- .lcda_solution(res)
  if (!is.na(seed)) set.seed(seed)          # reproducible layout
  lay <- NULL
  if (isTRUE(plot)) {
    drawn <- plot_partition(sol$graph, sol$membership, sol$leaders, ...)
    lay <- drawn$layout
  }
  attr(res, "lcda_layout") <- lay
  invisible(res)
}

# --- S3 plot methods ---------------------------------------------------------

.lcda_plot_result <- function(x, y, ...) {
  sol <- .lcda_solution(x, graph = y)
  plot_partition(sol$graph, sol$membership, sol$leaders, ...)
}

#' Plot a fitted LCDA result as a community-and-leader map
#'
#' `plot()` methods for the objects returned by [lcda_grasp()], [lcda_gr()] and
#' [lcda_ecg()]. The result carries the (simplified) graph it was fitted on, so
#' `plot(res)` is enough; pass a graph as the second argument to draw the same
#' partition on a different vertex layout or on the unsimplified original.
#'
#' @param x a fitted `lcda_grasp_result`, `lcda_gr_result` or `lcda_ecg_result`.
#' @param y optional [igraph::igraph] to draw on; defaults to the graph stored
#'   in `x`.
#' @param ... further arguments passed to [plot_partition()].
#' @return invisibly, a list with the `layout` used and the community `colors`.
#' @seealso [plot_partition()], [lcda_plot_communities()], [lcda_metrics()].
#' @examples
#' g <- igraph::make_graph("Zachary")
#' res <- lcda_grasp(g, B = 20, seed = 1)
#' set.seed(1)
#' plot(res)
#' @name plot.lcda
#' @export
plot.lcda_grasp_result <- function(x, y = NULL, ...) .lcda_plot_result(x, y, ...)

#' @rdname plot.lcda
#' @export
plot.lcda_gr_result <- function(x, y = NULL, ...) .lcda_plot_result(x, y, ...)

#' @rdname plot.lcda
#' @export
plot.lcda_ecg_result <- function(x, y = NULL, ...) .lcda_plot_result(x, y, ...)

# --- ggplot2 renderer (Suggests-guarded) -------------------------------------

.lcda_ggplot <- function(g, membership, leaders, layout = NULL,
                         leader_labels = TRUE, palette = "Dynamic",
                         vertex_size = 1.8, leader_size = 4.5,
                         title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    cli::cli_abort(c("{.pkg ggplot2} is required for {.fn autoplot}.",
                     i = "Use {.fn plot} for the base-graphics version."))
  .lcda_check_partition(g, membership, leaders)
  membership <- as.integer(membership)
  d <- max(membership)
  cols <- .lcda_palette(d, palette)
  lay <- .lcda_layout(g, layout)

  nodes <- data.frame(
    x = lay[, 1], y = lay[, 2],
    community = factor(membership, levels = seq_len(d)),
    is_leader = seq_len(igraph::vcount(g)) %in% leaders)
  nodes$label <- ifelse(nodes$is_leader,
                        .lcda_leader_labels(g, seq_len(nrow(nodes))), NA_character_)

  el <- igraph::as_edgelist(g, names = FALSE)
  edges <- data.frame(
    x = lay[el[, 1], 1], y = lay[el[, 1], 2],
    xend = lay[el[, 2], 1], yend = lay[el[, 2], 2],
    intra = membership[el[, 1]] == membership[el[, 2]])

  lead_df <- nodes[nodes$is_leader, , drop = FALSE]
  if (is.null(title))
    title <- sprintf("%d communities, %d leaders", d, length(leaders))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = intra),
      colour = "grey55", linewidth = 0.25, show.legend = FALSE) +
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.18, "TRUE" = 0.5)) +
    ggplot2::geom_point(
      data = nodes[!nodes$is_leader, , drop = FALSE],
      ggplot2::aes(x = x, y = y, colour = community),
      size = vertex_size, alpha = 0.85) +
    ggplot2::geom_point(
      data = lead_df,
      ggplot2::aes(x = x, y = y, fill = community),
      shape = 21, colour = "black", stroke = 0.7, size = leader_size) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, colour = "community") +
    ggplot2::theme_void(base_size = 11)

  if (isTRUE(leader_labels) && nrow(lead_df))
    p <- p + ggplot2::geom_text(
      data = lead_df, ggplot2::aes(x = x, y = y, label = label),
      vjust = -1.1, size = 3, fontface = "bold", colour = "grey15")
  p
}

#' Community-and-leader map as a ggplot object
#'
#' [ggplot2::autoplot()] methods for the fitted LCDA results: the same figure as
#' [plot.lcda_grasp_result()], returned as a `ggplot` object so it can be
#' themed, faceted, or saved with [ggplot2::ggsave()]. Requires the suggested
#' package \pkg{ggplot2}.
#'
#' @param object a fitted `lcda_grasp_result`, `lcda_gr_result` or
#'   `lcda_ecg_result`.
#' @param graph optional [igraph::igraph] to draw on; defaults to the graph
#'   stored in `object`.
#' @param layout optional layout matrix or layout function; defaults to
#'   Fruchterman-Reingold.
#' @param ... further arguments passed to the internal renderer
#'   (`leader_labels`, `palette`, `vertex_size`, `leader_size`, `title`).
#' @return a `ggplot` object.
#' @seealso [plot.lcda_grasp_result()], [lcda_plot_communities()].
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   g <- igraph::make_graph("Zachary")
#'   res <- lcda_grasp(g, B = 20, seed = 1)
#'   set.seed(1)
#'   ggplot2::autoplot(res)
#' }
#' @name autoplot.lcda
#' @method autoplot lcda_grasp_result
#' @export
autoplot.lcda_grasp_result <- function(object, graph = NULL, layout = NULL, ...) {
  s <- .lcda_solution(object, graph)
  .lcda_ggplot(s$graph, s$membership, s$leaders, layout = layout, ...)
}

#' @rdname autoplot.lcda
#' @method autoplot lcda_gr_result
#' @export
autoplot.lcda_gr_result <- autoplot.lcda_grasp_result

#' @rdname autoplot.lcda
#' @method autoplot lcda_ecg_result
#' @export
autoplot.lcda_ecg_result <- autoplot.lcda_grasp_result
