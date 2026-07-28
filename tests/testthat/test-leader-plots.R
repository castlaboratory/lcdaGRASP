# The community-and-leader map. Plots are checked for "runs and returns the
# right structure", not for pixels.

test_that("plot_partition draws and returns its layout and colours", {
  g <- igraph::make_graph("Zachary")
  res <- lcda_grasp(g, B = 10, seed = 1)
  pdf(NULL); on.exit(dev.off(), add = TRUE)

  out <- plot_partition(g, res$best$membership, res$best$leaders)
  expect_type(out, "list")
  expect_equal(dim(out$layout), c(igraph::vcount(g), 2L))
  expect_length(out$colors, max(res$best$membership))

  # explicit layout is honoured, and the optional flourishes run
  lay <- igraph::layout_in_circle(g)
  out2 <- plot_partition(g, res$best$membership, res$best$leaders,
                         layout = lay, mark_communities = TRUE,
                         leader_labels = FALSE, legend = FALSE)
  expect_equal(out2$layout, lay)
  expect_error(
    plot_partition(g, res$best$membership, res$best$leaders,
                   layout = igraph::layout_with_kk, legend_max = 1L), NA)
})

test_that("plot_partition validates its arguments", {
  g <- igraph::make_graph("Zachary")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot_partition(g, rep(1L, 5), 1L), "one entry per vertex")
  expect_error(plot_partition(g, rep(1L, 34), 999L), "1-based")
  expect_error(plot_partition(g, rep(1L, 34), 1L, layout = matrix(0, 2, 2)),
               "matrix")
})

test_that("plot() methods work on all three result classes", {
  g <- igraph::make_graph("Zachary")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  for (res in list(lcda_grasp(g, B = 8, seed = 1),
                   lcda_gr(g, B = 10, seed = 1),
                   lcda_ecg(g, B = 8, seed = 1))) {
    out <- plot(res)                       # graph taken from the result
    expect_equal(dim(out$layout), c(igraph::vcount(g), 2L))
    expect_error(plot(res, g), NA)         # explicit graph also accepted
  }
})

test_that("lcda_plot_communities runs the full graph -> communities -> leaders pipeline", {
  g <- igraph::make_graph("Zachary")
  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- lcda_plot_communities(g, method = "grasp", args = list(B = 8), seed = 1)
  expect_s3_class(res, "lcda_grasp_result")
  expect_equal(dim(attr(res, "lcda_layout")), c(igraph::vcount(g), 2L))
  # the fitted object feeds straight into the metric surface
  expect_s3_class(lcda_metrics(res, level = "leader"), "tbl_df")

  res_gr <- lcda_plot_communities(g, args = list(B = 10), seed = 1)  # default "gr"
  expect_s3_class(res_gr, "lcda_gr_result")
  res_ecg <- lcda_plot_communities(g, method = "ecg", args = list(B = 8),
                                   seed = 1, plot = FALSE)
  expect_s3_class(res_ecg, "lcda_ecg_result")
  expect_null(attr(res_ecg, "lcda_layout"))
  expect_error(lcda_plot_communities(g, args = "not a list"), "list")
})

test_that("autoplot() returns a ggplot for every result class", {
  skip_if_not_installed("ggplot2")
  g <- igraph::make_graph("Zachary")
  for (res in list(lcda_grasp(g, B = 8, seed = 1),
                   lcda_gr(g, B = 10, seed = 1),
                   lcda_ecg(g, B = 8, seed = 1))) {
    p <- ggplot2::autoplot(res)
    expect_s3_class(p, "ggplot")
  }
  p <- ggplot2::autoplot(lcda_grasp(g, B = 8, seed = 1),
                         layout = igraph::layout_in_circle(g),
                         leader_labels = FALSE)
  expect_s3_class(p, "ggplot")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(print(p), NA)
})
