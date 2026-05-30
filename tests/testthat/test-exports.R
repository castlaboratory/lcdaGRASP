# Smoke + invariant tests for the remaining exported API, to exercise the
# public surface (and the documented examples' code paths).

test_that("lcda_gr runs and tracks the reactive pool", {
  set.seed(1)
  g <- igraph::make_graph("Zachary")
  res <- lcda_gr(g, B = 30, seed = 1)
  expect_s3_class(res, "lcda_gr_result")
  expect_equal(length(res$best$membership), igraph::vcount(g))
  expect_true(all(res$best$leaders %in% seq_len(igraph::vcount(g))))
  expect_equal(nrow(res$pk_history), 30L)
  expect_true(res$best$Q > 0)
  expect_invisible(print(res))
})

test_that("scores agree with references and respect bounds", {
  g <- igraph::make_graph("Zachary")
  cl <- igraph::cluster_louvain(g); mem <- as.integer(igraph::membership(cl))
  # modularity_score matches igraph::modularity
  expect_equal(modularity_score(g, mem), igraph::modularity(g, mem), tolerance = 1e-8)
  # NCE scores are bounded in [0, 1]
  h <- nce_score(g, leaders = c(1, 34))
  expect_true(h >= 0 && h <= 1)
  hl <- nce_local_score(g, mem, leaders = c(1, 34))
  expect_true(hl >= 0 && hl <= 1)
})

test_that("lex_dominates implements lexicographic (Q, then H) order", {
  expect_true(lex_dominates(0.5, 0.1, 0.4, 0.9))   # higher Q wins
  expect_true(lex_dominates(0.5, 0.9, 0.5, 0.5))   # tie in Q, higher H wins
  expect_false(lex_dominates(0.4, 0.9, 0.5, 0.1))  # lower Q loses
  expect_false(lex_dominates(0.5, 0.5, 0.5, 0.9))  # tie in Q, lower H loses
})

test_that("centrality wrappers return one finite score per node", {
  g <- igraph::make_graph("Zachary")
  csr <- lcdaGRASP:::.as_csr(g)
  for (f in list(centrality_eigen, centrality_betweenness, centrality_closeness)) {
    v <- f(csr)
    expect_equal(length(v), igraph::vcount(g))
    expect_true(all(is.finite(v)))
  }
})

test_that("construction phases produce a valid partition with one leader per community", {
  g <- igraph::make_graph("Zachary")
  csr <- lcdaGRASP:::.as_csr(g)
  sol <- lcda_construct(csr, 0.1, 0.3)
  sol <- lcda_repair(csr, sol)
  sol <- lcda_local_search(csr, sol)
  expect_equal(length(sol$membership), igraph::vcount(g))
  expect_equal(length(sol$leaders), length(unique(sol$membership)))
})

test_that("statistical comparison helpers run and return sane structures", {
  set.seed(1)
  x <- rnorm(20, 0); y <- rnorm(20, 1)
  pt <- permutation_test(x, y, B = 200)
  expect_true(pt$p.value >= 0 && pt$p.value <= 1)
  val <- c(x, y); grp <- rep(c("a", "b"), each = 20)
  kp <- kruskal_then_permute(val, grp, B_perm = 200)
  expect_true(is.list(kp))
})

test_that("variant 2 and the VNMI path (n > 300) run and return valid output", {
  set.seed(7)
  # n > 300 triggers the VNMI local search; variant 2 recomputes centrality.
  g <- igraph::sample_gnp(360, p = 0.02)
  res <- lcda_grasp(g, variant = 2, B = 3, centrality = "eigen", similarity = "hpi", seed = 1)
  expect_s3_class(res, "lcda_grasp_result")
  expect_equal(length(res$best$membership), igraph::vcount(g))
  expect_true(all(res$best$leaders %in% seq_len(igraph::vcount(g))))
  expect_true(res$best$Q > 0)
})

test_that("graph constructors and plotters run without error", {
  el <- matrix(c(1, 2, 2, 3, 3, 1, 3, 4), ncol = 2, byrow = TRUE)
  g <- graph_from_edgelist_simple(el)
  expect_s3_class(g, "igraph")
  expect_s3_class(as_graph(el), "igraph")
  expect_identical(as_graph(g), g)
  # plotters are called for their side effect; just ensure they do not error
  res <- lcda_grasp(igraph::make_graph("Zachary"), B = 10, seed = 1)
  pdf(NULL)
  expect_error(plot_partition(igraph::make_graph("Zachary"), res$best$membership, res$best$leaders), NA)
  expect_error(plot_grasp_trajectory(res), NA)
  gr <- lcda_gr(igraph::make_graph("Zachary"), B = 12, seed = 1)
  expect_error(plot_reactive_pk(gr), NA)
  dev.off()
})
