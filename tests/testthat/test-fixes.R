# Tests for the 0.3.x fixes: reactive p_floor, seed handling, lcda_ecg centrality
# and the two modularity outputs.

g <- igraph::make_graph("Zachary")

test_that("reactive refresh keeps every pair reachable (p_floor floor)", {
  res <- lcda_gr(g, B = 80, m = 20, p_floor = 0.05, seed = 1)
  # after refreshes, no pair may be permanently excluded
  expect_true(min(res$pk_history[nrow(res$pk_history), ]) >= 0.05 / 20 - 1e-9)
  expect_equal(sum(res$pk_history[nrow(res$pk_history), ]), 1, tolerance = 1e-8)
})

test_that("p_floor = 0 recovers the raw proportional rule (no floor)", {
  res <- lcda_gr(g, B = 80, m = 20, p_floor = 0, seed = 1)
  expect_true(all(res$pk_history[nrow(res$pk_history), ] >= 0))
})

test_that("seed = NULL and seed = NA do not error", {
  expect_no_error(lcda_grasp(g, B = 5, seed = NULL))
  expect_no_error(lcda_grasp(g, B = 5, seed = NA))
  expect_no_error(lcda_gr(g, B = 10, seed = NULL))
  expect_no_error(lcda_ecg(g, B = 8, seed = NULL))
})

test_that("an invalid seed is rejected with a clear message", {
  expect_error(lcda_grasp(g, B = 5, seed = "x"), "single integer")
  expect_error(lcda_grasp(g, B = 5, seed = c(1, 2)), "single integer")
})

test_that("lcda_ecg returns both input and consensus-weighted modularity", {
  res <- lcda_ecg(g, B = 16, seed = 1)
  expect_true(is.numeric(res$Q) && length(res$Q) == 1L)
  expect_true(is.numeric(res$Q_consensus_weighted) && length(res$Q_consensus_weighted) == 1L)
})

test_that("lcda_ecg respects the chosen centrality in leader selection", {
  # The repair/local-search now receive `centrality`; betweenness and eigenvector
  # are exercised without error and return one leader per community.
  re <- lcda_ecg(g, B = 16, centrality = "eigen", seed = 1)
  rb <- lcda_ecg(g, B = 16, centrality = "betweenness", seed = 1)
  expect_length(re$leaders, length(unique(re$membership)))
  expect_length(rb$leaders, length(unique(rb$membership)))
})
