# Invariants for lcda_metrics(). Metaheuristic output varies with the seed, so
# every assertion below is a *range*, an *identity*, or a *conservation law* --
# never a hard-coded number.

zachary <- function() igraph::make_graph("Zachary")

test_that("level = 'overall' returns tidy long metrics inside their valid ranges", {
  g <- zachary()
  res <- lcda_grasp(g, B = 15, seed = 1)
  m <- lcda_metrics(res)

  expect_s3_class(m, "tbl_df")
  expect_named(m, c("algorithm", "scope", "metric", "value"))
  expect_true(all(is.finite(m$value)))
  expect_false(anyDuplicated(paste(m$scope, m$metric)) > 0)

  val <- function(k) m$value[m$metric == k]
  expect_true(val("Q") >= -1 && val("Q") <= 1)
  expect_true(val("H_global") >= 0 && val("H_global") <= 1)
  expect_true(val("H_local") >= 0 && val("H_local") <= 1)
  expect_equal(val("n_nodes"), igraph::vcount(g))
  expect_equal(val("n_edges"), igraph::ecount(g))
  expect_true(val("n_communities") >= 1)
  expect_equal(val("n_leaders"), val("n_communities"))
  expect_true(val("leader_coverage") > 0 && val("leader_coverage") <= 1)
  expect_true(val("size_min") <= val("size_max"))
  expect_equal(val("leaders_derived"), 0)   # LCDA results carry real leaders
})

test_that("overall metrics agree with the values the result object already stores", {
  g <- zachary()
  res <- lcda_gr(g, B = 20, seed = 3)
  m <- lcda_metrics(res)
  val <- function(k) m$value[m$metric == k]

  # lcda_metrics must not keep a second set of books
  expect_equal(val("Q"), res$best$Q, tolerance = 1e-10)
  expect_equal(val("H_global"), res$best$H, tolerance = 1e-10)
  expect_equal(val("elapsed_sec"), res$elapsed, tolerance = 1e-10)
  expect_equal(val("best_iter"), res$best$iter)
  expect_equal(val("B"), res$params$B)
  expect_equal(val("lex_decisive_n"), length(res$lex_decisive_iters))
  expect_equal(val("n_communities"), length(unique(res$best$membership)))
  # and Q must agree with the independent igraph reference
  expect_equal(val("Q"), igraph::modularity(res$graph, res$best$membership),
               tolerance = 1e-8)
})

test_that("recovery metrics appear only with a truth and match igraph::compare", {
  g <- zachary()
  res <- lcda_grasp(g, B = 10, seed = 5)
  truth <- c(rep(1L, 17), rep(2L, 17))

  expect_false("recovery" %in% lcda_metrics(res)$scope)
  m <- lcda_metrics(res, truth = truth)
  val <- function(k) m$value[m$metric == k]
  expect_true("recovery" %in% m$scope)
  expect_true(val("nmi") >= 0 && val("nmi") <= 1)
  expect_true(val("ari") >= -1 && val("ari") <= 1)
  expect_true(val("vi") >= 0)
  expect_equal(val("nmi"),
               igraph::compare(truth, res$best$membership, method = "nmi"),
               tolerance = 1e-10)
  expect_equal(val("n_communities_truth"), 2)

  # a perfect partition scores 1; truth given as a factor works too
  perfect <- lcda_metrics(list(membership = truth), g, truth = factor(truth))
  expect_equal(perfect$value[perfect$metric == "nmi"], 1, tolerance = 1e-8)
  expect_equal(perfect$value[perfect$metric == "ari"], 1, tolerance = 1e-8)
})

test_that("community-level metrics conserve the partition and decompose Q", {
  g <- zachary()
  res <- lcda_grasp(g, B = 15, seed = 7)
  cm <- lcda_metrics(res, level = "community")

  expect_equal(nrow(cm), length(unique(res$best$membership)))
  expect_equal(sum(cm$size), igraph::vcount(g))
  # the additive modularity decomposition must sum to Q
  expect_equal(sum(cm$q_contribution), res$best$Q, tolerance = 1e-8)
  # every edge is either internal to a community or on a boundary
  expect_equal(sum(cm$internal_edges) + sum(cm$boundary_edges) / 2,
               igraph::ecount(g), tolerance = 1e-8)
  expect_false(anyNA(cm$leader))
  expect_true(all(cm$leader %in% res$best$leaders))
  # each community's leader really belongs to that community
  expect_equal(res$best$membership[cm$leader], cm$community)
  ok <- !is.na(cm$conductance)
  expect_true(all(cm$conductance[ok] >= 0 & cm$conductance[ok] <= 1))
  ok <- !is.na(cm$internal_density)
  expect_true(all(cm$internal_density[ok] >= 0 & cm$internal_density[ok] <= 1))
})

test_that("leader-level metrics are consistent with the graph", {
  g <- zachary()
  res <- lcda_ecg(g, B = 12, overlap = TRUE, seed = 2)
  lm <- lcda_metrics(res, level = "leader")

  expect_equal(nrow(lm), length(res$leaders))
  expect_false(anyDuplicated(lm$leader) > 0)
  expect_equal(sort(lm$community), sort(unique(res$membership)))
  expect_equal(lm$degree, as.numeric(igraph::degree(g, lm$leader)))
  expect_true(all(lm$degree_within <= lm$degree))
  expect_true(all(lm$degree_between >= 0))
  expect_true(all(lm$nce_node >= 0 & lm$nce_node <= 1))
  expect_true(all(lm$degree_pctile_in_community > 0 &
                    lm$degree_pctile_in_community <= 1))
  expect_true(all(lm$degree_rank_in_community >= 1))
  # ECG-only columns
  expect_true(all(c("pool_lead_freq", "confidence") %in% names(lm)))
  expect_true(all(lm$pool_lead_freq >= 0 & lm$pool_lead_freq <= 1))
  expect_true(all(lm$confidence >= 0 & lm$confidence <= 1))
})

test_that("baselines without leaders are scored the same way, and say so", {
  g <- zachary()
  cl <- igraph::cluster_louvain(g)
  m <- lcda_metrics(cl, g)
  val <- function(k) m$value[m$metric == k]

  expect_equal(val("Q"), igraph::modularity(cl), tolerance = 1e-8)
  expect_equal(val("leaders_derived"), 1)     # flagged as derived, not joint
  expect_equal(val("n_leaders"), val("n_communities"))
  lm <- lcda_metrics(cl, g, level = "leader")
  expect_true(all(grepl("derived", lm$source)))

  # a bare membership vector is accepted too
  m2 <- lcda_metrics(as.integer(igraph::membership(cl)), g)
  expect_equal(m2$value[m2$metric == "Q"], val("Q"), tolerance = 1e-10)
})

test_that("weighted graphs are handled by the same code path", {
  set.seed(11)
  g <- igraph::make_graph("Zachary")
  igraph::E(g)$weight <- stats::runif(igraph::ecount(g), 0.5, 2)
  res <- lcda_grasp(g, B = 10, seed = 1)
  m <- lcda_metrics(res)
  cm <- lcda_metrics(res, level = "community")

  expect_equal(m$value[m$metric == "Q"],
               igraph::modularity(g, res$best$membership, weights = igraph::E(g)$weight),
               tolerance = 1e-8)
  expect_equal(sum(cm$q_contribution), res$best$Q, tolerance = 1e-8)
  expect_equal(m$value[m$metric == "total_edge_weight"],
               sum(igraph::E(g)$weight), tolerance = 1e-8)
})

test_that("lcda_metrics covers the quantities the paper's benchmark tables report", {
  # `repro_benchmarks` is the per-run table behind the paper's benchmark
  # tables. Every one of its measured columns must be obtainable from a fresh
  # run through lcda_metrics(), otherwise a user cannot rebuild those tables.
  skip_if_not(nzchar(system.file("extdata", "repro_benchmarks.rds",
                                package = "lcdaGRASP")))
  d <- lcda_data("repro_benchmarks")
  skip_if(is.null(d))

  paper_to_metric <- c(Q = "Q", H = "H_global", g_comms = "n_communities",
                       secs = "elapsed_sec")
  expect_true(all(names(paper_to_metric) %in% names(d$results)))

  res <- lcda_grasp(zachary(), B = 10, seed = 1)
  m <- lcda_metrics(res)
  expect_true(all(paper_to_metric %in% m$metric))

  # and the recovery columns used by the LFR / real-network tables
  m2 <- lcda_metrics(res, truth = c(rep(1L, 17), rep(2L, 17)))
  expect_true(all(c("nmi", "ari") %in% m2$metric))
})

test_that("the participation coefficient is bounded and zero on a clique partition", {
  g <- zachary()
  res <- lcda_grasp(g, B = 10, seed = 1)
  lm <- lcda_metrics(res, level = "leader")
  expect_true(all(lm$participation >= 0 & lm$participation < 1))

  # everything in one community => no edge leaves it => P = 0 for every node
  one <- lcda_metrics(list(membership = rep(1L, igraph::vcount(g))), g,
                      level = "leader")
  expect_equal(one$participation, 0, tolerance = 1e-12)

  # all singletons => every edge leaves => P = 1 - 1/k for a degree-k node
  sing <- lcda_metrics(list(membership = seq_len(igraph::vcount(g))), g,
                       level = "leader")
  k <- as.numeric(igraph::degree(g))[sing$leader]
  expect_equal(sing$participation, 1 - 1 / k, tolerance = 1e-10)
})

test_that("node_score reproduces the paper's external leader validation", {
  g <- zachary()
  res <- lcda_grasp(g, B = 10, seed = 1)

  # scoring leaders by their own within-community degree must put them at the
  # top of their community, by construction of the repair step
  lm <- lcda_metrics(res, level = "leader", node_score = NULL)
  wd <- rep(0, igraph::vcount(g))
  el <- igraph::as_edgelist(g, names = FALSE)
  mem <- res$best$membership
  intra <- el[mem[el[, 1]] == mem[el[, 2]], , drop = FALSE]
  tb <- table(c(intra))
  wd[as.integer(names(tb))] <- as.numeric(tb)
  scored <- lcda_metrics(res, level = "leader", node_score = wd)
  expect_equal(scored$score, lm$degree_within)
  expect_equal(scored$score_pctile_in_community, lm$degree_pctile_in_community)

  # a constant score puts every leader at percentile 1 and rank 1
  const <- lcda_metrics(res, level = "leader",
                        node_score = rep(1, igraph::vcount(g)))
  expect_true(all(const$score_pctile_in_community == 1))
  expect_true(all(const$score_rank_in_community == 1L))

  # the summary scope mirrors the per-leader columns
  ov <- lcda_metrics(res, node_score = wd)
  ls <- ov[ov$scope == "leader_score", ]
  expect_equal(ls$value[ls$metric == "score_pctile_mean"],
               mean(scored$score_pctile_in_community))
  expect_true(all(ls$value >= 0 & ls$value <= 1))
  expect_false("leader_score" %in% lcda_metrics(res)$scope)
  expect_error(lcda_metrics(res, node_score = 1:3), "length")
})

test_that("lcda_metrics rejects malformed input", {
  g <- zachary()
  expect_error(lcda_metrics("nonsense", g), "solution")
  expect_error(lcda_metrics(rep(1L, 5), g), "length")
  expect_error(lcda_metrics(rep(1L, 34)), "graph")
  expect_error(lcda_metrics(list(membership = rep(1L, 34), leaders = 999L), g),
               "1-based")
})
