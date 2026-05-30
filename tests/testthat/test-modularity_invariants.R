# tests/testthat/test-modularity_invariants.R
#
# These tests are *not* cosmetic. They check the two invariants whose
# violation would silently destroy the algorithm:
#
#   T1) Incremental delta-Q exactly equals modularity recomputed from
#       scratch after the move. If this fails, the whole local-search
#       phase optimises a function different from Q.
#
#   T2) The affinity update keeps S in sync with a from-scratch
#       community_affinity_cpp recomputation after each move. If this
#       fails, subsequent delta-Q evaluations are corrupted.
#
#   T3) Lemma 2(ii) — |RCL(alpha)| is non-decreasing in alpha.
#
#   T4) NCE bounded in [0, 1] and reaches 1 only when p = 1/2 (P5).

library(testthat)
library(igraph)

# -----------------------------------------------------------------
.toy_graph <- function(seed = 1, n = 50, p = 0.15) {
  set.seed(seed)
  igraph::sample_gnp(n, p, directed = FALSE)
}
.as_csr_local <- lcdaGRASP:::.as_csr   # internal helper

test_that("T1: incremental delta-Q equals from-scratch delta-Q", {
  g <- .toy_graph()
  csr <- .as_csr_local(g)
  n <- csr$n
  # arbitrary partition into 4 communities
  set.seed(7)
  mem <- sample.int(4, n, replace = TRUE) - 1L

  S  <- community_affinity_cpp(csr$indptr, csr$indices, mem, 4L, csr$weights)
  deg <- as.numeric(igraph::degree(g))   # strength == degree on an unweighted graph
  two_m <- sum(deg)

  Q0 <- modularity_cpp(csr$indptr, csr$indices, mem, csr$weights)
  for (i in sample.int(n, 10)) {
    ci <- mem[i]
    candidates <- setdiff(0:3, ci)
    for (t in candidates) {
      dq <- delta_modularity_cpp(S, i - 1L, t, ci, deg, two_m)
      mem_new <- mem; mem_new[i] <- t
      Q1 <- modularity_cpp(csr$indptr, csr$indices, mem_new, csr$weights)
      expect_equal(Q1 - Q0, dq, tolerance = 1e-10,
                   info = sprintf("i=%d ci=%d t=%d", i, ci, t))
    }
  }
})

test_that("T2: affinity update is consistent after a move", {
  g <- .toy_graph(seed = 3)
  csr <- .as_csr_local(g); n <- csr$n
  set.seed(11)
  mem <- sample.int(3, n, replace = TRUE) - 1L
  S  <- community_affinity_cpp(csr$indptr, csr$indices, mem, 3L, csr$weights)
  deg <- as.numeric(igraph::degree(g))   # strength == degree on an unweighted graph
  two_m <- sum(deg)

  i <- 5L; c_from <- mem[i + 1L]; c_to <- (c_from + 1L) %% 3L
  update_affinities_inplace_cpp(S, csr$indptr, csr$indices,
                                i, c_from, c_to, deg, two_m, csr$weights)
  mem[i + 1L] <- c_to
  S_ref <- community_affinity_cpp(csr$indptr, csr$indices, mem, 3L, csr$weights)
  expect_lt(max(abs(S - S_ref)), 1e-10)
})

test_that("T3: |RCL(alpha)| is non-decreasing in alpha (Lemma 2)", {
  set.seed(1)
  scores <- runif(100)
  alphas <- seq(0, 1, by = 0.05)
  sizes  <- sapply(alphas, function(a) length(rcl_max_cpp(scores, a)))
  expect_true(all(diff(sizes) >= 0))
  expect_equal(sizes[1], sum(scores == max(scores)))
  expect_equal(sizes[length(sizes)], length(scores))
})

test_that("T4: NCE bounded in [0,1] and maximised at p=1/2 (Prop 5)", {
  n <- 1000
  # Sweep degrees from 0 to n-1 and check H_b shape
  ks <- seq(0, n - 1, by = 50)
  vals <- vapply(ks, function(k) nce_node_cpp(k, n), numeric(1))
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1 + 1e-12))
  # max should be near k = n/2
  k_max <- ks[which.max(vals)]
  expect_lt(abs(k_max - n/2), n/10)
})

test_that("T5: full local search increases modularity monotonically", {
  g <- .toy_graph(seed = 42, n = 80)
  csr <- .as_csr_local(g); n <- csr$n
  set.seed(123)
  mem0 <- sample.int(5, n, replace = TRUE) - 1L
  Q_before <- modularity_cpp(csr$indptr, csr$indices, mem0, csr$weights)
  res <- local_search_cpp(csr$indptr, csr$indices, mem0, 5L, csr$weights)
  expect_gte(res$modularity, Q_before - 1e-10)
})

test_that("T6: end-to-end LCDA-GRASP runs and returns valid partition", {
  set.seed(2026)
  # Karate-club analogue
  g <- igraph::make_graph("Zachary")
  res <- lcda_grasp(g, alpha_c = 0.1, alpha_s = 0.3, B = 5,
                    centrality = "eigen", similarity = "hpi")
  expect_s3_class(res, "lcda_grasp_result")
  expect_equal(length(res$best$membership), igraph::vcount(g))
  expect_true(all(res$best$leaders %in% seq_len(igraph::vcount(g))))
  expect_true(res$best$Q > 0)
})

test_that("T7: LCDA-ECG returns a valid consensus partition + one leader per community", {
  set.seed(2026)
  g <- igraph::make_graph("Zachary")
  res <- lcda_ecg(g, B = 16, seed = 7)
  expect_s3_class(res, "lcda_ecg_result")
  expect_equal(length(res$membership), igraph::vcount(g))
  # exactly one leader per consensus community, all valid 1-based vertices
  expect_equal(length(res$leaders), length(unique(res$membership)))
  expect_true(all(res$leaders %in% seq_len(igraph::vcount(g))))
  expect_true(all(res$leaders_central %in% seq_len(igraph::vcount(g))))
  # each consensus leader lies inside its own community
  comms <- sort(unique(res$membership))
  expect_true(all(res$membership[res$leaders] == comms))
  # confidence is a per-node value bounded in [0, 1]
  expect_equal(length(res$confidence), igraph::vcount(g))
  expect_true(all(res$confidence >= 0 & res$confidence <= 1 + 1e-9))
  expect_true(res$Q > 0)
})

test_that("T8: LCDA-ECG overlap returns valid overlapping memberships", {
  set.seed(2026)
  g <- igraph::make_graph("Zachary")
  res <- lcda_ecg(g, B = 16, overlap = TRUE, tau = 0.7, seed = 7)
  n <- igraph::vcount(g)
  expect_equal(length(res$overlap_membership), n)
  expect_equal(length(res$is_overlap), n)
  # every node belongs to at least its hard community; all ids valid
  expect_true(all(vapply(res$overlap_membership, length, integer(1)) >= 1))
  expect_true(all(unlist(res$overlap_membership) %in% sort(unique(res$membership))))
  expect_true(all(mapply(function(s, h) h %in% s, res$overlap_membership, res$membership)))
  # overlap flag matches multi-membership
  expect_equal(res$is_overlap, vapply(res$overlap_membership, length, integer(1)) > 1L)
})

test_that("T9: weighted modularity matches igraph (and unweighted is the w=1 case)", {
  set.seed(5)
  g <- igraph::sample_gnp(70, 0.12, directed = FALSE)
  igraph::E(g)$weight <- runif(igraph::ecount(g), 0.3, 4)
  csr <- lcdaGRASP:::.as_csr(g)
  cl  <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  mem <- as.integer(igraph::membership(cl)) - 1L
  q_cpp <- modularity_cpp(csr$indptr, csr$indices, mem, csr$weights)
  q_ig  <- igraph::modularity(g, igraph::membership(cl), weights = igraph::E(g)$weight)
  expect_equal(q_cpp, q_ig, tolerance = 1e-8)
  # end-to-end LCDA on a weighted graph runs and optimises weighted Q
  res <- lcda_grasp(g, B = 15, seed = 1)
  expect_s3_class(res, "lcda_grasp_result")
  expect_true(res$best$Q > 0)
})
