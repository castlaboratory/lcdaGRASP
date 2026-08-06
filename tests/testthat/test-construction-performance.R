# Guards for the issue-#46 construction optimisations. Both rewrites promise
# BIT-IDENTICAL results, so these tests assert exact equality, not tolerance.

test_that("two-hop similarity kernels equal the brute-force intersection", {
  # The kernels count |Γ(L) ∩ Γ(v)| via a two-hop walk from L (O(deg²) per
  # leader) instead of a per-pair sorted merge. Same integers, same doubles.
  set.seed(4601)
  for (rep in 1:5) {
    g <- igraph::sample_gnp(60, 0.12)
    csr <- as_csr(g)
    nbr <- lapply(seq_len(csr$n), function(v) {
      lo <- csr$indptr[v]; hi <- csr$indptr[v + 1L]
      if (hi > lo) csr$indices[(lo + 1L):hi] else integer(0)
    })
    pool <- as.integer(sample(seq_len(csr$n), 40) - 1L)   # 0-based
    L <- sample(seq_len(csr$n), 1) - 1L
    kL <- length(nbr[[L + 1L]])
    common <- vapply(pool, function(v)
      length(intersect(nbr[[L + 1L]], nbr[[v + 1L]])), integer(1))
    kv <- vapply(pool, function(v) length(nbr[[v + 1L]]), integer(1))

    ref_hpi <- ifelse(pool == L, 1, ifelse(pmin(kL, kv) == 0, 0, common / pmin(kL, kv)))
    ref_dice <- ifelse(pool == L, 1, ifelse(kL + kv == 0, 0, 2 * common / (kL + kv)))
    ref_jac <- ifelse(pool == L, 1,
                      ifelse(kL + kv - common == 0, 0, common / (kL + kv - common)))

    expect_identical(similarity_hpi_cpp(csr$indptr, csr$indices, L, pool), ref_hpi)
    expect_identical(similarity_dice_cpp(csr$indptr, csr$indices, L, pool), ref_dice)
    expect_identical(similarity_jaccard_cpp(csr$indptr, csr$indices, L, pool), ref_jac)
  }
})

test_that("hoisted centrality gives the same construction as the internal one", {
  g <- igraph::sample_gnp(200, 0.05)
  csr <- as_csr(g)
  cent <- lcdaGRASP:::.dispatch_centrality("eigen")(csr)
  set.seed(77)
  a <- lcda_construct(csr, 0.1, 0.3, variant = 1)
  set.seed(77)
  b <- lcda_construct(csr, 0.1, 0.3, variant = 1, cent = cent)
  expect_identical(a, b)
  expect_error(lcda_construct(csr, 0.1, 0.3, variant = 1, cent = c(1, 2)),
               "length")
})
