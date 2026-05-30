// similarity_and_quality.cpp
//
// Implements:
//   - HPI, Dice, Jaccard similarity (eqs. 5-7) vectorised against a fixed
//     leader L: returns sim(v, L) for all v in the unclassified pool.
//   - NCE leader-quality score (eqs. 13-14), both the global definition
//     of the paper AND a community-conditioned alternative I_loc that
//     addresses the conceptual gap discussed in the critical review.
//   - Eigenvector centrality via power iteration, O(m log(n/eps)).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// ---------- common neighbourhood -------------------------------------

// Returns |Γ(u) ∩ Γ(v)| using a sorted-merge over the CSR neighbour lists.
// Both lists are assumed sorted ascending (the caller is responsible).
static int common_neighbours(int u, int v,
                             const IntegerVector& indptr,
                             const IntegerVector& indices) {
  int pu = indptr[u], pu_end = indptr[u+1];
  int pv = indptr[v], pv_end = indptr[v+1];
  int c = 0;
  while (pu < pu_end && pv < pv_end) {
    int a = indices[pu], b = indices[pv];
    if      (a < b) ++pu;
    else if (a > b) ++pv;
    else { ++c; ++pu; ++pv; }
  }
  return c;
}

// ---------- similarity vs leader L -----------------------------------

// HPI similarity (eq. 5) of every node in `pool` against leader L.
// [[Rcpp::export]]
NumericVector similarity_hpi_cpp(IntegerVector adj_indptr,
                                 IntegerVector adj_indices,
                                 int L,
                                 IntegerVector pool) {
  int kL = adj_indptr[L+1] - adj_indptr[L];
  int n_pool = pool.size();
  NumericVector out(n_pool);
  for (int i = 0; i < n_pool; ++i) {
    int v = pool[i];
    if (v == L) { out[i] = 1.0; continue; }
    int kv = adj_indptr[v+1] - adj_indptr[v];
    int c  = common_neighbours(L, v, adj_indptr, adj_indices);
    int denom = std::min(kL, kv);
    out[i] = (denom == 0) ? 0.0 : (double)c / (double)denom;
  }
  return out;
}

// Dice similarity (eq. 7)
// [[Rcpp::export]]
NumericVector similarity_dice_cpp(IntegerVector adj_indptr,
                                  IntegerVector adj_indices,
                                  int L,
                                  IntegerVector pool) {
  int kL = adj_indptr[L+1] - adj_indptr[L];
  int n_pool = pool.size();
  NumericVector out(n_pool);
  for (int i = 0; i < n_pool; ++i) {
    int v = pool[i];
    if (v == L) { out[i] = 1.0; continue; }
    int kv = adj_indptr[v+1] - adj_indptr[v];
    int c  = common_neighbours(L, v, adj_indptr, adj_indices);
    int denom = kL + kv;
    out[i] = (denom == 0) ? 0.0 : 2.0 * c / (double)denom;
  }
  return out;
}

// Jaccard similarity (eq. 6)
// [[Rcpp::export]]
NumericVector similarity_jaccard_cpp(IntegerVector adj_indptr,
                                     IntegerVector adj_indices,
                                     int L,
                                     IntegerVector pool) {
  int kL = adj_indptr[L+1] - adj_indptr[L];
  int n_pool = pool.size();
  NumericVector out(n_pool);
  for (int i = 0; i < n_pool; ++i) {
    int v = pool[i];
    if (v == L) { out[i] = 1.0; continue; }
    int kv = adj_indptr[v+1] - adj_indptr[v];
    int c  = common_neighbours(L, v, adj_indptr, adj_indices);
    int denom = kL + kv - c;
    out[i] = (denom == 0) ? 0.0 : (double)c / (double)denom;
  }
  return out;
}

// ---------- NCE leader score -----------------------------------------

// Binary entropy in bits.
static inline double Hb(double p) {
  if (p <= 0.0 || p >= 1.0) return 0.0;
  return -(p * std::log2(p) + (1.0 - p) * std::log2(1.0 - p));
}

// Node-Connection Entropy I(v) — global definition (eq. 13).
//
// I(v) = H_b(k_v / n) if 0 < k_v < n; else (1 + log2 n)/n.
//
// Note (see critical review): in a simple graph k_v <= n-1, so p1 = 1
// is unreachable; the fallback fires only for isolated nodes (k_v = 0).
// [[Rcpp::export]]
double nce_node_cpp(int kv, int n) {
  double p = (double)kv / (double)n;
  if (p > 0.0 && p < 1.0) return Hb(p);
  return (1.0 + std::log2((double)n)) / (double)n;
}

// Global NCE leader-quality score H(C) (eq. 14).
//
// @param leader_degrees integer vector of leader degrees, length g
// @param n total number of nodes
// [[Rcpp::export]]
double nce_score_cpp(IntegerVector leader_degrees, int n) {
  int g = leader_degrees.size();
  if (g == 0) return 0.0;
  double s = 0.0;
  for (int j = 0; j < g; ++j) s += nce_node_cpp(leader_degrees[j], n);
  return s / (double)g;
}

// Community-conditioned NCE — proposed alternative (see review §4).
//
// I_loc(v) = H_b( |Γ(v) ∩ C(v)| / (|C(v)| - 1) )
// i.e., entropy of the leader's connection probability *within its own
// community*. Penalises both hubs that span everywhere and isolated
// leaders inside their group.
// [[Rcpp::export]]
double nce_local_cpp(IntegerVector adj_indptr,
                     IntegerVector adj_indices,
                     IntegerVector membership,
                     IntegerVector leaders) {
  const int g = leaders.size();
  if (g == 0) return 0.0;

  // community sizes
  int d = 0;
  for (int v = 0; v < membership.size(); ++v)
    if (membership[v] + 1 > d) d = membership[v] + 1;
  std::vector<int> sz(d, 0);
  for (int v = 0; v < membership.size(); ++v) sz[membership[v]]++;

  double s = 0.0;
  for (int j = 0; j < g; ++j) {
    int L  = leaders[j];
    int cL = membership[L];
    if (sz[cL] <= 1) {
      s += 0.0;  // singleton community — no within-community structure
      continue;
    }
    int kin = 0;
    for (int p = adj_indptr[L]; p < adj_indptr[L+1]; ++p) {
      if (membership[adj_indices[p]] == cL) ++kin;
    }
    double p = (double)kin / (double)(sz[cL] - 1);
    if (p < 0.0) p = 0.0; if (p > 1.0) p = 1.0;
    s += Hb(p);
  }
  return s / (double)g;
}

// ---------- Eigenvector centrality via power iteration --------------

// Eigenvector centrality (eq. 2), power iteration on the adjacency.
//
// Returns the unit-L1-normalised positive eigenvector corresponding to
// the spectral radius. For disconnected graphs the result corresponds
// to the largest connected component's eigenvector embedded in the
// full space; isolated nodes get score 0.
// [[Rcpp::export]]
NumericVector eigen_centrality_cpp(IntegerVector adj_indptr,
                                   IntegerVector adj_indices,
                                   int max_iter = 1000,
                                   double tol = 1e-10) {
  const int n = adj_indptr.size() - 1;
  NumericVector x(n, 1.0 / std::sqrt((double)n));
  NumericVector y(n);
  for (int it = 0; it < max_iter; ++it) {
    // y = A x
    double max_y = 0.0;
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int p = adj_indptr[i]; p < adj_indptr[i+1]; ++p)
        s += x[adj_indices[p]];
      y[i] = s;
      if (std::fabs(s) > max_y) max_y = std::fabs(s);
    }
    if (max_y == 0.0) return x;        // empty graph: stay at uniform
    // normalise to ||y||_inf = 1 to keep numbers bounded
    double diff = 0.0;
    for (int i = 0; i < n; ++i) {
      double yi = y[i] / max_y;
      diff += std::fabs(yi - x[i]);
      x[i] = yi;
    }
    if (diff < tol) break;
  }
  // L1-normalise to a probability-like score
  double s1 = 0.0;
  for (int i = 0; i < n; ++i) s1 += x[i];
  if (s1 > 0.0) for (int i = 0; i < n; ++i) x[i] /= s1;
  return x;
}

// ---------- RCL builders --------------------------------------------

// RCL of maximum-type: keep e with f(e) >= fmax - alpha (fmax - fmin)
// (eq. 8 / 9). Returns the 0-based indices into `scores`.
// [[Rcpp::export]]
IntegerVector rcl_max_cpp(NumericVector scores, double alpha) {
  if (alpha < 0.0) alpha = 0.0;
  if (alpha > 1.0) alpha = 1.0;
  int n = scores.size();
  if (n == 0) return IntegerVector();
  double fmax = scores[0], fmin = scores[0];
  for (int i = 1; i < n; ++i) {
    if (scores[i] > fmax) fmax = scores[i];
    if (scores[i] < fmin) fmin = scores[i];
  }
  double thr = fmax - alpha * (fmax - fmin);
  std::vector<int> rcl;
  rcl.reserve(n);
  for (int i = 0; i < n; ++i) if (scores[i] >= thr) rcl.push_back(i);
  return IntegerVector(rcl.begin(), rcl.end());
}

// RCL of minimum-threshold type for similarities (eq. 10).
// Keep v with g(v) >= gmin + alpha (gmax - gmin).
// [[Rcpp::export]]
IntegerVector rcl_min_cpp(NumericVector scores, double alpha) {
  if (alpha < 0.0) alpha = 0.0;
  if (alpha > 1.0) alpha = 1.0;
  int n = scores.size();
  if (n == 0) return IntegerVector();
  double gmax = scores[0], gmin = scores[0];
  for (int i = 1; i < n; ++i) {
    if (scores[i] > gmax) gmax = scores[i];
    if (scores[i] < gmin) gmin = scores[i];
  }
  double thr = gmin + alpha * (gmax - gmin);
  std::vector<int> rcl;
  rcl.reserve(n);
  for (int i = 0; i < n; ++i) if (scores[i] >= thr) rcl.push_back(i);
  return IntegerVector(rcl.begin(), rcl.end());
}
