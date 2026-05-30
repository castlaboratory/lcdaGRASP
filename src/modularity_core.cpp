// modularity_core.cpp
// Performance-critical kernels for modularity, incremental delta-Q,
// community affinities, first-improvement local search, and the VNMI
// extension for n > 300. Mirrors equations (1), (11), (12) of
// Ospina et al. (2026) "A GRASP Framework for Community and Leader
// Detection in Complex Networks".
//
// Conventions:
//   - The graph is passed as a sparse adjacency representation:
//     `adj_indptr` (length n+1) and `adj_indices` (length 2m for an
//     undirected simple graph stored symmetrically). This is the same
//     CSR layout produced by Matrix::sparseMatrix.
//   - Communities are passed as an IntegerVector `membership` of length n,
//     values 0..d-1.
//   - All indexing is 0-based at the C++ level; the R wrappers translate.
//   - Memory layout is deliberately cache-friendly: affinities are kept
//     as a dense n x d matrix because d is typically small (< 50) and
//     access is by row in the hot loop.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

// ---------- helpers ----------------------------------------------------

// Returns the degree vector from CSR.
static inline std::vector<int> degree_from_csr(const IntegerVector& indptr) {
  int n = indptr.size() - 1;
  std::vector<int> deg(n);
  for (int i = 0; i < n; ++i) deg[i] = indptr[i+1] - indptr[i];
  return deg;
}

// Returns 2m (sum of degrees in undirected simple graph stored symmetrically).
static inline double two_m_from_degrees(const std::vector<int>& deg) {
  double tm = 0.0;
  for (int d : deg) tm += d;
  return tm;
}

// Weighted node strength s_i = sum of incident edge weights (the weighted
// analogue of degree). For an all-ones weight vector this reduces exactly to
// the integer degree, so the unweighted kernels are the special case w = 1.
static inline std::vector<double> strength_from_csr(const IntegerVector& indptr,
                                                    const NumericVector& weights) {
  int n = indptr.size() - 1;
  std::vector<double> s(n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int p = indptr[i]; p < indptr[i+1]; ++p) s[i] += weights[p];
  return s;
}

// = 2W (each undirected edge contributes its weight to both endpoints).
static inline double two_W_from_strength(const std::vector<double>& s) {
  double tw = 0.0;
  for (double v : s) tw += v;
  return tw;
}

// ---------- (1) full modularity ---------------------------------------

// Full modularity Q (eq. 1, paper)
//
// @param adj_indptr CSR row pointer (length n+1)
// @param adj_indices CSR neighbour indices (length 2m), unweighted
// @param membership integer vector, length n, values in [0, d-1]
// @return double scalar
// @param adj_weights CSR edge weights (length 2m); all-ones for unweighted.
// [[Rcpp::export]]
double modularity_cpp(IntegerVector adj_indptr,
                      IntegerVector adj_indices,
                      IntegerVector membership,
                      NumericVector adj_weights) {
  const int n = adj_indptr.size() - 1;
  if (membership.size() != n) stop("membership length must equal n");

  std::vector<double> deg = strength_from_csr(adj_indptr, adj_weights);  // weighted strength
  const double two_m = two_W_from_strength(deg);                          // 2W
  if (two_m <= 0.0) return 0.0;

  // Q = (1/2W) sum_ij [ w_ij - s_i s_j / 2W ] delta(c_i, c_j)
  //   = sum_c [ L_c / 2W  -  ( D_c / 2W )^2 ]
  // where m == W is the total edge weight (two_m == 2W), L_c is the intra-
  // community weight counted BOTH orientations (= 2 * e_in(c)), so the intra
  // term is L_c/2W (NOT L_c/W), and D_c is the community's total strength.
  // Self-loops are absent in simple graphs; w_ii = 0. Reduces to the
  // unweighted formula when all weights are 1.

  // discover d
  int d = 0;
  for (int v = 0; v < n; ++v) if (membership[v] + 1 > d) d = membership[v] + 1;

  std::vector<double> L(d, 0.0);   // sum of a_ij over (i,j) in same community (both orientations)
  std::vector<double> D(d, 0.0);   // sum of k_i over i in community

  for (int i = 0; i < n; ++i) {
    int ci = membership[i];
    D[ci] += deg[i];
    for (int p = adj_indptr[i]; p < adj_indptr[i+1]; ++p) {
      int j = adj_indices[p];
      if (membership[j] == ci) L[ci] += adj_weights[p];  // weighted, both orientations
    }
  }

  double Q = 0.0;
  for (int c = 0; c < d; ++c) {
    // L[c] above is 2 * (intra edges)
    double l_c = L[c];                    // = 2 * e_in(c)
    double d_c = D[c];
    Q += l_c / two_m - (d_c / two_m) * (d_c / two_m);
  }
  return Q;
}

// ---------- (11) community affinities --------------------------------

// Community affinity matrix s_i^{(t)} = sum_{l in C_t} b_il (eq. 11)
//
// Returns an n x d numeric matrix S where S(i, t) is the surprise sum
// between node i and community t under the configuration null model.
//
// Time complexity: O(2m + n*d). The 2m comes from the sparse edge scan;
// the n*d comes from the degree-product correction term.
//
// [[Rcpp::export]]
NumericMatrix community_affinity_cpp(IntegerVector adj_indptr,
                                     IntegerVector adj_indices,
                                     IntegerVector membership,
                                     int d,
                                     NumericVector adj_weights) {
  const int n = adj_indptr.size() - 1;
  std::vector<double> deg = strength_from_csr(adj_indptr, adj_weights);  // strength
  const double two_m = two_W_from_strength(deg);                          // 2W

  // First pass: sum of w_il over l in community t for each i.
  NumericMatrix S(n, d);

  // Per-community total strength (for the null-model term).
  std::vector<double> Dc(d, 0.0);
  for (int v = 0; v < n; ++v) Dc[membership[v]] += deg[v];

  for (int i = 0; i < n; ++i) {
    // Edge contribution: a_il = w_il for neighbours.
    for (int p = adj_indptr[i]; p < adj_indptr[i+1]; ++p) {
      int l = adj_indices[p];
      S(i, membership[l]) += adj_weights[p];
    }
    // Subtract s_i s_l / 2W, summed over l in community t = s_i * Dc[t] / 2W.
    for (int t = 0; t < d; ++t) {
      S(i, t) -= deg[i] * Dc[t] / two_m;
    }
  }
  return S;
}

// ---------- (12) incremental delta-Q ---------------------------------

// Modularity change from moving node i to community t (eq. 12).
//
// delta Q = (1/W) * [ S(i, t)  -  ( S(i, c_i) - b_ii ) ]   (m == W below)
// with b_ii = -s_i^2 / (2W) in a simple graph (w_ii=0). The factor 2 from the
// symmetric pair (i,l)/(l,i) cancels the 1/(2W), leaving the 1/W prefactor.
//
// @param S precomputed affinity matrix from community_affinity_cpp
// @param i 0-based node index
// @param t 0-based candidate community index (must differ from current)
// @param current_comm 0-based current community of i
// @param deg degree vector
// @param two_m sum of degrees
// [[Rcpp::export]]
double delta_modularity_cpp(NumericMatrix S,
                            int i, int t, int current_comm,
                            NumericVector deg, double two_m) {
  if (t == current_comm) return 0.0;
  const double m = two_m / 2.0;                       // = W
  const double k_i = deg[i];                          // strength s_i
  const double b_ii = - (k_i * k_i) / two_m;          // a_ii = 0
  const double affinity_to_dest = S(i, t);
  const double affinity_to_curr_minus_self = S(i, current_comm) - b_ii;
  return (affinity_to_dest - affinity_to_curr_minus_self) / m;
}

// ---------- update affinities after a move ----------------------------

// After node i moves from c_from to c_to, only S(l, c_from) and S(l, c_to)
// change for l a neighbour of i (plus a global correction in the null
// term whenever a community's total degree changes). We bundle the
// degree-shift into the correction.
//
// Returns the updated S in-place (matrix passed by reference).
// Caller must also update the membership vector.
//
// This is *the* hot function in Algorithm 3 of the paper.

// [[Rcpp::export]]
void update_affinities_inplace_cpp(NumericMatrix S,
                                   IntegerVector adj_indptr,
                                   IntegerVector adj_indices,
                                   int i, int c_from, int c_to,
                                   NumericVector deg, double two_m,
                                   NumericVector adj_weights) {
  const int n = adj_indptr.size() - 1;
  const double k_i = deg[i];   // strength s_i

  // 1) Edge part: for each neighbour l of i, l's affinity to c_from loses the
  //    edge weight w_il and c_to gains it.
  for (int p = adj_indptr[i]; p < adj_indptr[i+1]; ++p) {
    int l = adj_indices[p];
    S(l, c_from) -= adj_weights[p];
    S(l, c_to)   += adj_weights[p];
  }

  // 2) Null part: D_{c_from} decreases by k_i, D_{c_to} increases by k_i.
  //    Therefore for *every* node l, the null term k_l * D_{c}/2m shifts:
  //      S(l, c_from) is *minus* the null term => add k_l * k_i / 2m
  //      S(l, c_to)   gets -k_l * k_i / 2m
  // (S is defined as edge-count - null-term)
  const double k_i_over_2m = k_i / two_m;
  for (int l = 0; l < n; ++l) {
    S(l, c_from) += deg[l] * k_i_over_2m;
    S(l, c_to)   -= deg[l] * k_i_over_2m;
  }
}

// Note: the global O(n) sweep above is the price we pay for using a dense
// S matrix. For very large n a sparse + lazy strategy is preferable. The
// VNMI variant below restricts updates to a sampled subset of size n'.

// ---------- (Algorithm 3) first-improvement local search -------------

// First-improvement local search on modularity.
//
// Implements Algorithm 3 of the paper. Iterates until no improving move
// is found. Returns the updated membership vector and the final Q.
//
// For n > n_threshold the caller should dispatch to VNMI instead.
// [[Rcpp::export]]
List local_search_cpp(IntegerVector adj_indptr,
                      IntegerVector adj_indices,
                      IntegerVector membership_in,
                      int d,
                      NumericVector adj_weights,
                      double tol = 1e-12,
                      int max_iter = 1000000) {
  const int n = adj_indptr.size() - 1;
  IntegerVector membership = clone(membership_in);
  std::vector<double> deg = strength_from_csr(adj_indptr, adj_weights);
  const double two_m = two_W_from_strength(deg);

  NumericMatrix S = community_affinity_cpp(adj_indptr, adj_indices, membership, d, adj_weights);
  NumericVector deg_iv = NumericVector(deg.begin(), deg.end());

  int moves = 0;
  bool improved = true;
  while (improved && moves < max_iter) {
    improved = false;
    for (int i = 0; i < n; ++i) {
      int ci = membership[i];
      double best_gain = tol;
      int best_t = -1;
      for (int t = 0; t < d; ++t) {
        if (t == ci) continue;
        double g = delta_modularity_cpp(S, i, t, ci, deg_iv, two_m);
        if (g > best_gain) { best_gain = g; best_t = t; }
      }
      if (best_t >= 0) {
        update_affinities_inplace_cpp(S, adj_indptr, adj_indices,
                                      i, ci, best_t, deg_iv, two_m, adj_weights);
        membership[i] = best_t;
        moves++;
        improved = true;
        break;   // restart outer loop (first-improvement semantics)
      }
    }
  }
  double Q = modularity_cpp(adj_indptr, adj_indices, membership, adj_weights);
  return List::create(
    _["membership"] = membership,
    _["modularity"] = Q,
    _["moves"] = moves
  );
}

// ---------- VNMI: multi-improvement on a sampled subset -------------

// Variable Neighbourhood Multi-Improvement local search (n > n').
//
// Restricts evaluation to a random subset V' of size `n_prime`, and
// accepts *all* profitable moves in a single pass before recomputing
// affinities. Repeats passes until no profitable move exceeds `eps`.
// The V' sample uses R's RNG, so reproducibility is controlled by set.seed()
// at the R level (no separate seed argument).
// [[Rcpp::export]]
List vnmi_local_search_cpp(IntegerVector adj_indptr,
                           IntegerVector adj_indices,
                           IntegerVector membership_in,
                           int d,
                           NumericVector adj_weights,
                           int n_prime = 300,
                           double eps = 1e-4,
                           int max_passes = 200) {
  const int n = adj_indptr.size() - 1;
  IntegerVector membership = clone(membership_in);
  std::vector<double> deg = strength_from_csr(adj_indptr, adj_weights);
  const double two_m = two_W_from_strength(deg);

  if (n_prime > n) n_prime = n;

  // Use R's RNG to get a reproducible sample.
  GetRNGstate();
  std::vector<int> V_prime(n_prime);
  {
    // sample without replacement via reservoir
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    for (int k = 0; k < n_prime; ++k) {
      double u = unif_rand();
      int swap = k + (int)std::floor(u * (n - k));
      if (swap >= n) swap = n - 1;
      std::swap(idx[k], idx[swap]);
      V_prime[k] = idx[k];
    }
  }
  PutRNGstate();

  NumericVector deg_iv = NumericVector(deg.begin(), deg.end());
  int passes = 0;
  while (passes < max_passes) {
    NumericMatrix S = community_affinity_cpp(adj_indptr, adj_indices, membership, d, adj_weights);

    // For each node v' in V', determine its best target community.
    std::vector<int>     proposed_to(n_prime, -1);
    std::vector<double>  proposed_gain(n_prime, 0.0);
    double max_gain = 0.0;
    for (int k = 0; k < n_prime; ++k) {
      int v = V_prime[k];
      int cv = membership[v];
      double best_gain = eps;
      int best_t = -1;
      for (int t = 0; t < d; ++t) {
        if (t == cv) continue;
        double g = delta_modularity_cpp(S, v, t, cv, deg_iv, two_m);
        if (g > best_gain) { best_gain = g; best_t = t; }
      }
      if (best_t >= 0) {
        proposed_to[k]   = best_t;
        proposed_gain[k] = best_gain;
        if (best_gain > max_gain) max_gain = best_gain;
      }
    }
    if (max_gain < eps) break;

    // Apply all proposed moves at once.
    // (Concurrent moves may interact, but the multi-improvement strategy
    //  is justified empirically by Rios et al. 2017; the affinity matrix
    //  is recomputed from scratch in the next pass.)
    for (int k = 0; k < n_prime; ++k) {
      if (proposed_to[k] >= 0) membership[V_prime[k]] = proposed_to[k];
    }
    passes++;
  }
  double Q = modularity_cpp(adj_indptr, adj_indices, membership, adj_weights);
  return List::create(
    _["membership"] = membership,
    _["modularity"] = Q,
    _["passes"] = passes,
    _["n_prime"] = n_prime
  );
}
