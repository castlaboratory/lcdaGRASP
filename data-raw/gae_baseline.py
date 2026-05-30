#!/usr/bin/env python3
"""Graph Auto-Encoder (GAE; Kipf & Welling 2016) community-detection baseline.

The paper argues GNN methods are not directly applicable to attribute-free
graphs. To make the comparison honest rather than deflected, we run a canonical
GAE on the SAME topology, supplying *synthesised structural features* (degree,
clustering, k-core, eigenvector centrality) as node attributes -- the standard
way to apply GNNs to plain graphs. A 2-layer GCN encoder reconstructs the
adjacency; KMeans on the learned embedding yields communities.

Implemented in plain PyTorch (the graph convolution is a sparse-ish matmul), so
it needs only torch + scikit-learn, not the full PyG stack.

Usage: gae_baseline.py --edges E --n N --k K --out OUT [--dim 16 --epochs 200 --seed 1]
Writes one 1-based community label per node (node order 1..n).
"""
import argparse
import numpy as np
import networkx as nx
import torch
import torch.nn as nn
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def structural_features(G, n):
    deg = dict(G.degree())
    clu = nx.clustering(G)
    core = nx.core_number(G)
    try:
        eig = nx.eigenvector_centrality_numpy(G)
    except Exception:
        eig = {v: 0.0 for v in G.nodes()}
    X = np.zeros((n, 4), dtype=np.float32)
    for v in G.nodes():
        X[v] = [deg.get(v, 0), clu.get(v, 0.0), core.get(v, 0), eig.get(v, 0.0)]
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    return X


def norm_adj(A):
    A = A + np.eye(A.shape[0], dtype=np.float32)
    d = A.sum(1)
    dinv = np.power(d, -0.5, where=d > 0)
    dinv[d == 0] = 0.0
    return (dinv[:, None] * A) * dinv[None, :]


class GAE(nn.Module):
    def __init__(self, in_dim, hid=32, z=16):
        super().__init__()
        self.w0 = nn.Linear(in_dim, hid)
        self.w1 = nn.Linear(hid, z)

    def encode(self, Ahat, X):
        h = torch.relu(Ahat @ self.w0(X))
        return Ahat @ self.w1(h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dim", type=int, default=16)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    torch.manual_seed(a.seed); np.random.seed(a.seed)

    G = nx.empty_graph(a.n)
    with open(a.edges) as fe:
        for line in fe:
            p = line.split()
            if len(p) >= 2:
                u, v = int(p[0]) - 1, int(p[1]) - 1
                if u != v:
                    G.add_edge(u, v)
    A = nx.to_numpy_array(G, nodelist=range(a.n), dtype=np.float32)
    X = torch.tensor(structural_features(G, a.n))
    Ahat = torch.tensor(norm_adj(A))
    At = torch.tensor(A)
    pos_w = torch.tensor((A.size - A.sum()) / (A.sum() + 1e-8))

    model = GAE(X.shape[1], hid=2 * a.dim, z=a.dim)
    opt = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pos_w)
    for _ in range(a.epochs):
        opt.zero_grad()
        Z = model.encode(Ahat, X)
        logits = Z @ Z.t()
        loss = lossf(logits, At)
        loss.backward(); opt.step()

    Z = model.encode(Ahat, X).detach().numpy()
    if a.k > 0:
        # ground-truth k supplied (generous baseline)
        lab = KMeans(n_clusters=a.k, n_init=10, random_state=a.seed).fit_predict(Z)
    else:
        # auto-k: no ground-truth count -- pick k maximising the silhouette
        # over a candidate range (the realistic attribute-free setting).
        kmax = max(2, min(40, a.n // 5))
        best_lab, best_s = None, -2.0
        for kk in range(2, kmax + 1):
            l = KMeans(n_clusters=kk, n_init=10, random_state=a.seed).fit_predict(Z)
            if len(set(l)) < 2:
                continue
            s = silhouette_score(Z, l)
            if s > best_s:
                best_s, best_lab = s, l
        lab = best_lab if best_lab is not None else np.zeros(a.n, dtype=int)
    with open(a.out, "w") as fo:
        for l in lab:
            fo.write(f"{int(l) + 1}\n")


if __name__ == "__main__":
    main()
