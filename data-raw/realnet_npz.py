#!/usr/bin/env python3
"""Fetch a gnn-benchmark .npz graph (Shchur et al. 2018) and emit an undirected
edge list + a ground-truth label per node, in the format the R side consumes
(1-based "u v" per line; one label per line in node order). The .npz is cached.

Datasets (--dataset):
  amazon_computers  -> co-purchase, 13,752 nodes, 10 product categories
  coauthor_physics  -> co-authorship, 34,493 nodes, 5 research fields
  coauthor_cs       -> co-authorship, 18,333 nodes, 15 research fields

Usage:
  realnet_npz.py --dataset coauthor_physics --edges E.txt --labels L.txt [--cache P.npz]
"""
import argparse
import io
import os
import sys
import urllib.request

import numpy as np
import scipy.sparse as sp

BASE = "https://github.com/shchur/gnn-benchmark/raw/master/data/npz/"
URLS = {
    "amazon_computers": BASE + "amazon_electronics_computers.npz",
    "coauthor_physics": BASE + "ms_academic_phy.npz",
    "coauthor_cs":      BASE + "ms_academic_cs.npz",
}


def load_npz(url, cache):
    if cache and os.path.exists(cache):
        with open(cache, "rb") as f:
            raw = f.read()
    else:
        raw = urllib.request.urlopen(url, timeout=120).read()
        if cache:
            with open(cache, "wb") as f:
                f.write(raw)
    return np.load(io.BytesIO(raw), allow_pickle=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="amazon_computers", choices=list(URLS))
    ap.add_argument("--edges", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--cache", default="")
    a = ap.parse_args()

    z = load_npz(URLS[a.dataset], a.cache)
    n = int(z["adj_shape"][0])
    adj = sp.csr_matrix(
        (z["adj_data"], z["adj_indices"], z["adj_indptr"]),
        shape=tuple(z["adj_shape"]),
    )
    adj = adj.maximum(adj.T)                 # symmetrise (undirected)
    upper = sp.triu(adj, k=1).tocoo()        # i<j, drop self-loops
    labels = np.asarray(z["labels"]).astype(int)

    with open(a.edges, "w") as fe:
        for u, v in zip(upper.row.tolist(), upper.col.tolist()):
            fe.write(f"{u + 1} {v + 1}\n")   # 1-based, matches the LFR/R side
    with open(a.labels, "w") as fl:
        for lab in labels.tolist():
            fl.write(f"{lab}\n")

    sys.stderr.write(
        f"{a.dataset} n={n} m={upper.nnz} classes={len(np.unique(labels))}\n")


if __name__ == "__main__":
    main()
