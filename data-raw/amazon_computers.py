#!/usr/bin/env python3
"""Fetch the Amazon-Computers co-purchase network (Shchur et al. 2018) and emit
an undirected edge list + a ground-truth label per node, in the same format the
R side already consumes for LFR (1-based "u v" per line; one label per line in
node order). The 32 MB .npz is cached so re-runs do not re-download.

Source: gnn-benchmark/data/npz/amazon_electronics_computers.npz
  (Shchur, Mumme, Bojchevski & Gunnemann, "Pitfalls of Graph Neural Network
   Evaluation", 2018). 13,752 nodes, 10 product-category classes, undirected.

Usage:
  amazon_computers.py --edges E.txt --labels L.txt [--cache PATH.npz]
"""
import argparse
import io
import os
import sys
import urllib.request

import numpy as np
import scipy.sparse as sp

URL = ("https://github.com/shchur/gnn-benchmark/raw/master/"
       "data/npz/amazon_electronics_computers.npz")


def load_npz(cache):
    if cache and os.path.exists(cache):
        with open(cache, "rb") as f:
            raw = f.read()
    else:
        raw = urllib.request.urlopen(URL, timeout=120).read()
        if cache:
            with open(cache, "wb") as f:
                f.write(raw)
    return np.load(io.BytesIO(raw), allow_pickle=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--cache", default="")
    a = ap.parse_args()

    z = load_npz(a.cache)
    n = int(z["adj_shape"][0])
    adj = sp.csr_matrix(
        (z["adj_data"], z["adj_indices"], z["adj_indptr"]),
        shape=tuple(z["adj_shape"]),
    )
    # symmetrise (datasets are undirected) and drop self-loops, then take i<j
    adj = adj.maximum(adj.T)
    upper = sp.triu(adj, k=1).tocoo()
    labels = np.asarray(z["labels"]).astype(int)

    with open(a.edges, "w") as fe:
        for u, v in zip(upper.row.tolist(), upper.col.tolist()):
            fe.write(f"{u + 1} {v + 1}\n")          # 1-based, matches LFR side
    with open(a.labels, "w") as fl:
        for lab in labels.tolist():                  # one per node, node order
            fl.write(f"{lab}\n")

    sys.stderr.write(
        f"Amazon-Computers n={n} m={upper.nnz} "
        f"classes={len(np.unique(labels))}\n")


if __name__ == "__main__":
    main()
