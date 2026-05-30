#!/usr/bin/env python3
"""ECG (Ensemble Clustering for Graphs; Poulin & Theberge 2019) membership.

ECG is not in core igraph; it lives in the author's `partition-igraph`
package. This helper reads a 1-based edge list, runs community_ecg, and writes
one community label per node (node order 1..n, 1-based labels) so the R side
can compute NMI/ARI/Q. Used as a recent, runnable partition-quality baseline
that reportedly beats Louvain/Leiden on LFR recovery.

Run via the .venv-lfr interpreter (pip install partition-igraph).
"""
import argparse
import igraph
import partition_igraph  # noqa: F401  (registers Graph.community_ecg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--edges", required=True)   # 1-based "u v" per line
    ap.add_argument("--out", required=True)     # 1-based label per node, order 1..n
    ap.add_argument("--ens", type=int, default=16)  # ensemble size
    a = ap.parse_args()

    edges = []
    with open(a.edges) as fe:
        for line in fe:
            line = line.split()
            if len(line) >= 2:
                u, v = int(line[0]) - 1, int(line[1]) - 1   # -> 0-based
                if u != v:
                    edges.append((u, v))
    g = igraph.Graph(n=a.n, edges=edges, directed=False)
    g.simplify()
    ec = g.community_ecg(ens_size=a.ens)
    mem = ec.membership                          # 0-based labels, node order 0..n-1
    with open(a.out, "w") as fo:
        for lab in mem:
            fo.write(f"{int(lab) + 1}\n")        # -> 1-based labels


if __name__ == "__main__":
    main()
