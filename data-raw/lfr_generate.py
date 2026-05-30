#!/usr/bin/env python3
"""Canonical LFR benchmark generator (Lancichinetti, Fortunato & Radicchi 2008).

Replaces the degree-corrected SBM proxy used previously. Generates one LFR
graph for a given mixing parameter mu and writes:
  - <edges>: one "u v" per line, 1-based node ids (for the R/igraph side);
  - <truth>: one community label per line, node order 1..n.

networkx's LFR_benchmark_graph can fail to converge for some (mu, seed); we
retry with shifted seeds a few times before giving up (nonzero exit).

Used by data-raw/30_lfr_robustness.R via the .venv-lfr interpreter.
"""
import argparse
import sys
import networkx as nx
from networkx.generators.community import LFR_benchmark_graph


def generate(n, mu, seed, tau1, tau2, avg_degree, max_degree, min_comm, max_comm,
             retries=12):
    last = None
    for k in range(retries):
        try:
            G = LFR_benchmark_graph(
                n=n, tau1=tau1, tau2=tau2, mu=mu,
                average_degree=avg_degree, max_degree=max_degree,
                min_community=min_comm, max_community=max_comm,
                tol=1e-7, max_iters=2000, seed=seed + 100 * k,
            )
            G.remove_edges_from(nx.selfloop_edges(G))
            return G
        except Exception as e:  # nx.ExceededMaxIterations and friends
            last = e
    raise RuntimeError(f"LFR failed after {retries} retries (mu={mu}): {last}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tau1", type=float, default=2.5)   # degree power law
    ap.add_argument("--tau2", type=float, default=1.5)   # community-size power law
    ap.add_argument("--avg-degree", type=float, default=12.0)
    ap.add_argument("--max-degree", type=int, default=50)
    ap.add_argument("--min-comm", type=int, default=20)
    ap.add_argument("--max-comm", type=int, default=60)
    ap.add_argument("--edges", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args()

    G = generate(a.n, a.mu, a.seed, a.tau1, a.tau2, a.avg_degree, a.max_degree,
                 a.min_comm, a.max_comm)

    # community membership: each node carries a frozenset "community" attribute
    comm_of = {}
    for node in G.nodes():
        comm = G.nodes[node]["community"]
        comm_of[node] = min(comm)  # canonical label = smallest member id
    labels = sorted(set(comm_of.values()))
    relabel = {c: i + 1 for i, c in enumerate(labels)}  # 1-based contiguous

    nodes = sorted(G.nodes())
    idx = {nd: i + 1 for i, nd in enumerate(nodes)}      # 1-based node ids

    with open(a.edges, "w") as fe:
        for u, v in G.edges():
            if u != v:
                fe.write(f"{idx[u]} {idx[v]}\n")
    with open(a.truth, "w") as ft:
        for nd in nodes:
            ft.write(f"{relabel[comm_of[nd]]}\n")

    sys.stderr.write(
        f"LFR n={G.number_of_nodes()} m={G.number_of_edges()} "
        f"comms={len(labels)} mu={a.mu}\n")


if __name__ == "__main__":
    main()
