#!/usr/bin/env python3
"""Build a co-authorship graph from OpenAlex with an EXTERNAL per-author impact
signal (citations, #works) for validating leader designation.

Fetches works matching an OpenAlex `--filter` (cursor-paginated, capped at
`--max-works`); for each work it records co-authorship edges and adds the work's
`cited_by_count` to each of its authors. Authors with at least `--min-works`
works in the slice are kept; the rest is denoise. Emits, in the format the R
side consumes:
  --edges : 1-based "u v" undirected co-authorship edges (one per line)
  --impact: per node "citations n_works" in node order (1..K)
The raw works JSON is cached so re-runs do not re-query the API.

Usage:
  openalex_coauthor.py --filter '...' --max-works 4000 --min-works 2 \
      --edges E.txt --impact I.txt [--cache C.json] [--mailto you@x.org]
"""
import argparse
import itertools
import json
import os
import sys
import urllib.parse
import urllib.request

API = "https://api.openalex.org/works"


def fetch_works(filt, max_works, mailto, cache):
    if cache and os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    out, cursor, got = [], "*", 0
    while got < max_works:
        q = {"filter": filt, "per-page": 200, "cursor": cursor,
             "select": "id,authorships,cited_by_count", "mailto": mailto}
        url = API + "?" + urllib.parse.urlencode(q)
        page = json.load(urllib.request.urlopen(url, timeout=120))
        res = page.get("results", [])
        if not res:
            break
        out.extend(res)
        got += len(res)
        cursor = page.get("meta", {}).get("next_cursor")
        if not cursor:
            break
    out = out[:max_works]
    if cache:
        with open(cache, "w") as f:
            json.dump(out, f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filter", required=True)
    ap.add_argument("--max-works", type=int, default=4000)
    ap.add_argument("--min-works", type=int, default=2)
    ap.add_argument("--edges", required=True)
    ap.add_argument("--impact", required=True)
    ap.add_argument("--cache", default="")
    ap.add_argument("--mailto", default="leite@de.ufpe.br")
    a = ap.parse_args()

    works = fetch_works(a.filter, a.max_works, a.mailto, a.cache)

    citations, nworks, pairs = {}, {}, {}
    for w in works:
        auth = [au["author"]["id"] for au in w.get("authorships", [])
                if au.get("author", {}).get("id")]
        auth = list(dict.fromkeys(auth))            # unique, keep order
        c = int(w.get("cited_by_count", 0) or 0)
        for aid in auth:
            citations[aid] = citations.get(aid, 0) + c
            nworks[aid] = nworks.get(aid, 0) + 1
        for u, v in itertools.combinations(auth, 2):
            e = (u, v) if u < v else (v, u)
            pairs[e] = pairs.get(e, 0) + 1

    keep = {aid for aid, k in nworks.items() if k >= a.min_works}
    ids = sorted(keep)
    idx = {aid: i + 1 for i, aid in enumerate(ids)}   # 1-based node ids

    n_edges = 0
    with open(a.edges, "w") as fe:
        for (u, v) in pairs:
            if u in keep and v in keep:
                fe.write(f"{idx[u]} {idx[v]}\n"); n_edges += 1
    with open(a.impact, "w") as fi:
        for aid in ids:                                # node order 1..K
            fi.write(f"{citations[aid]} {nworks[aid]}\n")

    sys.stderr.write(
        f"OpenAlex works={len(works)} authors_kept={len(ids)} edges={n_edges} "
        f"(min_works={a.min_works})\n")


if __name__ == "__main__":
    main()
