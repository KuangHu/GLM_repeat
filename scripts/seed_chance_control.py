#!/usr/bin/env python3
"""Chance calibration for the explainer's 12-nt seed in each candidate region.

Many candidates receive a `repeat_pair` model built on a 12-nt seed that occurs
twice, often with both copies OUTSIDE the candidate envelope. Before such a model
can be read as evidence of a repeat, we need to know how often a region of this
length and composition yields a repeated 12-mer by chance.

For each candidate the explainer's own input region (candidate +/- flank) is
shuffled many times with a dinucleotide-preserving shuffle (Altschul-Erikson
Eulerian-path method), because DNA dinucleotide bias -- not just base
composition -- drives short-word recurrence. The statistic is the number of
distinct 12-mers occurring at least twice, using the same >=3-distinct-bases
filter the probe applies.

Diagnostic only. It never changes a candidate, a model or any frozen output; it
tells a reviewer whether a given seed is surprising.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path


def read_fasta(path):
    out, name = {}, None
    for line in Path(path).open():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v).upper() for k, v in out.items()}


def tsv(path):
    p = Path(path)
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def dinuc_shuffle(seq, rng):
    """Altschul-Erikson dinucleotide-preserving shuffle."""
    if len(seq) < 4:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    last = seq[-1]
    # pick a random last-edge per vertex forming a tree rooted at `last`
    for _ in range(12):
        tree, ok = {}, True
        for v in edges:
            if v != last:
                tree[v] = rng.choice(edges[v])
        for v in list(tree):                       # check connectivity to `last`
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in tree:
                    ok = False
                    break
                seen.add(cur)
                cur = tree[cur]
            if not ok:
                break
        if ok:
            break
    if not ok:
        return ''.join(rng.sample(seq, len(seq)))  # fall back, flagged by caller
    rest = {}
    for v, outs in edges.items():
        o = list(outs)
        if v in tree:
            o.remove(tree[v])
        rng.shuffle(o)
        rest[v] = o + ([tree[v]] if v in tree else [])
    out, cur = [seq[0]], seq[0]
    idx = defaultdict(int)
    for _ in range(len(seq) - 1):
        nxt = rest[cur][idx[cur]]
        idx[cur] += 1
        out.append(nxt)
        cur = nxt
    return ''.join(out)


def repeated_kmers(seq, k=12):
    c = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    return sum(1 for m, n in c.items() if n > 1 and len(set(m)) > 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='append', required=True, metavar='LABEL=DIR')
    ap.add_argument('--selection', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--flank', type=int, default=500)
    ap.add_argument('--k', type=int, default=12)
    ap.add_argument('--shuffles', type=int, default=300)
    ap.add_argument('--seed', type=int, default=20261017)
    a = ap.parse_args()
    runs = dict(s.partition('=')[::2] for s in a.run)
    G = {lab: read_fasta(Path(d) / 'input.fa') for lab, d in runs.items()}
    cand = {}
    for lab, d in runs.items():
        for c in tsv(Path(d) / 'calls/candidates.tsv'):
            cand[(lab, c['candidate_id'])] = c
    rng = random.Random(a.seed)
    rows = []
    for r in tsv(a.selection):
        lab, cid = r['run'], r['candidate_id']
        c = cand[(lab, cid)]
        g = G[lab][c['contig']]
        s0, e0 = int(c['start0']), int(c['end0'])
        lo, hi = max(0, s0 - a.flank), min(len(g), e0 + a.flank)
        region = g[lo:hi]
        obs = repeated_kmers(region, a.k)
        null = [repeated_kmers(dinuc_shuffle(region, rng), a.k) for _ in range(a.shuffles)]
        ge = sum(1 for x in null if x >= obs)
        rows.append({
            'run': lab, 'candidate_id': cid, 'region_nt': hi - lo,
            'observed_repeated_12mers': obs,
            'null_mean': round(sum(null) / len(null), 2),
            'null_max': max(null),
            'p_empirical': round((ge + 1) / (len(null) + 1), 4),
            'surprising_at_0.05': (ge + 1) / (len(null) + 1) < 0.05,
        })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    n_sig = sum(r['surprising_at_0.05'] for r in rows)
    print(f'{len(rows)} regions; {n_sig} have more repeated {a.k}-mers than '
          f'{1 - 0.05:.0%} of dinucleotide-preserving shuffles')
    print(f'-> {a.out}')


if __name__ == '__main__':
    main()
