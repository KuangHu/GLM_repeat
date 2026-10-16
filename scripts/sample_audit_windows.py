#!/usr/bin/env python3
"""Sample ~24 held-out windows for the exact audit, in their ORIGINAL contexts.

Three strata, deterministic:
  high          windows carrying the highest-scoring primary-tier candidates
  adjacent      windows with the LOWEST surviving cluster score, i.e. nearest the
                cluster-FORMATION boundary. NOTE: this is not distance to tau.
                tau is a per-base track threshold; a cluster score is
                median(peak_mean)*log2(1+n_peaks) and lives on a different scale
                (floor ~tau*log2(3) for a minimal 2-peak cluster). The stratum is
                one-sided: windows below the boundary have no cluster at all and
                form the separate `non_hit` stratum.
  non_hit       windows with no candidate at the primary tier, drawn at random

The emitted FASTA contains each window's original 3,000-nt sequence exactly as the
sliding-window pipeline produced it; no candidate interval is excised and no
context is trimmed.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import peak_cluster as pc  # noqa: E402
from compare_tracks import load_runs  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--windows', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--threshold', type=float, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--per-stratum', type=int, default=8)
    p.add_argument('--seed', type=int, default=20261015)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    windows = {r['window_id']: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    tracks = load_runs(a.scores)

    best = {}
    for (rid, orient), (v, _, _) in tracks.items():
        peaks, _, _, _ = pc.call_peaks(v, 'background_calibrated', a.threshold)
        cl = pc.form_clusters(peaks)
        s = max((c['score'] for c in cl), default=0.0)
        best[rid] = max(best.get(rid, 0.0), s)
    hits = {k: v for k, v in best.items() if v > 0}
    misses = sorted(k for k, v in best.items() if v == 0)

    rng = np.random.default_rng(a.seed)
    high = [k for k, _ in sorted(hits.items(), key=lambda kv: -kv[1])[:a.per_stratum]]
    rest = sorted(set(hits) - set(high))
    # lowest surviving cluster scores = nearest the cluster-formation boundary.
    # NOT abs(cluster_score - tau): those quantities are not on the same scale.
    adjacent = [k for k, _ in sorted(((k, hits[k]) for k in rest),
                                     key=lambda kv: kv[1])[:a.per_stratum]]
    non_hit = sorted(rng.choice(misses, size=min(a.per_stratum, len(misses)),
                                replace=False).tolist())
    strata = {'high': high, 'adjacent': adjacent, 'non_hit': non_hit}

    seqs, name = {}, None
    for line in a.fasta.open():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line)
    seqs = {k: ''.join(v) for k, v in seqs.items()}

    rows = []
    with (a.out / 'audit_windows.fa').open('w') as fa:
        for stratum, ids in strata.items():
            for rid in ids:
                w = windows[rid]
                fa.write(f'>{rid}\n{seqs[rid]}\n')
                rows.append({'window_id': rid, 'stratum': stratum,
                             'rp_best_cluster_score': round(best[rid], 4),
                             'contig': w['contig'], 'start0': w['start0'],
                             'end0': w['end0'], 'length_nt': len(seqs[rid]),
                             'gc_fraction': w['gc_fraction'],
                             'window_sha256': w['window_sha256']})
    with (a.out / 'audit_windows.tsv').open('w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        wr.writeheader()
        wr.writerows(rows)
    (a.out / 'records.txt').write_text('\n'.join(r['window_id'] for r in rows) + '\n')
    (a.out / 'sampling.json').write_text(json.dumps({
        'strata_sizes': {k: len(v) for k, v in strata.items()},
        'primary_threshold_used_for_sampling': a.threshold,
        'n_windows_with_candidate': len(hits), 'n_windows_without': len(misses),
        'seed': a.seed,
        'adjacent_stratum_selection_rule':
            'the 8 windows with the lowest surviving cluster score (nearest the '
            'cluster-formation boundary), NOT the 8 nearest in value to tau; '
            'cluster score and tau are different scales',
        'note': 'sequences are the original 3,000-nt sliding windows, untrimmed; '
                'the RP threshold is used only to STRATIFY, and is not transferred '
                'to exact scoring'}, indent=2) + '\n')
    for k, v in strata.items():
        print(f'{k}: {len(v)} windows')
    print(f'{len(rows)} windows -> {a.out}')


if __name__ == '__main__':
    main()
