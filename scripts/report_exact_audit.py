#!/usr/bin/env python3
"""Exact audit of sampled held-out windows: compare exact and RP on the same windows.

Compares three things per window: the per-base tracks, the peak-support structure,
and candidate localization.

Two things this deliberately does NOT do:
  * It does not transfer the RP operating point to exact as a calibrated threshold.
    Peaks are called on both modes at a MATCHED tau purely so the peak structure is
    comparable; that matched value is NOT an exact operating point, because exact
    has no calibration at this density (see results/calibration/thresholds/FROZEN.md).
  * It does not treat exact output as biological ground truth. Exact is a reference
    measurement of the same model, not an independent observation of the genome.
"""
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import peak_cluster as pc  # noqa: E402
from compare_tracks import load_runs  # noqa: E402


def jacc(a, b):
    ia = np.zeros(max(len(a), 1), bool)
    return ia


def span_overlap(x, y):
    """Jaccard of two interval lists on a common coordinate length."""
    if not x and not y:
        return 1.0
    n = max([e for _, e in x + y] + [1])
    ax, ay = np.zeros(n, bool), np.zeros(n, bool)
    for s, e in x:
        ax[s:e] = True
    for s, e in y:
        ay[s:e] = True
    u = (ax | ay).sum()
    return float((ax & ay).sum() / u) if u else 1.0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exact', type=Path, nargs='+', required=True)
    p.add_argument('--rp', type=Path, nargs='+', required=True)
    p.add_argument('--sample', type=Path, required=True)
    p.add_argument('--matched-tau', type=float, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sample = {r['window_id']: r for r in csv.DictReader(a.sample.open(), delimiter='\t')}
    ex, rp = load_runs(a.exact), load_runs(a.rp)
    keys = sorted(set(ex) & set(rp))
    rows = []
    for rid, orient in keys:
        e, r = ex[(rid, orient)][0], rp[(rid, orient)][0]
        se, sr = pc.smooth(e), pc.smooth(r)
        pe, _, _, _ = pc.call_peaks(e, 'background_calibrated', a.matched_tau)
        pr, _, _, _ = pc.call_peaks(r, 'background_calibrated', a.matched_tau)
        ce, cr = pc.form_clusters(pe), pc.form_clusters(pr)
        rows.append({
            'window_id': rid, 'stratum': sample[rid]['stratum'], 'orientation': orient,
            'rp_best_cluster_score': sample[rid]['rp_best_cluster_score'],
            'bias_rp_minus_exact': round(float((r - e).mean()), 4),
            'mae': round(float(np.abs(r - e).mean()), 4),
            'p99_abs': round(float(np.quantile(np.abs(r - e), .99)), 4),
            'pearson': round(float(np.corrcoef(e, r)[0, 1]), 4),
            'pearson_smoothed16': round(float(np.corrcoef(se, sr)[0, 1]), 4),
            'n_peaks_exact': len(pe), 'n_peaks_rp': len(pr),
            'peak_span_jaccard': round(span_overlap(
                [(x['start'], x['end']) for x in pe],
                [(x['start'], x['end']) for x in pr]), 4),
            'n_clusters_exact': len(ce), 'n_clusters_rp': len(cr),
            'cluster_span_jaccard': round(span_overlap(
                [(x['start'], x['end']) for x in ce],
                [(x['start'], x['end']) for x in cr]), 4),
            'best_cluster_exact': round(max([c['score'] for c in ce], default=0.0), 4),
            'best_cluster_rp': round(max([c['score'] for c in cr], default=0.0), 4)})
    with (a.out / 'per_track.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    by = defaultdict(list)
    for r in rows:
        by[r['stratum']].append(r)
    med = lambda xs: float(np.median(xs))
    L = ['# Exact audit of sampled held-out windows', '',
         f'24 windows x 2 orientations = {len(rows)} tracks, scored exactly in their '
         'original untrimmed 3,000-nt contexts.', '',
         '**The RP operating point was not transferred to exact.** Peaks below are '
         f'called on both modes at a matched tau = {a.matched_tau} solely so that peak '
         'structure is comparable. That is not a calibrated exact operating point: '
         'exact has no calibration at this density, and establishing one would need '
         '~1 Mb of exact background (~150 GPU-hours), which was explicitly declined.',
         '', '**Exact is not treated as ground truth.** It is a reference measurement '
         'of the same model under a different masking scheme, not an independent '
         'observation of the genome.', '',
         '## Agreement by stratum (median over tracks)', '',
         '| stratum | tracks | bias RP-exact | MAE | p99 abs | Pearson | Pearson 16-nt |'
         ' peak-span Jaccard | cluster-span Jaccard |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for s in ('high', 'adjacent', 'non_hit'):
        rs = by.get(s, [])
        if not rs:
            continue
        L.append(f"| {s} | {len(rs)} | {med([x['bias_rp_minus_exact'] for x in rs]):+.4f} "
                 f"| {med([x['mae'] for x in rs]):.4f} "
                 f"| {med([x['p99_abs'] for x in rs]):.4f} "
                 f"| {med([x['pearson'] for x in rs]):.4f} "
                 f"| {med([x['pearson_smoothed16'] for x in rs]):.4f} "
                 f"| {med([x['peak_span_jaccard'] for x in rs]):.3f} "
                 f"| {med([x['cluster_span_jaccard'] for x in rs]):.3f} |")
    L += ['', '## Cluster-level correspondence at matched tau', '',
          '| stratum | tracks | exact has cluster | RP has cluster | both | neither |',
          '|---|---:|---:|---:|---:|---:|']
    for s in ('high', 'adjacent', 'non_hit'):
        rs = by.get(s, [])
        if not rs:
            continue
        eh = sum(x['n_clusters_exact'] > 0 for x in rs)
        rh = sum(x['n_clusters_rp'] > 0 for x in rs)
        both = sum(x['n_clusters_exact'] > 0 and x['n_clusters_rp'] > 0 for x in rs)
        nei = sum(x['n_clusters_exact'] == 0 and x['n_clusters_rp'] == 0 for x in rs)
        L.append(f'| {s} | {len(rs)} | {eh} | {rh} | {both} | {nei} |')
    disagree = [r for r in rows if (r['n_clusters_exact'] > 0) != (r['n_clusters_rp'] > 0)]
    L += ['', f'Tracks where only one mode forms a cluster: {len(disagree)} of {len(rows)}.']
    for r in sorted(disagree, key=lambda x: -abs(x['best_cluster_exact'] - x['best_cluster_rp']))[:10]:
        L.append(f"- {r['window_id']} ({r['orientation']}, {r['stratum']}): "
                 f"exact {r['n_clusters_exact']} cluster(s) best {r['best_cluster_exact']}, "
                 f"RP {r['n_clusters_rp']} best {r['best_cluster_rp']}")
    L += ['', '## Reading', '',
          'The `adjacent` stratum is the one that matters: it was drawn at the '
          'cluster-formation boundary (lowest surviving RP cluster scores, '
          '1.305-1.337), so small track differences can flip a detection either way. '
          'Disagreement there is expected and is the quantity of interest, not a '
          'defect in either mode.', '',
          'Sampling-rule note: the `adjacent` stratum is NOT "closest in value to '
          'tau"; see results/pilot/heldout/exact_audit/sampling.json.']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
