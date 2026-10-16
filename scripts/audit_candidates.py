#!/usr/bin/env python3
"""Candidate geometry and merging audit. CPU only: no new model inference.

Recomputes peaks, per-window clusters and merged candidates from existing score
tracks, and reports for each final candidate: span, contributing pre-merge cluster
count, distinct source windows and orientations, peak-support coverage and internal
peak gaps. Also distinguishes the two ways a large cluster count can collapse into
few candidates:

  MIN_PEAKS rejection : clusters never formed because a window produced < 2 peaks
  transitive merging  : clusters formed, then chained into one interval because
                        each overlapped the next (50% window stride makes this
                        possible across a whole contig)
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import peak_cluster as pc  # noqa: E402
from compare_tracks import load_runs  # noqa: E402
from pilot_scan import unique_scanned_nt  # noqa: E402


def collect(tracks, windows, tau):
    """Peaks and clusters in genomic coordinates, plus MIN_PEAKS accounting."""
    raw, peaks_g = [], []
    n_single_peak_windows = n_zero_peak = n_peaks_total = 0
    for (rid, orient), (v, _, _) in tracks.items():
        pks, _, _, _ = pc.call_peaks(v, 'background_calibrated', tau)
        w0 = int(windows[rid]['start0'])
        contig = windows[rid]['contig']
        n_peaks_total += len(pks)
        if not pks:
            n_zero_peak += 1
        elif len(pks) == 1:
            n_single_peak_windows += 1
        for pk in pks:
            peaks_g.append((contig, w0 + pk['start'], w0 + pk['end']))
        for c in pc.form_clusters(pks):
            raw.append({'contig': contig, 'window_id': rid, 'orientation': orient,
                        'start0': w0 + c['start'], 'end0': w0 + c['end'],
                        'score': c['score'], 'n_peaks': c['n_peaks'],
                        'median_peak_mean': c['median_peak_mean'],
                        'peak_spans': [(w0 + p['start'], w0 + p['end'])
                                       for p in c['peaks']],
                        'peak_lengths': '', 'gaps': ''})
    stats = {'tracks': len(tracks), 'peaks_total': n_peaks_total,
             'tracks_with_0_peaks': n_zero_peak,
             'tracks_with_1_peak_rejected_by_MIN_PEAKS': n_single_peak_windows,
             'per_window_clusters': len(raw)}
    return raw, peaks_g, stats


def union_spans(spans):
    if not spans:
        return []
    spans = sorted(spans)
    out = [list(spans[0])]
    for s, e in spans[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--windows', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--taus', type=float, nargs='+', required=True)
    p.add_argument('--labels', nargs='+', required=True)
    p.add_argument('--set-label', required=True)
    a = p.parse_args()
    if len(a.taus) != len(a.labels):
        p.error('--taus and --labels must match')
    a.out.mkdir(parents=True, exist_ok=True)
    windows = {r['window_id']: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    tracks = load_runs(a.scores)
    ids = sorted({rid for rid, _ in tracks})
    scanned, _ = unique_scanned_nt(windows, ids)

    summary, tiers = {}, {}
    for tau, label in zip(a.taus, a.labels):
        raw, peaks_g, stats = collect(tracks, windows, tau)
        merged = pc.merge_genomic(raw) if raw else []
        by_contig_peaks = {}
        for contig, s, e in peaks_g:
            by_contig_peaks.setdefault(contig, []).append((s, e))
        by_contig_peaks = {k: union_spans(v) for k, v in by_contig_peaks.items()}
        rows = []
        for m in merged:
            s, e = m['start0'], m['end0']
            contrib = [r for r in raw if r['contig'] == m['contig']
                       and r['start0'] < e and s < r['end0']]
            pk = [(max(s, ps), min(e, pe)) for ps, pe in by_contig_peaks.get(m['contig'], [])
                  if ps < e and s < pe]
            pk = union_spans(pk)
            cov = sum(b - aa for aa, b in pk)
            gaps = [pk[i + 1][0] - pk[i][1] for i in range(len(pk) - 1)]
            rows.append({
                'contig': m['contig'], 'start0': s, 'end0': e, 'span_nt': e - s,
                'score': round(m['score'], 4), 'n_peaks_best_cluster': m['n_peaks'],
                'pre_merge_clusters': len(contrib),
                'source_windows': len({r['window_id'] for r in contrib}),
                'orientations': ''.join(sorted({r['orientation'] for r in contrib})),
                'peak_blocks_in_span': len(pk),
                'peak_support_nt': cov,
                'peak_support_fraction': round(cov / (e - s), 4) if e > s else 0.0,
                'max_internal_peak_gap_nt': max(gaps) if gaps else 0,
                'median_internal_peak_gap_nt': int(np.median(gaps)) if gaps else 0})
        rows.sort(key=lambda r: -r['span_nt'])
        with (a.out / f'geometry_{label}.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        spans = [r['span_nt'] for r in rows]
        union_cov = sum(b - aa for c in {r['contig'] for r in rows}
                        for aa, b in union_spans([(r['start0'], r['end0'])
                                                  for r in rows if r['contig'] == c]))
        summary[label] = {
            'tau': tau, **stats, 'merged_candidates': len(rows),
            'candidates_per_Mb': round(len(rows) / (scanned / 1e6), 2),
            'union_candidate_span_nt': union_cov,
            'union_candidate_span_fraction_of_scanned': round(union_cov / scanned, 4),
            'span_median': int(np.median(spans)) if spans else 0,
            'span_p95': int(np.percentile(spans, 95)) if spans else 0,
            'span_max': max(spans) if spans else 0,
            'pre_merge_clusters_per_candidate_median':
                float(np.median([r['pre_merge_clusters'] for r in rows])) if rows else 0,
            'pre_merge_clusters_per_candidate_max':
                max([r['pre_merge_clusters'] for r in rows], default=0),
            'source_windows_per_candidate_median':
                float(np.median([r['source_windows'] for r in rows])) if rows else 0,
            'source_windows_per_candidate_max':
                max([r['source_windows'] for r in rows], default=0),
            'peak_support_fraction_median':
                float(np.median([r['peak_support_fraction'] for r in rows])) if rows else 0,
            'max_internal_peak_gap_overall':
                max([r['max_internal_peak_gap_nt'] for r in rows], default=0)}
        tiers[label] = rows
        print(f"{label} tau={tau}: {len(rows)} candidates, median span "
              f"{summary[label]['span_median']} nt, p95 {summary[label]['span_p95']}, "
              f"max {summary[label]['span_max']}")

    # cross-tier correspondence (no nesting assumed)
    if len(a.labels) == 2:
        A, B = a.labels
        def olap(x, ys):
            return [y for y in ys if y['contig'] == x['contig']
                    and y['start0'] < x['end0'] and x['start0'] < y['end0']]
        a_in_b = sum(1 for x in tiers[A] if olap(x, tiers[B]))
        b_in_a = sum(1 for x in tiers[B] if olap(x, tiers[A]))
        summary['cross_tier'] = {
            'primary': A, 'stringent': B,
            'both_tiers_from_primary_side': a_in_b,
            'primary_only': len(tiers[A]) - a_in_b,
            'both_tiers_from_stringent_side': b_in_a,
            'stringent_only': len(tiers[B]) - b_in_a,
            'nesting_not_assumed': 'counts are computed from each side separately; '
                                   'a stringent candidate need not lie inside a '
                                   'primary one'}
    out = {'set': a.set_label, 'unique_scanned_nt': scanned,
           'unique_scanned_Mb': round(scanned / 1e6, 4),
           'peak_cluster_params': pc.PARAMS, 'tiers': summary,
           'no_new_inference': True}
    (a.out / 'audit.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(summary.get('cross_tier', {}), indent=2))


if __name__ == '__main__':
    main()
