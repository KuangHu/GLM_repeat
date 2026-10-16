#!/usr/bin/env python3
"""Label-blind peak-cluster scorer (design: results/scorer/PEAK_CLUSTER_DESIGN.md).

Reads per-base score tracks and emits peaks, clusters and window scores. Labels,
annotations, motifs and family names are never inputs. Clusters are emitted in
genomic coordinates and merged across overlapping windows, so overlapping windows
never contribute independent peaks.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_tracks import load_runs  # noqa: E402

# --- frozen parameters (see PEAK_CLUSTER_DESIGN.md; not re-tuned this round) ---
SMOOTH_NT = 16
THRESH_K = 3.0
MIN_PEAK_NT = 8
MAX_GAP_NT = 500
MIN_PEAKS = 2
PARAMS = dict(SMOOTH_NT=SMOOTH_NT, THRESH_K=THRESH_K, MIN_PEAK_NT=MIN_PEAK_NT,
              MAX_GAP_NT=MAX_GAP_NT, MIN_PEAKS=MIN_PEAKS,
              cluster_score='median(peak_mean) * log2(1 + n_peaks)',
              window_score='max over clusters and orientations')


def smooth(v, w=SMOOTH_NT):
    """Centred running mean, same length as the input."""
    if len(v) < w:
        return np.full(len(v), float(v.mean()))
    pad = np.pad(v.astype(np.float64), (w // 2, w - 1 - w // 2), mode='edge')
    return np.convolve(pad, np.ones(w) / w, mode='valid')


THRESHOLD_POLICIES = ('window_local', 'background_calibrated')
LN4 = float(np.log(4))          # ceiling of the above-uniform measurement


def call_peaks(track, policy='window_local', tau_fixed=None):
    """Maximal runs above the threshold set by `policy`.

    window_local         : median + THRESH_K * robust SD, computed per window.
                           Preserved unchanged as the original scorer.
    background_calibrated: a single `tau_fixed` derived from calibration background
                           only (see calibrate_threshold.py). Never fitted on
                           positives.
    """
    s = smooth(track)
    med = float(np.median(s))
    sd = 1.4826 * float(np.median(np.abs(s - med)))
    if policy == 'background_calibrated':
        if tau_fixed is None:
            raise ValueError('background_calibrated requires tau_fixed')
        tau = float(tau_fixed)
    elif policy == 'window_local':
        tau = med + THRESH_K * sd
    else:
        raise ValueError(f'Unknown threshold policy: {policy}')
    above = s > tau
    peaks, start = [], None
    for i, b in enumerate(above):
        if b and start is None:
            start = i
        elif not b and start is not None:
            if i - start >= MIN_PEAK_NT:
                peaks.append((start, i))
            start = None
    if start is not None and len(above) - start >= MIN_PEAK_NT:
        peaks.append((start, len(above)))
    return [{'start': a, 'end': b, 'mean': float(s[a:b].mean()),
             'summit': int(a + np.argmax(s[a:b]))} for a, b in peaks], tau, med, sd


def threshold_is_unreachable(tau):
    """A tau at or above ln(4) can never be exceeded, so the window yields zero
    peaks by construction. This is a property of the threshold, NEVER evidence
    that the window lacks repetitive organization."""
    return tau >= LN4


def form_clusters(peaks):
    """Group peaks separated by <= MAX_GAP_NT; keep those with >= MIN_PEAKS."""
    if not peaks:
        return []
    groups, cur = [], [peaks[0]]
    for p in peaks[1:]:
        if p['start'] - cur[-1]['end'] <= MAX_GAP_NT:
            cur.append(p)
        else:
            groups.append(cur)
            cur = [p]
    groups.append(cur)
    out = []
    for g in groups:
        if len(g) < MIN_PEAKS:
            continue
        med = float(np.median([p['mean'] for p in g]))
        out.append({'start': g[0]['start'], 'end': g[-1]['end'], 'n_peaks': len(g),
                    'median_peak_mean': med,
                    'score': med * float(np.log2(1 + len(g))),
                    'peak_lengths': [p['end'] - p['start'] for p in g],
                    'gaps': [g[i + 1]['start'] - g[i]['end'] for i in range(len(g) - 1)],
                    'peaks': g})
    return out


def merge_genomic(cands):
    """Merge overlapping candidates on the same contig across windows."""
    out = []
    for contig in sorted({c['contig'] for c in cands}):
        rows = sorted([c for c in cands if c['contig'] == contig],
                      key=lambda c: (c['start0'], c['end0']))
        cur = None
        for c in rows:
            if cur and c['start0'] < cur['end0']:
                cur['end0'] = max(cur['end0'], c['end0'])
                cur['n_windows'] += 1
                cur['source_windows'].append(c['window_id'])
                if c['score'] > cur['score']:
                    for k in ('score', 'n_peaks', 'median_peak_mean', 'orientation',
                              'peak_lengths', 'gaps', 'best_window'):
                        cur[k] = c[k] if k != 'best_window' else c['window_id']
            else:
                if cur:
                    out.append(cur)
                cur = {**c, 'n_windows': 1, 'source_windows': [c['window_id']],
                       'best_window': c['window_id']}
        if cur:
            out.append(cur)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--windows', type=Path, required=True,
                   help='TSV with record_id/window_id, accession/contig and start0')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--id-col', default='record_id')
    p.add_argument('--threshold-policy', choices=THRESHOLD_POLICIES,
                   default='window_local')
    p.add_argument('--threshold', type=float,
                   help='tau for background_calibrated (from calibrate_threshold.py)')
    p.add_argument('--threshold-label', default='',
                   help='provenance tag, e.g. bg10_per_Mb')
    p.add_argument('--contig-col', default='accession_version')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    win = {r[a.id_col]: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    tracks = load_runs(a.scores)

    if a.threshold_policy == 'background_calibrated' and a.threshold is None:
        p.error('background_calibrated requires --threshold')
    peak_rows, cand_rows, window_rows, unreachable = [], [], {}, []
    for (rid, orient), (v, meta, rec) in sorted(tracks.items()):
        w = win[rid]
        w0, contig = int(w['start0']), w[a.contig_col]
        peaks, tau, med, sd = call_peaks(v, a.threshold_policy, a.threshold)
        if threshold_is_unreachable(tau):
            unreachable.append({'window_id': rid, 'orientation': orient,
                                'threshold': round(tau, 4), 'ln4': round(LN4, 4)})
        for i, pk in enumerate(peaks, 1):
            peak_rows.append({'window_id': rid, 'orientation': orient, 'peak_index': i,
                              'contig': contig, 'start0': w0 + pk['start'],
                              'end0': w0 + pk['end'], 'length_nt': pk['end'] - pk['start'],
                              'peak_mean': round(pk['mean'], 4),
                              'summit0': w0 + pk['summit'],
                              'threshold': round(tau, 4), 'window_median': round(med, 4),
                              'window_robust_sd': round(sd, 4)})
        for c in form_clusters(peaks):
            cand_rows.append({'window_id': rid, 'orientation': orient, 'contig': contig,
                              'start0': w0 + c['start'], 'end0': w0 + c['end'],
                              'span_nt': c['end'] - c['start'], 'n_peaks': c['n_peaks'],
                              'median_peak_mean': round(c['median_peak_mean'], 4),
                              'score': round(c['score'], 4),
                              'peak_lengths': ','.join(map(str, c['peak_lengths'])),
                              'gaps': ','.join(map(str, c['gaps']))})
        best = max([c['score'] for c in form_clusters(peaks)], default=0.0)
        window_rows[rid] = max(window_rows.get(rid, 0.0), best)

    merged = merge_genomic(cand_rows)
    for m in merged:
        m['source_windows'] = ','.join(sorted(set(m['source_windows'])))
    for name, rows in [('peaks.tsv', peak_rows), ('clusters_per_window.tsv', cand_rows)]:
        with (a.out / name).open('w', newline='') as fh:
            if rows:
                wr = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
                wr.writeheader()
                wr.writerows(rows)
    if merged:
        cols = ['contig', 'start0', 'end0', 'span_nt', 'score', 'n_peaks',
                'median_peak_mean', 'orientation', 'n_windows', 'best_window',
                'source_windows', 'peak_lengths', 'gaps']
        with (a.out / 'candidates_merged.tsv').open('w', newline='') as fh:
            wr = csv.DictWriter(fh, cols, delimiter='\t', extrasaction='ignore')
            wr.writeheader()
            wr.writerows(sorted(merged, key=lambda m: -m['score']))
        with (a.out / 'candidates_merged.bed').open('w') as fh:
            for m in sorted(merged, key=lambda x: (x['contig'], x['start0'])):
                fh.write(f"{m['contig']}\t{m['start0']}\t{m['end0']}\t"
                         f"glmrep_cluster_n{m['n_peaks']}\t{m['score']:.4f}\t.\n")
    with (a.out / 'window_scores.tsv').open('w', newline='') as fh:
        wr = csv.writer(fh, delimiter='\t')
        wr.writerow(['record_id', 'peak_cluster_score'])
        for rid in sorted(window_rows):
            wr.writerow([rid, f'{window_rows[rid]:.4f}'])
    any_meta = next(iter(tracks.values()))[1]
    (a.out / 'run.json').write_text(json.dumps({
        'params': PARAMS, 'threshold_policy': a.threshold_policy,
        'threshold': a.threshold, 'threshold_label': a.threshold_label,
        'ln4_ceiling': LN4,
        'unreachable_threshold_tracks': unreachable,
        'unreachable_note': 'a threshold >= ln(4) cannot be exceeded; such a window '
                            'yields zero peaks by construction and that is NEVER '
                            'evidence that no repetitive organization exists',
        'design': 'results/scorer/PEAK_CLUSTER_DESIGN.md',
        'authorship': 'agent-designed, pre-declared before any track was scored',
        'source_scores': [str(s) for s in a.scores],
        'track_mode': any_meta['mode'], 'model_revision': any_meta['revision'],
        'attention_mask_policy': any_meta.get('attention_mask_policy', 'always'),
        'n_tracks': len(tracks), 'n_windows': len(window_rows),
        'n_peaks': len(peak_rows), 'n_clusters_per_window': len(cand_rows),
        'n_candidates_merged': len(merged),
        'label_free': 'no label, annotation, motif or family information is an input',
    }, indent=2) + '\n')
    print(f'{len(tracks)} tracks -> {len(peak_rows)} peaks, {len(cand_rows)} clusters, '
          f'{len(merged)} merged candidates -> {a.out}')


if __name__ == '__main__':
    main()
