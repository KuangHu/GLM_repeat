#!/usr/bin/env python3
"""Sliding-window pilot scan: candidates and yield per unique scanned genomic Mb.

Used identically for calibration and held-out material. Candidates are merged
across strands AND across overlapping windows, so the 50% stride overlap and the
two orientations never produce independent candidates. Density is reported per
unique scanned genomic Mb (union of window spans per contig, counted once), not
per summed window length and not per orientation.

Modes:
  calibrate : sweep tau on CALIBRATION material and report yield per Mb
  scan      : apply a frozen tau and emit candidates
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

TARGETS = [100.0, 10.0]          # primary discovery point first, stricter tier second
GRID = np.round(np.arange(0.30, 1.3501, 0.01), 4)


def unique_scanned_nt(windows, ids):
    """Union of genomic spans of the emitted windows, per contig, counted once."""
    by = {}
    for rid in ids:
        w = windows[rid]
        by.setdefault(w['contig'], []).append((int(w['start0']), int(w['end0'])))
    total, per_contig = 0, {}
    for contig, spans in by.items():
        spans.sort()
        cur_s, cur_e, acc = spans[0][0], spans[0][1], 0
        for s, e in spans[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                acc += cur_e - cur_s
                cur_s, cur_e = s, e
        acc += cur_e - cur_s
        per_contig[contig] = acc
        total += acc
    return total, per_contig


def candidates_at(tracks, windows, tau):
    """Per-window clusters from both strands, merged into genomic candidates."""
    raw = []
    for (rid, orient), (v, _, _) in tracks.items():
        peaks, _, _, _ = pc.call_peaks(v, 'background_calibrated', tau)
        w = windows[rid]
        w0 = int(w['start0'])
        for c in pc.form_clusters(peaks):
            raw.append({'contig': w['contig'], 'window_id': rid, 'orientation': orient,
                        'start0': w0 + c['start'], 'end0': w0 + c['end'],
                        'score': c['score'], 'n_peaks': c['n_peaks'],
                        'median_peak_mean': c['median_peak_mean'],
                        'peak_lengths': ','.join(map(str, c['peak_lengths'])),
                        'gaps': ','.join(map(str, c['gaps']))})
    merged = pc.merge_genomic(raw) if raw else []
    for m in merged:
        m['span_nt'] = m['end0'] - m['start0']
        m['source_windows'] = ','.join(sorted(set(m['source_windows'])))
    return merged, raw


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['calibrate', 'scan'])
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--windows', type=Path, required=True, help='make_windows windows.tsv')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--threshold', type=float, help='frozen tau, required for scan')
    p.add_argument('--threshold-label', default='')
    p.add_argument('--set-label', required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    windows = {r['window_id']: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    tracks = load_runs(a.scores)
    ids = sorted({rid for rid, _ in tracks})
    missing = [w for w in windows if w not in set(ids)]
    scanned, per_contig = unique_scanned_nt(windows, ids)
    meta = next(iter(tracks.values()))[1]
    base = {'set': a.set_label, 'mode': meta['mode'], 'groups': meta.get('groups'),
            'repeats': meta.get('repeats'), 'seed': meta.get('seed'),
            'attention_mask_policy': meta.get('attention_mask_policy'),
            'model_revision': meta['revision'], 'peak_cluster_params': pc.PARAMS,
            'n_windows_scored': len(ids), 'n_windows_missing': len(missing),
            'unique_scanned_nt': scanned, 'unique_scanned_Mb': round(scanned / 1e6, 4),
            'per_contig_unique_nt': per_contig,
            'density_definition': 'merged candidates per unique scanned genomic Mb; '
                                  'strands and overlapping windows merged first'}
    if missing:
        base['missing_windows_note'] = ('windows with no track are NOT counted as '
                                        'negatives; they are excluded from scanned nt')
        scanned, per_contig = unique_scanned_nt(windows, ids)

    if a.mode == 'calibrate':
        rows = []
        for tau in GRID:
            merged, raw = candidates_at(tracks, windows, float(tau))
            rows.append({'tau': float(tau), 'merged_candidates': len(merged),
                         'per_window_clusters': len(raw),
                         'per_Mb': len(merged) / (scanned / 1e6)})
            if not merged:
                break
        with (a.out / 'sweep.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        # The density/tau curve is NOT monotonic: at low tau the whole window
        # exceeds the threshold, collapses into a single peak and fails MIN_PEAKS,
        # so yield is low for a degenerate reason. Operating points must be taken
        # from the stringent branch, at or above the density maximum.
        peak_row = max(rows, key=lambda r: r['per_Mb'])
        stringent = [r for r in rows if r['tau'] >= peak_row['tau']]
        chosen = {'_density_maximum': {'tau': peak_row['tau'],
                                       'per_Mb': round(peak_row['per_Mb'], 2),
                                       'note': 'operating points are taken at or '
                                               'above this tau (stringent branch)'}}
        for t in TARGETS:
            ok = [r for r in stringent if r['per_Mb'] <= t]
            pick = min(ok, key=lambda r: r['tau']) if ok else None
            chosen[f'bg{t:g}_per_Mb'] = (
                {'tau': pick['tau'], 'achieved_per_Mb': round(pick['per_Mb'], 2),
                 'merged_candidates': pick['merged_candidates'],
                 'unreachable': pc.threshold_is_unreachable(pick['tau'])}
                if pick else {'tau': None, 'note': 'not reachable on this grid'})
        out = {**base, 'target_densities_per_Mb': TARGETS, 'thresholds': chosen,
               'calibrated_on': [str(s) for s in a.scores],
               'no_heldout_or_positive_used': True}
        (a.out / 'thresholds.json').write_text(json.dumps(out, indent=2) + '\n')
        print(json.dumps(chosen, indent=2))
        print(f'calibration: {len(ids)} windows, {scanned / 1e6:.3f} unique Mb')
    else:
        if a.threshold is None:
            p.error('scan requires --threshold')
        merged, raw = candidates_at(tracks, windows, a.threshold)
        cols = ['contig', 'start0', 'end0', 'span_nt', 'score', 'n_peaks',
                'median_peak_mean', 'orientation', 'n_windows', 'best_window',
                'source_windows', 'peak_lengths', 'gaps']
        with (a.out / 'candidates.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(sorted(merged, key=lambda m: -m['score']))
        with (a.out / 'candidates.bed').open('w') as fh:
            for m in sorted(merged, key=lambda x: (x['contig'], x['start0'])):
                fh.write(f"{m['contig']}\t{m['start0']}\t{m['end0']}\t"
                         f"glmrep_n{m['n_peaks']}\t{m['score']:.4f}\t.\n")
        dens = {c: round(sum(1 for m in merged if m['contig'] == c) / (n / 1e6), 2)
                for c, n in per_contig.items()}
        out = {**base, 'threshold': a.threshold, 'threshold_label': a.threshold_label,
               'n_candidates': len(merged), 'n_per_window_clusters': len(raw),
               'candidates_per_Mb': round(len(merged) / (scanned / 1e6), 2),
               'per_contig_candidates_per_Mb': dens,
               'yield_not_fdr': 'these are candidate yields, not false-discovery rates'}
        (a.out / 'scan.json').write_text(json.dumps(out, indent=2) + '\n')
        print(f"{a.set_label}: {len(merged)} candidates over "
              f"{scanned / 1e6:.3f} unique Mb = {out['candidates_per_Mb']}/Mb")
        for c, d in sorted(dens.items(), key=lambda kv: -kv[1]):
            print(f'  {c}: {d}/Mb ({per_contig[c] / 1e6:.3f} Mb)')


if __name__ == '__main__':
    main()
