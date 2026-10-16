#!/usr/bin/env python3
"""Derive peak thresholds from calibration background ONLY.

For a grid of candidate thresholds, run the frozen peak/cluster definition over
calibration background tracks and measure merged candidates per Mb of background
scanned. The operating thresholds are the ones achieving the pre-declared target
densities (WORKFLOW_SETUP.md section 8: 1, 10 and 100 merged background hits/Mb).

No positive, no seed19 window and no validation window is read here. Thresholds are
written once and are meant to be frozen before any evaluation.
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

TARGET_DENSITIES = [1.0, 10.0, 100.0]   # merged candidates per Mb of background
GRID = np.round(np.arange(0.30, 1.3501, 0.01), 4)


def density_at(tracks, win, tau, contig_col):
    """Merged background candidates per Mb at this threshold."""
    cands, scanned = [], 0
    seen = set()
    for (rid, orient), (v, _, _) in tracks.items():
        if rid not in seen:
            scanned += len(v)
            seen.add(rid)
        peaks, _, _, _ = pc.call_peaks(v, 'background_calibrated', tau)
        w = win[rid]
        for c in pc.form_clusters(peaks):
            cands.append({'contig': w[contig_col], 'window_id': rid,
                          'start0': int(w['start0']) + c['start'],
                          'end0': int(w['start0']) + c['end'],
                          'score': c['score'], 'n_peaks': c['n_peaks'],
                          'median_peak_mean': c['median_peak_mean'],
                          'orientation': orient, 'peak_lengths': '', 'gaps': ''})
    merged = pc.merge_genomic(cands) if cands else []
    return len(merged), scanned, len(cands)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True,
                   help='calibration BACKGROUND score dirs only')
    p.add_argument('--windows', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--contig-col', default='accession_version')
    p.add_argument('--mode-label', required=True, help='exact | rp_g20_r2')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    win = {r['record_id']: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    tracks = load_runs(a.scores)
    fams = {win[rid].get('family', '') for rid, _ in tracks}
    if fams - {'calibration_background', 'background', ''}:
        raise SystemExit(f'Non-background material in calibration input: {fams}')

    rows = []
    for tau in GRID:
        n_merged, scanned, n_raw = density_at(tracks, win, float(tau), a.contig_col)
        rows.append({'tau': float(tau), 'merged_candidates': n_merged,
                     'per_window_clusters': n_raw, 'background_nt': scanned,
                     'per_Mb': n_merged / (scanned / 1e6)})
        if n_merged == 0:
            break
    with (a.out / f'sweep_{a.mode_label}.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    chosen = {}
    for target in TARGET_DENSITIES:
        ok = [r for r in rows if r['per_Mb'] <= target]
        pick = min(ok, key=lambda r: r['tau']) if ok else None
        chosen[f'bg{target:g}_per_Mb'] = (
            {'tau': pick['tau'], 'achieved_per_Mb': round(pick['per_Mb'], 3),
             'merged_candidates': pick['merged_candidates'],
             'unreachable': pc.threshold_is_unreachable(pick['tau'])}
            if pick else {'tau': None, 'note': 'target not reachable on this grid'})
    out = {'mode': a.mode_label, 'target_densities_per_Mb': TARGET_DENSITIES,
           'grid': [float(GRID[0]), float(GRID[-1]), 0.01],
           'background_nt_scanned': rows[0]['background_nt'],
           'background_windows': len({rid for rid, _ in tracks}),
           'peak_cluster_params': pc.PARAMS,
           'calibrated_on': [str(s) for s in a.scores],
           'no_positives_used': True, 'thresholds': chosen}
    (a.out / f'thresholds_{a.mode_label}.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['thresholds'], indent=2))
    print(f"background: {out['background_windows']} windows, "
          f"{rows[0]['background_nt'] / 1e6:.3f} Mb (single orientation count)")


if __name__ == '__main__':
    main()
