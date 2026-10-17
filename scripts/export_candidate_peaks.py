#!/usr/bin/env python3
"""Export the exact peaks the frozen detector called, per candidate.

The candidate tables record peak LENGTHS and gaps but not peak coordinates, so a
reviewer cannot see where inside a candidate the signal actually sits. This
replays `peak_cluster.call_peaks` on each window track at the frozen threshold --
the same call `pilot_scan.candidates_at` makes -- and writes the resulting peaks
in genomic coordinates.

These are the detector's own peaks, not a recomputation on a merged or averaged
track: peaks are called per window and per orientation exactly as during the
scan, then attributed to whichever candidates they overlap. A peak may therefore
appear more than once for a candidate, from different windows or orientations,
which is what the detector saw.

Read-only with respect to every frozen artefact.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import peak_cluster as pc


def tsv(path):
    p = Path(path)
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--threshold-file', type=Path, required=True)
    ap.add_argument('--threshold-key', default='bg100_per_Mb')
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    tau = float(json.loads(a.threshold_file.read_text())['thresholds'][a.threshold_key]['tau'])
    if pc.threshold_is_unreachable(tau):
        raise ValueError(f'tau={tau} >= ln(4); peaks would be impossible by construction')

    wins = {w['window_id']: w for w in tsv(a.run / 'prepared/windows.tsv')}
    cands = tsv(a.run / 'calls/candidates.tsv')
    by_contig = defaultdict(list)
    for c in cands:
        by_contig[c['contig']].append(c)

    idx = defaultdict(list)
    for rj in sorted((a.run / 'scores').glob('**/run.json')):
        meta = json.loads(rj.read_text())
        for rec in meta['records']:
            idx[rec['id']].append((rj.parent / rec['track'], rec['scoring_orientation']))

    rows = []
    for wid, w in wins.items():
        ws = int(w['start0'])
        local = [c for c in by_contig.get(w['contig'], [])
                 if int(c['start0']) < int(w['end0']) and ws < int(c['end0'])]
        if not local:
            continue
        for path, orient in idx.get(wid, []):
            v = np.load(path)['above_uniform'].astype(float)
            peaks, used_tau, med, sd = pc.call_peaks(v, 'background_calibrated', tau)
            for p in peaks:
                gs, ge = ws + p['start'], ws + p['end']
                for c in local:
                    if gs < int(c['end0']) and int(c['start0']) < ge:
                        rows.append({
                            'candidate_id': c['candidate_id'], 'contig': w['contig'],
                            'peak_start0': gs, 'peak_end0': ge, 'peak_len': ge - gs,
                            'summit0': ws + p['summit'],
                            'peak_mean_above_uniform': round(p['mean'], 5),
                            'window_id': wid, 'orientation': orient,
                            'window_median': round(med, 5), 'window_robust_sd': round(sd, 5),
                            'tau': used_tau,
                            'contained_in_candidate':
                                gs >= int(c['start0']) and ge <= int(c['end0'])})
    rows.sort(key=lambda r: (r['candidate_id'], r['peak_start0'], r['orientation']))
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    per = defaultdict(set)
    for r in rows:
        per[r['candidate_id']].add((r['peak_start0'], r['peak_end0']))
    print(f'{len(rows)} peak records over {len(per)} of {len(cands)} candidates '
          f'(tau={tau}, ln4={pc.LN4:.4f}) -> {a.out}')


if __name__ == '__main__':
    main()
