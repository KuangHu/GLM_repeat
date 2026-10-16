#!/usr/bin/env python3
"""Reconcile called peaks with annotated repeat/core intervals, coordinate by coordinate.

Distinguishes two quantities that must not be conflated:
  edge gap        previous peak end -> next peak start (inter-peak spacer)
  start-to-start  previous peak start -> next peak start (periodicity/pitch)
An annotated unit spacing is a start-to-start quantity; an edge gap is smaller than
it by roughly the unit length, so the two are never directly comparable.

Each peak is classified against the annotated units it overlaps:
  matched    overlaps exactly one annotated unit
  merged     overlaps two or more annotated units (one peak spanning several)
  extra      inside the locus span but overlapping no annotated unit
  off_label  outside the locus span entirely
and each annotated unit is marked covered or missing.
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


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--windows', type=Path, required=True)
    p.add_argument('--units', type=Path, required=True,
                   help='TSV: contig,start0,end0,unit_id,locus_id (annotated units)')
    p.add_argument('--loci', type=Path, required=True)
    p.add_argument('--threshold', type=float, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    windows = {r['window_id']: r for r in csv.DictReader(a.windows.open(), delimiter='\t')}
    units = list(csv.DictReader(a.units.open(), delimiter='\t'))
    loci = {r['locus_id']: r for r in csv.DictReader(a.loci.open(), delimiter='\t')}
    tracks = load_runs(a.scores)

    # peaks in genomic coordinates, unioned per contig across windows/strands
    raw = {}
    for (rid, orient), (v, _, _) in tracks.items():
        w = windows[rid]
        w0 = int(w['start0'])
        pks, _, _, _ = pc.call_peaks(v, 'background_calibrated', a.threshold)
        for pk in pks:
            raw.setdefault(w['contig'], []).append((w0 + pk['start'], w0 + pk['end']))
    merged_peaks = {}
    for contig, spans in raw.items():
        spans.sort()
        out = [list(spans[0])]
        for s, e in spans[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged_peaks[contig] = [tuple(x) for x in out]

    peak_rows, locus_rows, unit_rows = [], [], []
    for locus_id, lo in loci.items():
        contig = lo['contig']
        ls, le = int(lo['start0']), int(lo['end0'])
        us = sorted([(int(u['start0']), int(u['end0']), u['unit_id'])
                     for u in units if u['locus_id'] == locus_id])
        pk = [(s, e) for s, e in merged_peaks.get(contig, []) if ov(ls, le, s, e) > 0]
        pk.sort()
        classed = []
        for s, e in pk:
            hits = [u for u in us if ov(s, e, u[0], u[1]) > 0]
            if not hits:
                kind = 'extra_inside_locus' if ov(ls, le, s, e) > 0 else 'off_label'
            elif len(hits) == 1:
                kind = 'matched'
            else:
                kind = f'merged_{len(hits)}_units'
            classed.append({'contig': contig, 'locus_id': locus_id, 'start0': s,
                            'end0': e, 'length_nt': e - s, 'kind': kind,
                            'units_overlapped': ','.join(h[2] for h in hits)})
        peak_rows += classed
        covered = {h for c in classed for h in c['units_overlapped'].split(',') if h}
        for s, e, uid in us:
            unit_rows.append({'contig': contig, 'locus_id': locus_id, 'unit_id': uid,
                              'start0': s, 'end0': e,
                              'status': 'covered' if uid in covered else 'MISSING'})
        starts = [c['start0'] for c in classed]
        edge_gaps = [classed[i + 1]['start0'] - classed[i]['end0']
                     for i in range(len(classed) - 1)]
        s2s = [starts[i + 1] - starts[i] for i in range(len(starts) - 1)]
        u_s2s = [us[i + 1][0] - us[i][0] for i in range(len(us) - 1)]
        u_edge = [us[i + 1][0] - us[i][1] for i in range(len(us) - 1)]
        med = lambda x: int(np.median(x)) if x else None
        locus_rows.append({
            'locus_id': locus_id, 'family': lo['family'], 'contig': contig,
            'locus_nt': le - ls, 'n_annotated_units': len(us),
            'unit_len_median': med([e - s for s, e, _ in us]),
            'annotated_start_to_start_median': med(u_s2s),
            'annotated_start_to_start_range': f'{min(u_s2s)}-{max(u_s2s)}' if u_s2s else '',
            'annotated_edge_gap_median': med(u_edge),
            'n_peaks_in_locus': len(classed),
            'peak_len_median': med([c['length_nt'] for c in classed]),
            'peak_start_to_start_median': med(s2s),
            'peak_start_to_start_range': f'{min(s2s)}-{max(s2s)}' if s2s else '',
            'peak_edge_gap_median': med(edge_gaps),
            'peak_edge_gap_max': max(edge_gaps) if edge_gaps else None,
            'n_matched': sum(c['kind'] == 'matched' for c in classed),
            'n_merged': sum(c['kind'].startswith('merged') for c in classed),
            'n_extra_inside_locus': sum(c['kind'] == 'extra_inside_locus' for c in classed),
            'n_units_covered': len(covered), 'n_units_missing': len(us) - len(covered)})

    for name, rows in [('peaks_classified.tsv', peak_rows),
                       ('units_status.tsv', unit_rows),
                       ('locus_reconciliation.tsv', locus_rows)]:
        if rows:
            with (a.out / name).open('w', newline='') as fh:
                w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
                w.writeheader()
                w.writerows(rows)
    (a.out / 'run.json').write_text(json.dumps({
        'threshold': a.threshold, 'peak_cluster_params': pc.PARAMS,
        'definitions': {
            'edge_gap': 'previous peak end -> next peak start',
            'start_to_start': 'previous peak start -> next peak start (pitch)',
            'warning': 'annotated unit spacing is start-to-start; an edge gap is '
                       'smaller by roughly the unit length. Never compare directly.'},
        'n_loci': len(locus_rows), 'n_peaks': len(peak_rows),
        'n_units': len(unit_rows)}, indent=2) + '\n')
    for r in locus_rows:
        print(f"{r['locus_id'][:38]:38s} units {r['n_annotated_units']:2d} "
              f"s2s {r['annotated_start_to_start_median']} | peaks "
              f"{r['n_peaks_in_locus']:2d} s2s {r['peak_start_to_start_median']} "
              f"edge {r['peak_edge_gap_median']} | matched {r['n_matched']} "
              f"merged {r['n_merged']} extra {r['n_extra_inside_locus']} "
              f"missing {r['n_units_missing']}")


if __name__ == '__main__':
    main()
