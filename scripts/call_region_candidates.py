#!/usr/bin/env python3
"""Call frozen detector on selected windows, validating completion and provenance.

All score tracks must be in input-forward coordinates, as timing_score writes.
Candidates in incidental CDS context are retained and tagged, not silently lost.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from compare_tracks import load_runs
from pilot_scan import candidates_at
from peak_cluster import PARAMS
from region_io import SpanIndex, new_dir, provenance, table, write_table

META_KEYS = ['revision','mode','groups','repeats','seed','attention_mask_policy']
FIELDS = ['candidate_id','contig','start0','end0','span_nt','score','n_peaks',
          'median_peak_mean','orientation','n_windows','best_window','source_windows',
          'peak_lengths','gaps','target_overlap_nt','selection_class','audit_contig',
          'merge_geometry','n_member_clusters','biological_locus_status']


def validate_tracks(tracks, windows):
    expected = {(wid, ori) for wid in windows for ori in ('+','-')}
    if set(tracks) != expected:
        raise ValueError(f'Incomplete/unexpected tracks: missing={len(expected-set(tracks))}, '
                         f'extra={len(set(tracks)-expected)}. Never count missing tracks as negatives.')
    reference = None
    for (wid, _), (v, meta, rec) in tracks.items():
        w = windows[wid]
        if len(v) != int(w['end0'])-int(w['start0']) or not np.isfinite(v).all():
            raise ValueError(f'Invalid score length or non-finite values: {wid}')
        if rec['sequence_sha256'] != w['window_sha256']:
            raise ValueError(f'Sequence hash mismatch: {wid}')
        if rec.get('coordinates') != 'input_fasta_record_0_based':
            raise ValueError(f'Unknown score coordinate convention: {wid}')
        signature = {k: meta[k] for k in META_KEYS}
        if reference is not None and signature != reference:
            raise ValueError('Mixed model/scoring settings across score tracks')
        reference = signature
    return reference


def decorate(merged, raw, targets, audit):
    members = []
    by = defaultdict(list)
    for r in raw:
        by[r['contig']].append(r)
    for m in merged:
        c,s,e = m['contig'],m['start0'],m['end0']
        mid = f'{c}:{s}-{e}'
        sub = [r for r in by[c] if r['start0'] < e and r['end0'] > s]
        common = min(r['end0'] for r in sub) > max(r['start0'] for r in sub)
        overlap = targets[c].overlap(s,e)
        m.update(candidate_id=mid, target_overlap_nt=overlap,
                 selection_class='target_associated' if overlap else 'audit_only' if c in audit else 'context_only',
                 audit_contig=c in audit, merge_geometry='common_overlap' if common else 'transitive_chain',
                 n_member_clusters=len(sub), biological_locus_status='unresolved_envelope')
        m['n_windows'] = len({r['window_id'] for r in sub})
        members.extend(dict(candidate_id=mid, **r) for r in sub)
    return members


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--scores', type=Path, nargs='*', default=[])
    p.add_argument('--out', type=Path, required=True)
    threshold = p.add_mutually_exclusive_group(required=True)
    threshold.add_argument('--threshold-file', type=Path, help='Frozen pilot thresholds.json')
    threshold.add_argument('--threshold', type=float, help='Exploratory tau, no density/FDR guarantee')
    p.add_argument('--threshold-key', default='bg100_per_Mb')
    a = p.parse_args()
    rows = table(a.prepared/'windows.tsv')
    windows = {r['window_id']: r for r in rows}
    if len(windows) != len(rows):
        raise ValueError('Duplicate window IDs')
    tracks = load_runs(a.scores)
    signature = validate_tracks(tracks, windows)
    calibration = json.loads(a.threshold_file.read_text()) if a.threshold_file else None
    tau = calibration['thresholds'][a.threshold_key]['tau'] if calibration else a.threshold
    if tau is None or not math.isfinite(tau) or not 0 <= tau < math.log(4):
        raise ValueError('Threshold must be finite, >=0 and < ln(4)')
    if calibration and signature:
        expected = {k: calibration['model_revision' if k=='revision' else k] for k in META_KEYS}
        if expected != signature or calibration['peak_cluster_params'] != PARAMS:
            raise ValueError('Frozen threshold model/scorer settings do not match these tracks')
        selection = json.loads((a.prepared/'provenance.json').read_text(encoding='utf-8'))['parameters']
        if (selection['window'],selection['stride'],selection['min_segment']) != (3000,1500,40):
            raise ValueError('This frozen pilot threshold requires the original 3000/1500/40 window setup; '
                             'changed context needs independent recalibration')
    target_rows = table(a.prepared/'targets.tsv')
    by = defaultdict(list)
    for r in target_rows:
        by[r['contig']].append((int(r['start0']),int(r['end0'])))
    contigs = table(a.prepared/'contigs.tsv')
    targets = {r['contig']: SpanIndex(by[r['contig']]) for r in contigs}
    audit = {r['contig'] for r in contigs if r['audit_contig']=='True'}
    merged, raw = candidates_at(tracks, windows, tau)
    members = decorate(merged, raw, targets, audit)
    # Counterfactual on the EXACT same audit tracks, with all non-policy windows removed.
    audit_tracks = {k:v for k,v in tracks.items() if windows[k[0]]['contig'] in audit}
    policy_tracks = {k:v for k,v in audit_tracks.items() if windows[k[0]]['policy_selected']=='True'}
    audit_full, _ = candidates_at(audit_tracks, windows, tau)
    audit_policy, _ = candidates_at(policy_tracks, windows, tau)
    for r in audit_full:
        r['overlaps_policy_candidate'] = any(
            x['contig']==r['contig'] and x['start0']<r['end0'] and x['end0']>r['start0'] for x in audit_policy)
    out = new_dir(a.out)
    write_table(out/'candidates.tsv', sorted(merged,key=lambda r:-r['score']), FIELDS)
    write_table(out/'members.tsv', members, ['candidate_id','contig','start0','end0','window_id',
                                            'orientation','score','n_peaks','peak_lengths','gaps'])
    write_table(out/'audit_full.tsv', audit_full, FIELDS[:13]+['overlaps_policy_candidate'])
    write_table(out/'audit_policy.tsv', audit_policy, FIELDS[:13])
    scanned = {}
    for r in rows:
        scanned.setdefault(r['contig'], []).append((int(r['start0']),int(r['end0'])))
    unique_nt = sum(SpanIndex(v).prefix[-1] for v in scanned.values())
    target_scanned = sum(targets[c].overlap(s,e) for c,sp in scanned.items() for s,e in SpanIndex(sp).spans)
    n_target = sum(r['selection_class']=='target_associated' for r in merged)
    result = dict(n_candidates=len(merged), n_windows=len(windows), unique_scanned_nt=unique_nt,
                  target_scanned_nt=target_scanned, n_target_associated=n_target,
                  all_candidates_per_scanned_Mb=len(merged)*1e6/unique_nt if unique_nt else None,
                  target_candidates_per_target_Mb=n_target*1e6/target_scanned if target_scanned else None,
                  audit_full_candidates=len(audit_full), audit_policy_candidates=len(audit_policy),
                  audit_full_without_policy_overlap=sum(not r['overlaps_policy_candidate'] for r in audit_full),
                  threshold=tau, scoring=signature,
                  limitations=['yields are NOT FDR or biological probabilities',
                    'old full-genome calibration does not guarantee yield after CDS selection',
                    'overlap audit is candidate-envelope recovery, not repeat copy recall',
                    'transitively merged envelopes may contain multiple biological loci'])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    provenance(out,a,[a.prepared/'windows.tsv',a.prepared/'targets.tsv',a.threshold_file],
               score_manifests={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                                for d in a.scores for p in d.glob('**/run.json')})
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
