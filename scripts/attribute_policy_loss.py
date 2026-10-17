#!/usr/bin/env python3
"""Attribute CDS-policy copy differences to windows, candidate boundaries or models.

A copy difference between the full-grid and policy runs can arise from three
distinct causes, which must not be reported as one number:

  window_excluded    a window supporting the full-grid candidate is not in the
                     policy window set, so that sequence was never scored
  outside_explainer_input
                     the policy candidate's envelope moved, so the copy's
                     coordinates fall outside (policy candidate +/- flank) and
                     the explainer never saw them
  model_change       the sequence was available to both runs, but the hypothesis
                     set or repeat-unit boundaries differ

Copies are deduplicated by genomic coordinate and reported as retained / lost /
gained, under BOTH exact-coordinate matching and interval matching at an explicit
reciprocal-overlap threshold. Loose overlap is never used to declare retention.
All hypotheses are compared, not only the primary display model.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def load_copies(path):
    """{candidate_id: {(contig,start0,end0): set(model_id)}} over ALL hypotheses."""
    by_cand = defaultdict(lambda: defaultdict(set))
    allc = defaultdict(set)
    for r in csv.DictReader(Path(path).open(), delimiter='\t'):
        key = (r['contig'], int(r['start0']), int(r['end0']))
        by_cand[r['candidate_id']][key].add(r['model_id'])
        allc[key].add(r['model_id'])
    return by_cand, allc


def interval_match(key, pool, min_recip):
    c, s, e = key
    best = None
    for (c2, s2, e2) in pool:
        if c2 != c:
            continue
        o = ov(s, e, s2, e2)
        if not o:
            continue
        recip = min(o / (e - s), o / (e2 - s2))
        if recip >= min_recip and (best is None or recip > best[1]):
            best = ((c2, s2, e2), recip)
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--flank', type=int, default=500)
    p.add_argument('--min-reciprocal-overlap', type=float, default=0.9)
    p.add_argument('--loose-reciprocal-overlap', type=float, default=0.5,
                   help='permissive threshold for detecting phase/boundary shifts; '
                        'short copies lose reciprocal overlap fast under small shifts')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    R = a.run
    fcand = list(csv.DictReader((R / 'calls/audit_full.tsv').open(), delimiter='\t'))
    pcand = list(csv.DictReader((R / 'calls/audit_policy.tsv').open(), delimiter='\t'))
    fby, fall = load_copies(R / 'explain_full/copies.tsv')
    pby, pall = load_copies(R / 'explain_policy/copies.tsv')
    wins = list(csv.DictReader((R / 'prepared/windows_all.tsv').open(), delimiter='\t'))
    policy_win = {w['window_id'] for w in wins if w['policy_selected'].lower() in ('true', '1', 'yes')}

    cid = lambda c: f"{c['contig']}:{c['start0']}-{c['end0']}"

    # ---- global copy-set relation, deduplicated by genomic coordinate
    retained = sorted(set(fall) & set(pall))
    lost = sorted(set(fall) - set(pall))
    gained = sorted(set(pall) - set(fall))
    lost_iv = {k: interval_match(k, pall, a.min_reciprocal_overlap) for k in lost}
    lost_strict = [k for k, v in lost_iv.items() if v is None]
    summary = {
        'min_reciprocal_overlap': a.min_reciprocal_overlap,
        'full_copies': len(fall), 'policy_copies': len(pall),
        'retained_exact': len(retained), 'lost_exact': len(lost), 'gained_exact': len(gained),
        'net_change': len(pall) - len(fall),
        'lost_exact_but_interval_matched': len(lost) - len(lost_strict),
        'lost_under_interval_matching': len(lost_strict),
        'note': ('percentages describe explainer copy RECORDS, not verified '
                 'biological repeats; differences may reflect candidate boundary '
                 'or hypothesis changes, not only filtered sequence'),
    }

    with (a.out / 'copy_relation.tsv').open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['relation', 'contig', 'start0', 'end0', 'full_models',
                    'policy_models', 'interval_match', 'reciprocal_overlap'])
        for k in retained:
            w.writerow(['retained_exact', *k, len(fall[k]), len(pall[k]), '', ''])
        for k in lost:
            m = lost_iv[k]
            w.writerow(['lost', *k, len(fall[k]), 0,
                        f'{m[0][0]}:{m[0][1]}-{m[0][2]}' if m else '',
                        f'{m[1]:.3f}' if m else ''])
        for k in gained:
            w.writerow(['gained', *k, 0, len(pall[k]), '', ''])

    # ---- per-candidate attribution for affected candidates
    rows = []
    for f in fcand:
        f_id = cid(f)
        fs, fe = int(f['start0']), int(f['end0'])
        hits = [q for q in pcand if q['contig'] == f['contig']
                and ov(fs, fe, int(q['start0']), int(q['end0'])) > 0]
        fcopies = set(fby.get(f_id, {}))
        pcopies = set()
        for q in hits:
            pcopies |= set(pby.get(cid(q), {}))
        missing = fcopies - pcopies
        if hits and not missing:
            continue
        # windows supporting the full candidate region (candidate +/- flank)
        lo, hi = fs - a.flank, fe + a.flank
        sup = [w for w in wins if w['contig'] == f['contig']
               and ov(lo, hi, int(w['start0']), int(w['end0'])) > 0]
        sup_excluded = [w['window_id'] for w in sup if w['window_id'] not in policy_win]
        # policy explainer input range for the matched candidates
        pranges = [(int(q['start0']) - a.flank, int(q['end0']) + a.flank) for q in hits]
        for k in sorted(missing):
            c, s, e = k
            in_excluded_window = any(ov(s, e, int(w['start0']), int(w['end0'])) > 0
                                     for w in sup if w['window_id'] in sup_excluded)
            in_policy_input = any(ov(s, e, s2, e2) == (e - s) for s2, e2 in pranges)
            if in_excluded_window and not in_policy_input:
                cause = 'window_excluded'
            elif not in_policy_input:
                cause = 'outside_explainer_input'
            elif in_excluded_window:
                cause = 'window_excluded_but_within_policy_input'
            else:
                cause = 'model_change'
            # Representation outcome is a SEPARATE axis from the upstream trigger:
            # a copy may be re-represented with shifted boundaries, or have no
            # counterpart at all. Reported at a strict and a permissive threshold.
            iv = interval_match(k, pall, a.min_reciprocal_overlap)
            iv_loose = interval_match(k, pall, a.loose_reciprocal_overlap)
            any_ov = interval_match(k, pall, 1e-9)
            if iv:
                repr_out = 'boundary_shift_within_strict_threshold'
            elif iv_loose:
                repr_out = 'phase_or_boundary_shift'
            elif any_ov:
                repr_out = 'partial_overlap_only'
            else:
                repr_out = 'no_corresponding_copy'
            rows.append({'full_candidate': f_id, 'full_score': f['score'],
                         'representation_outcome': repr_out,
                         'envelope_status': 'lost' if not hits else 'overlapped_but_shrank',
                         'contig': c, 'copy_start0': s, 'copy_end0': e,
                         'cause': cause,
                         'supporting_windows': len(sup),
                         'supporting_windows_excluded_by_policy': len(sup_excluded),
                         'interval_match_elsewhere':
                             f'{iv_loose[0][0]}:{iv_loose[0][1]}-{iv_loose[0][2]}' if iv_loose else '',
                         'reciprocal_overlap': f'{iv_loose[1]:.3f}' if iv_loose else ''})
    if rows:
        with (a.out / 'attribution.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
    causes = defaultdict(int)
    repr_counts = defaultdict(int)
    cross = defaultdict(int)
    for r in rows:
        causes[r['cause']] += 1
        repr_counts[r['representation_outcome']] += 1
        cross[f"{r['cause']} | {r['representation_outcome']}"] += 1
    summary['affected_candidates'] = len({r['full_candidate'] for r in rows})
    summary['attributed_copies'] = len(rows)
    summary['upstream_trigger'] = dict(causes)
    summary['representation_outcome'] = dict(repr_counts)
    summary['trigger_x_outcome'] = dict(cross)
    summary['axis_note'] = (
        'upstream_trigger says which pipeline stage changed; representation_outcome '
        'says whether the copy still appears in some form. They are independent: a '
        'copy whose window was excluded may still be re-represented from another '
        'window with shifted boundaries. Only no_corresponding_copy means the '
        'evidence is absent from the policy output, and even that is a computational '
        'evidence difference, not a verified biological miss.')
    (a.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
