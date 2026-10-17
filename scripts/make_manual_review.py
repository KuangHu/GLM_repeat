#!/usr/bin/env python3
"""Build the stratified manual-inspection sheet for a scan round.

Strata, kept separate in the output so the random estimate is never mixed with the
enriched one:
  random      40 candidates drawn uniformly at random from ALL candidates; this is
              the only stratum usable for an unbiased yield estimate
  high        highest scoring
  adjacent    lowest surviving cluster score (nearest the formation boundary; this
              is NOT distance to tau, which is a per-base track threshold)
  structure   one per repeat architecture reported by the explainer
  neighbour   candidates overlapping annotated tRNA/rRNA/CDS neighbourhoods
  policy_lost ALL candidates the CDS policy would have dropped (never subsampled)

Random-sample IDs are retained in full: a candidate that also qualifies for an
enriched stratum keeps its `random_sample` flag and stays in the random
denominator. It is never removed from the random set.

Three independent judgement columns are provided and must not be collapsed into a
single true/false verdict:
  sequence_evidence     clear_repeat | suspected_repeat | no_repeat_found | cannot_judge
  explanation_coverage  model_adequate | partial | boundary_or_phase_problem |
                        structure_unsupported
  biological_context    known_common_structure | worth_followup | function_unknown
A candidate may qualify for several enriched strata; it is listed once with every
reason recorded.
"""
import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--n-random', type=int, default=40)
    p.add_argument('--n-per-stratum', type=int, default=10)
    p.add_argument('--seed', type=int, default=20261016)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cand = list(csv.DictReader((a.run / 'calls/audit_full.tsv').open(), delimiter='\t'))
    pol = list(csv.DictReader((a.run / 'calls/audit_policy.tsv').open(), delimiter='\t'))
    models = {}
    mp = a.run / 'explain_full/models.tsv'
    if mp.exists():
        for r in csv.DictReader(mp.open(), delimiter='\t'):
            models.setdefault(r['candidate_id'], r)
    cid = lambda c: f"{c['contig']}:{c['start0']}-{c['end0']}"
    ov = lambda a0, a1, b0, b1: max(0, min(a1, b1) - max(a0, b0))
    lost = [c for c in cand if not any(
        c['contig'] == q['contig'] and ov(int(c['start0']), int(c['end0']),
                                          int(q['start0']), int(q['end0'])) > 0 for q in pol)]
    rng = random.Random(a.seed)
    reasons = defaultdict(list)
    for c in rng.sample(cand, min(a.n_random, len(cand))):
        reasons[cid(c)].append('random')
    for c in sorted(cand, key=lambda r: -float(r['score']))[:a.n_per_stratum]:
        reasons[cid(c)].append('high')
    for c in sorted(cand, key=lambda r: float(r['score']))[:a.n_per_stratum]:
        reasons[cid(c)].append('adjacent')
    by_arch = defaultdict(list)
    for c in cand:
        m = models.get(cid(c))
        if m:
            by_arch[m.get('architecture', '')].append(c)
    for arch, cs in sorted(by_arch.items()):
        for c in sorted(cs, key=lambda r: -float(r['score']))[:2]:
            reasons[cid(c)].append(f'structure:{arch}')
    for c in lost:
        reasons[cid(c)].append('policy_lost')
    index = {cid(c): c for c in cand}
    rows = []
    for k, rs in reasons.items():
        c = index[k]
        m = models.get(k, {})
        rows.append({'candidate_id': k, 'strata': ','.join(sorted(set(rs))),
                     'random_sample': 'random' in rs,
                     'contig': c['contig'], 'start0': c['start0'], 'end0': c['end0'],
                     'score': c['score'], 'n_peaks': c.get('n_peaks', ''),
                     'explainer_status': m.get('status', 'no_model'),
                     'architecture': m.get('architecture', ''),
                     'consensus': m.get('consensus', ''),
                     'n_copies': m.get('n_copies', ''),
                     # three independent judgements; never collapse into true/false
                     'sequence_evidence': '',   # clear_repeat | suspected_repeat | no_repeat_found | cannot_judge
                     'explanation_coverage': '',# model_adequate | partial | boundary_or_phase_problem | structure_unsupported
                     'biological_context': '',  # known_common_structure | worth_followup | function_unknown
                     'notes': ''})
    rows.sort(key=lambda r: (not r['random_sample'], -float(r['score'])))
    with (a.out / 'manual_review.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    counts = defaultdict(int)
    for r in rows:
        for s in r['strata'].split(','):
            counts[s.split(':')[0]] += 1
    (a.out / 'summary.json').write_text(json.dumps({
        'total_candidates': len(cand), 'sheet_rows': len(rows),
        'random_stratum': sum(r['random_sample'] for r in rows),
        'policy_lost_all_included': len(lost),
        'stratum_counts': dict(counts), 'seed': a.seed,
        'note': 'only the random stratum supports an unbiased yield estimate; '
                'enriched strata are reported separately'}, indent=2) + '\n')
    print(f"{len(cand)} candidates -> {len(rows)} review rows "
          f"({sum(r['random_sample'] for r in rows)} random, {len(lost)} policy-lost)")
    for k, v in sorted(counts.items()):
        print(f'  {k:14s} {v}')


if __name__ == '__main__':
    main()
