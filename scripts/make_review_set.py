#!/usr/bin/env python3
"""First candidate review set, after deduplication.

Deduplicates explained candidates into repeat-pattern clusters by consensus similarity
(strand-canonical, substring or >=80% ungapped identity), then selects:
  * interpretable `unclassified` representatives spanning different architectures
  * several high-scoring candidates with NO sequence explanation

These are sequence-level descriptions only. An `unclassified` candidate is one
whose span carries no overlapping GenBank annotation in the records used here; that
is a statement about available annotation, NOT a claim of novelty, function, or a
new biological system.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans('ACGT', 'TGCA')


def canonical(s):
    rc = s.translate(COMP)[::-1]
    return min(s, rc)


def similar(a, b, ident=0.8):
    if a in b or b in a:
        return True
    if abs(len(a) - len(b)) > 2:
        return False
    n = min(len(a), len(b))
    best = 0
    for off in (0, 1, -1):
        x = a[max(0, off):max(0, off) + n]
        y = b[max(0, -off):max(0, -off) + n]
        m = min(len(x), len(y))
        if m:
            best = max(best, sum(p == q for p, q in zip(x[:m], y[:m])) / m)
    return best >= ident


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--unresolved', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--per-architecture', type=int, default=3)
    p.add_argument('--n-unexplained', type=int, default=8)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    models = list(csv.DictReader(a.models.open(), delimiter='\t'))
    unres = list(csv.DictReader(a.unresolved.open(), delimiter='\t'))

    # ---- deduplicate into repeat-pattern clusters by consensus similarity
    fam, reps = {}, []
    for m in sorted(models, key=lambda r: -float(r['detector_score'] or 0)):
        c = canonical(m['consensus'])
        hit = None
        for i, r in enumerate(reps):
            if similar(c, r):
                hit = i
                break
        if hit is None:
            reps.append(c)
            hit = len(reps) - 1
        fam[m['candidate_id']] = hit
    members = defaultdict(list)
    for m in models:
        members[fam[m['candidate_id']]].append(m)

    rows = []
    for fid, ms in members.items():
        best = max(ms, key=lambda r: float(r['detector_score'] or 0))
        rows.append({'pattern_cluster_id': fid, 'n_members': len(ms),
                     'representative_candidate': best['candidate_id'],
                     'architecture': best['architecture'],
                     'consensus': best['consensus'], 'unit_len': best['unit_len'],
                     'n_copies': best['n_copies'],
                     'median_start_to_start': best['median_start_to_start'],
                     'spacer_cv': best['spacer_cv'],
                     'information_bits_per_position': best['information_bits_per_position'],
                     'detector_score': best['detector_score'],
                     'annotation_category': best['annotation_category'],
                     'member_candidates': ','.join(m['candidate_id'] for m in ms[:8])})
    rows.sort(key=lambda r: -float(r['detector_score'] or 0))
    with (a.out / 'pattern_clusters.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # ---- review selection
    picked, by_arch = [], defaultdict(int)
    for r in rows:
        if r['annotation_category'] != 'unclassified':
            continue
        if by_arch[r['architecture']] >= a.per_architecture:
            continue
        by_arch[r['architecture']] += 1
        picked.append({**r, 'review_reason': 'interpretable unclassified representative'})
    unexpl = sorted(unres, key=lambda r: -float(r['score'] or 0))[:a.n_unexplained]
    for u in unexpl:
        picked.append({'pattern_cluster_id': '', 'n_members': 1,
                       'representative_candidate': u['candidate_id'],
                       'architecture': 'none_found', 'consensus': '',
                       'unit_len': '', 'n_copies': 0, 'median_start_to_start': '',
                       'spacer_cv': '', 'information_bits_per_position': '',
                       'detector_score': u['score'], 'annotation_category': '',
                       'member_candidates': '',
                       'review_reason': 'high-scoring, no sequence explanation; '
                                        'alternatives: ' + (u['alternatives'] or '[]')})
    with (a.out / 'review_set.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(picked[0]), delimiter='\t')
        w.writeheader()
        w.writerows(picked)
    (a.out / 'README.md').write_text(
        '# First candidate review set\n\n'
        f'{len(models)} explained candidates deduplicate into {len(rows)} '
        f'**repeat-pattern clusters** by consensus similarity (strand-canonical; '
        f'substring or >=80% ungapped identity). A repeat-pattern cluster is a group '
        f'of candidates whose consensus sequences resemble each other. It is NOT a '
        f'biological family, system or classification: sequence similarity alone '
        f'does not establish shared origin or function.\n\n'
        f'Selected for review: {len(picked) - len(unexpl)} interpretable '
        f'`unclassified` representatives spanning '
        f'{len(by_arch)} architectures, plus {len(unexpl)} high-scoring candidates '
        'with no sequence explanation.\n\n'
        '**These are sequence-level descriptions only.** `unclassified` means the '
        'candidate span carries no overlapping annotation in the GenBank records '
        'used here. It is not a claim of novelty, function, or a new biological '
        'system, and nothing here should be read as one.\n\n'
        'Alignment is ungapped: indel-containing arrays are under-reported.\n')
    print(f'{len(models)} explained -> {len(rows)} repeat-pattern clusters after dedup')
    print(f'review set: {len(picked)} entries '
          f'({len(picked) - len(unexpl)} unclassified representatives across '
          f'{len(by_arch)} architectures, {len(unexpl)} unexplained)')
    for k, v in by_arch.items():
        print(f'  {k:24s} {v}')


if __name__ == '__main__':
    main()
