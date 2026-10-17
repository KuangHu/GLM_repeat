#!/usr/bin/env python3
"""Attach pilercr CRISPR-array annotation to candidates.

A regularly spaced array with an invariant short core is exactly what a CRISPR
locus looks like to this pipeline, and the explainer's 12-nt consensus is often a
substring of the full CRISPR repeat. Without this layer a reviewer would have to
guess; with it the candidate is either inside an annotated array or it is not.

An overlap here is an ANNOTATION, not a claim that the pipeline discovered CRISPR.
A candidate that carries a genome's CRISPR repeat consensus but falls OUTSIDE every
pilercr array is reported as such, because that is a different and more interesting
situation than either plain case.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

HDR = re.compile(r'^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(?:(\d+)\s+)?([ACGTN-]+)\s*$')


def parse_pilercr(path):
    """Arrays from the SUMMARY BY POSITION block; positions are 1-based."""
    out, contig, inblock = defaultdict(list), None, False
    for line in Path(path).open():
        if 'SUMMARY BY POSITION' in line:
            inblock = True
            continue
        if 'SUMMARY BY SIMILARITY' in line:
            inblock = False
        if not inblock:
            continue
        if line.startswith('>'):
            contig = line[1:].split()[0]
            continue
        m = HDR.match(line.rstrip())
        if m and contig:
            _, _seq, pos, length, ncop, rep, spac, _dist, cons = m.groups()
            out[contig].append({'start0': int(pos) - 1, 'end0': int(pos) - 1 + int(length),
                                'n_copies': int(ncop), 'repeat_len': int(rep),
                                'spacer_len': int(spac), 'consensus': cons})
    return out


def revcomp(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--pilercr', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    arr = parse_pilercr(a.pilercr)
    # every CRISPR repeat consensus seen anywhere in this run, for the off-array check
    cons_all = [c['consensus'].replace('-', '') for v in arr.values() for c in v]

    models = {}
    mp = a.run / 'explain_full/models.tsv'
    if mp.exists():
        for m in csv.DictReader(mp.open(), delimiter='\t'):
            if m.get('primary', '').lower() in ('true', '1', 'yes'):
                models[m['candidate_id']] = m

    rows, n_in, n_off = [], 0, 0
    with (a.run / 'calls/candidates.tsv').open() as fh:
        for c in csv.DictReader(fh, delimiter='\t'):
            s0, e0 = int(c['start0']), int(c['end0'])
            ov = [x for x in arr.get(c['contig'], []) if x['start0'] < e0 and s0 < x['end0']]
            cov = sum(min(e0, x['end0']) - max(s0, x['start0']) for x in ov)
            mcons = (models.get(c['candidate_id'], {}) or {}).get('consensus', '')
            # does the explainer consensus sit inside a known CRISPR repeat?
            matches = ''
            if mcons and len(mcons) >= 8:
                for cs in cons_all:
                    if mcons in cs or revcomp(mcons) in cs:
                        matches = cs
                        break
            if ov:
                n_in += 1
            elif matches:
                n_off += 1
            rows.append({
                'candidate_id': c['candidate_id'], 'contig': c['contig'],
                'start0': s0, 'end0': e0, 'score': c['score'],
                'n_CRISPR_arrays_overlapping': len(ov),
                'CRISPR_covered_nt': cov,
                'CRISPR_coverage_fraction': round(cov / max(1, e0 - s0), 3),
                'CRISPR_array_copies': ','.join(str(x['n_copies']) for x in ov),
                'CRISPR_repeat_len': ','.join(str(x['repeat_len']) for x in ov),
                'CRISPR_spacer_len': ','.join(str(x['spacer_len']) for x in ov),
                'CRISPR_consensus': ','.join(x['consensus'] for x in ov),
                'model_consensus_in_CRISPR_repeat': matches,
                'CRISPR_repeat_without_called_array': bool(matches) and not ov,
            })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    tot = sum(len(v) for v in arr.values())
    print(f'{a.run.name}: {tot} pilercr arrays; {n_in}/{len(rows)} candidates inside one, '
          f'{n_off} carry a CRISPR repeat consensus WITHOUT a called array -> {a.out}')


if __name__ == '__main__':
    main()
