#!/usr/bin/env python3
"""Attach tRNA annotation to candidates.

pyrodigal predicts protein-coding genes only, so the candidate pages previously
carried no annotation for structured non-coding RNA. Many high-scoring candidates
sit on tRNA genes, which a CDS-only annotation layer makes invisible. This parses
aragorn output and records, per candidate, which annotated tRNAs it overlaps.

A tRNA overlap is an ANNOTATION, not a judgement: it says the detector fired on a
known structured ncRNA gene, which is a different thing from a repeat.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

ROW = re.compile(r'^\s*\d+\s+(\S+)\s+(c?)\[(-?\d+),(\d+)\]')


def parse_aragorn(path):
    out = defaultdict(list)
    contig = None
    for line in Path(path).open():
        if line.startswith('>'):
            contig = line[1:].split()[0]
            continue
        m = ROW.match(line)
        if m and contig:
            name, comp, s, e = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
            out[contig].append({'name': name, 'strand': '-' if comp else '+',
                                'start0': max(0, s - 1), 'end0': e})
    return out


def tsv(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--aragorn', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    ann = parse_aragorn(a.aragorn)
    rows, hits = [], 0
    for c in tsv(a.run / 'calls/candidates.tsv'):
        s0, e0 = int(c['start0']), int(c['end0'])
        ov = [t for t in ann.get(c['contig'], [])
              if t['start0'] < e0 and s0 < t['end0']]
        cov = sum(min(e0, t['end0']) - max(s0, t['start0']) for t in ov)
        if ov:
            hits += 1
        rows.append({
            'candidate_id': c['candidate_id'], 'contig': c['contig'],
            'start0': s0, 'end0': e0, 'span_nt': e0 - s0, 'score': c['score'],
            'n_tRNA_overlapping': len(ov),
            'tRNA_covered_nt': cov,
            'tRNA_coverage_fraction': round(cov / max(1, e0 - s0), 3),
            'tRNA_names': ','.join(t['name'] for t in ov),
            'tRNA_coords': ','.join(f"{t['start0']}-{t['end0']}{t['strand']}" for t in ov),
        })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    tot = sum(len(v) for v in ann.values())
    print(f'{a.run.name}: {tot} annotated tRNAs; {hits}/{len(rows)} candidates '
          f'({100*hits/len(rows):.1f}%) overlap at least one -> {a.out}')


if __name__ == '__main__':
    main()
