#!/usr/bin/env python3
"""Select round-3 phage genomes from the frozen viral metadata, excluding all prior rounds.

Families are interleaved rather than drained one at a time. An earlier round used
a naive round-robin that exhausted each stratum's genera before moving on and
produced 77 genomes all from one class; interleaving keeps the panel broad.

Selection uses taxonomy, completeness and length only -- never detector scores,
repeat content or neighbouring genes.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', type=Path, required=True,
                    help='normalized virus_records.tsv')
    ap.add_argument('--exclude-manifest', type=Path, action='append', default=[])
    ap.add_argument('--exclusions', type=Path)
    ap.add_argument('--target-nt', type=int, default=17_000_000)
    ap.add_argument('--min-nt', type=int, default=10_000)
    ap.add_argument('--max-nt', type=int, default=350_000)
    ap.add_argument('--max-family-share', type=float, default=0.20)
    ap.add_argument('--seed', type=int, default=20261017)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()

    used = set()
    for m in a.exclude_manifest:
        if m.exists():
            with m.open() as fh:
                for r in csv.DictReader(fh, delimiter='\t'):
                    acc = r.get('accession_version') or r.get('accession') or ''
                    if acc:
                        used.add(acc.split('.')[0])
    if a.exclusions and a.exclusions.exists():
        with a.exclusions.open() as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                acc = r.get('accession_version') or r.get('accession') or ''
                if acc:
                    used.add(acc.split('.')[0])

    by_family, seen = defaultdict(list), set()
    with a.records.open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r.get('accession_version', '')
            if not acc or acc.split('.')[0] in used or acc in seen:
                continue
            n = int(r.get('length_nt', 0) or 0)
            if not (a.min_nt <= n <= a.max_nt):
                continue
            if str(r.get('completeness', '')).upper() != 'COMPLETE':
                continue
            # prokaryotic hosts only; the eukaryotic-virus contamination found in
            # an earlier round is excluded here by taxonomy, not by name
            if r.get('host_domain_from_viral_taxonomy') not in (
                    'bacteria', 'archaea', 'bacteria_or_archaea'):
                continue
            seen.add(acc)
            by_family[r.get('family') or 'unclassified'].append(r)

    rng = random.Random(a.seed)
    for fam in by_family:
        rng.shuffle(by_family[fam])
    fams = sorted(by_family, key=lambda f: (-len(by_family[f]), f))
    cap = a.target_nt * a.max_family_share
    chosen, total, fam_nt = [], 0, defaultdict(int)
    idx = {f: 0 for f in fams}
    progress = True
    while total < a.target_nt and progress:
        progress = False
        for f in fams:                       # one genome per family per pass
            if total >= a.target_nt:
                break
            i = idx[f]
            while i < len(by_family[f]):
                r = by_family[f][i]
                i += 1
                n = int(r['length_nt'])
                if fam_nt[f] + n > cap and fam_nt[f] > 0:
                    continue
                idx[f] = i
                chosen.append(r)
                fam_nt[f] += n
                total += n
                progress = True
                break
            else:
                idx[f] = i

    cols = ['accession_version', 'virus_name', 'family', 'genus', 'branch',
            'length_nt', 'completeness', 'molecule_type_from_taxonomy',
            'host_domain_from_viral_taxonomy', 'release_date', 'virus_taxid']
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w') as fh:
        fh.write('\t'.join(cols + ['selection_reason']) + '\n')
        for r in chosen:
            fh.write('\t'.join(str(r.get(c, '')).replace('\t', ' ') for c in cols)
                     + '\tnot scanned in round-1 or round-2; family-interleaved draw\n')
    shares = {f: fam_nt[f] / total for f in fam_nt if fam_nt[f]}
    top = max(shares.values()) if shares else 0
    print(f'{len(chosen)} phage genomes, {total:,} nt ({total/1e6:.2f} Mb) '
          f'across {len(shares)} families; max family share {top:.1%}')
    print(f'pool after exclusions: {len(seen)} genomes in {len(by_family)} families')
    print(f'-> {a.out}')


if __name__ == '__main__':
    main()
