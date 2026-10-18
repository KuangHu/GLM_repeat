#!/usr/bin/env python3
"""Draw a uniform random fraction of the unscanned database.

This is deliberately NOT the family-interleaved draw used for rounds 1-3. Those
draws maximised taxonomic breadth, which is right for exploration but makes the
panel unrepresentative of the database as a whole. To decide whether a full
database run is worth its compute, the sample must be an unbiased preview of
that run, so genomes are drawn uniformly at random with a fixed seed and the
complete drawn list is written out.

Selection uses no detector score, no repeat content and no gene content. The only
filters are the ones that define the scannable database: prokaryotic host,
COMPLETE, length in range, and not already scanned.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records', type=Path, required=True)
    ap.add_argument('--exclude-manifest', type=Path, action='append', default=[])
    ap.add_argument('--exclude-fasta', type=Path, action='append', default=[])
    ap.add_argument('--exclusions', type=Path)
    ap.add_argument('--fraction', type=float, default=0.10)
    ap.add_argument('--min-nt', type=int, default=10_000)
    ap.add_argument('--max-nt', type=int, default=350_000)
    ap.add_argument('--seed', type=int, default=20261018)
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
    for f in a.exclude_fasta:
        if f.exists():
            for line in f.open():
                if line.startswith('>'):
                    used.add(line[1:].split()[0].split('.')[0])
    if a.exclusions and a.exclusions.exists():
        with a.exclusions.open() as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                acc = r.get('accession_version') or r.get('accession') or ''
                if acc:
                    used.add(acc.split('.')[0])

    pool = []
    with a.records.open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r.get('accession_version', '')
            if not acc or acc.split('.')[0] in used:
                continue
            n = int(r.get('length_nt', 0) or 0)
            if not (a.min_nt <= n <= a.max_nt):
                continue
            if str(r.get('completeness', '')).upper() != 'COMPLETE':
                continue
            if r.get('host_domain_from_viral_taxonomy') not in (
                    'bacteria', 'archaea', 'bacteria_or_archaea'):
                continue
            pool.append(r)

    pool.sort(key=lambda r: r['accession_version'])      # deterministic order first
    k = max(1, round(len(pool) * a.fraction))
    rng = random.Random(a.seed)
    drawn = rng.sample(pool, k)
    drawn.sort(key=lambda r: r['accession_version'])
    tot = sum(int(r['length_nt']) for r in drawn)
    pool_nt = sum(int(r['length_nt']) for r in pool)

    cols = ['accession_version', 'virus_name', 'family', 'genus', 'branch',
            'length_nt', 'completeness', 'molecule_type_from_taxonomy',
            'host_domain_from_viral_taxonomy', 'release_date', 'virus_taxid']
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w') as fh:
        fh.write('\t'.join(cols + ['selection_reason']) + '\n')
        for r in drawn:
            fh.write('\t'.join(str(r.get(c, '')).replace('\t', ' ') for c in cols)
                     + f'\tuniform random {a.fraction:.0%} of the unscanned database, '
                       f'seed {a.seed}\n')
    fam = Counter(r.get('family') or 'unclassified' for r in drawn)
    pool_fam = Counter(r.get('family') or 'unclassified' for r in pool)
    prov = {
        'draw': 'uniform random without replacement over the unscanned pool',
        'seed': a.seed, 'fraction_requested': a.fraction,
        'pool_genomes': len(pool), 'pool_Mb': round(pool_nt / 1e6, 2),
        'drawn_genomes': len(drawn), 'drawn_Mb': round(tot / 1e6, 2),
        'fraction_achieved_by_genomes': round(len(drawn) / len(pool), 4),
        'fraction_achieved_by_nt': round(tot / pool_nt, 4),
        'families_in_draw': len(fam), 'families_in_pool': len(pool_fam),
        'max_family_share_of_draw': round(max(fam.values()) / len(drawn), 4),
        'why_random_not_interleaved':
            'a diversity-maximising draw would overstate what a full database run '
            'yields; this sample is meant to be an unbiased preview of that run',
        'not_a_yield_claim': 'candidate yields from this panel are still yields, not '
                             'FDR or biological probabilities',
    }
    a.out.with_suffix('.provenance.json').write_text(json.dumps(prov, indent=2) + '\n')
    print(json.dumps(prov, indent=2))


if __name__ == '__main__':
    main()
