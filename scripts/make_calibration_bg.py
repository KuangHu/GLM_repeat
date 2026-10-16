#!/usr/bin/env python3
"""Sample additional background windows for threshold calibration.

Deterministic, label-blind sampling from the packaged reference genomes. Excludes
every seed19 window and every known positive locus with a margin, so calibration
material is disjoint from the development panel. No gLM2 score, repeat finder,
low-complexity filter or CDS call is used to select or reject a window.
"""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path

import numpy as np


def read_genome(path):
    return ''.join(l.strip() for l in gzip.open(path, 'rt')
                   if not l.startswith('>')).upper()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bench', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--n-per-genome', type=int, default=40)
    p.add_argument('--window', type=int, default=3000)
    p.add_argument('--margin', type=int, default=5000,
                   help='keep this far from seed19 windows and known positives')
    p.add_argument('--seed', type=int, default=20261015)
    a = p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    a.out.mkdir(parents=True, exist_ok=True)

    man = list(csv.DictReader((a.bench / 'data/seed19_manifest.tsv').open(), delimiter='\t'))
    pos = list(csv.DictReader((a.bench / 'data/positive_loci.tsv').open(), delimiter='\t'))
    excl = {}
    for r in man:
        excl.setdefault(r['accession_version'], []).append(
            (int(r['start0']) - a.margin, int(r['end0']) + a.margin))
    for r in pos:
        excl.setdefault(r['accession_version'], []).append(
            (int(r['label_start0']) - a.margin, int(r['label_end0']) + a.margin))

    rng = np.random.default_rng(a.seed)
    rows = []
    with (a.out / 'calibration_bg.fa').open('w') as fa:
        for path in sorted((a.bench / 'genomes').glob('*.fna.gz')):
            acc = path.name.replace('.fna.gz', '')
            g = read_genome(path)
            blocked = np.zeros(len(g), bool)
            for s, e in excl.get(acc, []):
                blocked[max(0, s):min(len(g), e)] = True
            starts = [s for s in range(0, len(g) - a.window + 1, a.window)
                      if not blocked[s:s + a.window].any()
                      and set(g[s:s + a.window]) <= set('ACGT')]
            if not starts:
                print(f'{acc}: no eligible window')
                continue
            take = sorted(rng.choice(starts, size=min(a.n_per_genome, len(starts)),
                                     replace=False).tolist())
            for s in take:
                w = g[s:s + a.window].lower()
                rid = f'cb_{acc}_{s}'
                fa.write(f'>{rid}\n{w}\n')
                rows.append({'record_id': rid, 'family': 'calibration_background',
                             'source_id': acc, 'accession_version': acc,
                             'start0': s, 'end0': s + a.window, 'length_nt': a.window,
                             'gc_fraction': f"{(w.count('g') + w.count('c')) / len(w):.4f}",
                             'cds_fraction': '', 'locus_id': '', 'label_start0': '',
                             'label_end0': '', 'label_start_in_window0': '',
                             'label_end_in_window0': '', 'paired_positive_id': '',
                             'split': 'calibration_only',
                             'evidence': 'unlabeled real DNA; deterministic sampling; '
                                         'no score/repeat/low-complexity filter',
                             'sequence_sha256': hashlib.sha256(w.encode()).hexdigest()})
            print(f'{acc}: {len(take)} of {len(starts)} eligible windows')
    with (a.out / 'calibration_bg_manifest.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} windows ({len(rows) * a.window / 1e3:.0f} kb) -> {a.out}')


if __name__ == '__main__':
    main()
