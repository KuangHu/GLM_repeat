#!/usr/bin/env python3
"""Frozen-parameter validation of the no-model repeat probe.

Measures probe sensitivity on a SYNTHETIC panel with known ground truth, stratified
by the three properties that matter for the indel-coverage question:

  unit length   8, 12, 20, 40, 80 nt
  degeneracy    0%, 5%, 10%, 20% per-base substitution in each copy
  indel type    none | 1 nt indel per copy | 3 nt indel per copy | mixed

Real modelled candidates cannot answer this: their unit lengths and divergence are
whatever the data happens to contain, and their indel content is unknown. A
synthetic panel with planted repeats is the only way to state what fraction of
indel-containing repeats the probe would recover.

Background is real genomic DNA (shuffled windows from the scanned panel) so base
composition is realistic; the planted repeat is the only known signal.

This validates the PROBE, not the detector or the explainer.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe_no_model import make_aligner, probe, read_fasta  # noqa: E402

UNIT_LENS = [8, 12, 20, 40, 80]
DIVERGENCE = [0.0, 0.05, 0.10, 0.20]
INDELS = ['none', 'indel1', 'indel3', 'mixed']
COPIES = 4
SPACER = 120


def mutate(unit, div, rng):
    return ''.join(rng.choice([c for c in 'ACGT' if c != b]) if rng.random() < div else b
                   for b in unit)


def apply_indel(seq, kind, rng):
    if kind == 'none':
        return seq
    n = 1 if kind == 'indel1' else 3 if kind == 'indel3' else rng.choice([1, 2, 3])
    if rng.random() < 0.5 and len(seq) > n + 2:          # deletion
        i = rng.randrange(1, len(seq) - n)
        return seq[:i] + seq[i + n:]
    i = rng.randrange(1, len(seq))                        # insertion
    return seq[:i] + ''.join(rng.choice('ACGT') for _ in range(n)) + seq[i:]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--background', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--replicates', type=int, default=20)
    p.add_argument('--seed', type=int, default=20261016)
    p.add_argument('--k', type=int, default=10)
    p.add_argument('--min-aln', type=int, default=20)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(a.background)
    pool = [s for s in seqs.values() if len(s) > 5000]
    rng = random.Random(a.seed)
    aligner = make_aligner()

    rows = []
    for ul in UNIT_LENS:
        for div in DIVERGENCE:
            for ind in INDELS:
                det = 0
                gapped_calls = 0
                for rep in range(a.replicates):
                    src = rng.choice(pool)
                    off = rng.randrange(0, len(src) - 2000)
                    bg = list(src[off:off + 2000])
                    rng.shuffle(bg)                       # destroy native repeats
                    bg = ''.join(bg)
                    unit = ''.join(rng.choice('ACGT') for _ in range(ul))
                    block = ''
                    for _ in range(COPIES):
                        c = apply_indel(mutate(unit, div, rng), ind, rng)
                        block += c + ''.join(rng.choice('ACGT') for _ in range(SPACER))
                    ins = rng.randrange(200, len(bg) - len(block) - 200)
                    seq = bg[:ins] + block + bg[ins + len(block):]
                    r = probe(seq, aligner, a.k, min_aln=a.min_aln)
                    if r['verdict'] in ('ungapped_similarity',
                                        'gapped_similarity_explainer_cannot_represent'):
                        det += 1
                    if r['verdict'] == 'gapped_similarity_explainer_cannot_represent':
                        gapped_calls += 1
                rows.append({'unit_len': ul, 'divergence': div, 'indel': ind,
                             'replicates': a.replicates, 'detected': det,
                             'sensitivity': round(det / a.replicates, 3),
                             'called_gapped': gapped_calls})
    with (a.out / 'sensitivity.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    by_unit = defaultdict(lambda: [0, 0])
    by_indel = defaultdict(lambda: [0, 0])
    by_div = defaultdict(lambda: [0, 0])
    for r in rows:
        for d, key in ((by_unit, r['unit_len']), (by_indel, r['indel']), (by_div, r['divergence'])):
            d[key][0] += r['detected']
            d[key][1] += r['replicates']
    agg = lambda d: {str(k): {'detected': v[0], 'total': v[1],
                              'sensitivity': round(v[0] / v[1], 3)} for k, v in sorted(d.items())}
    out = {
        'frozen_probe_parameters': {
            'seed_kmer_k': a.k, 'seed_requires_exact_match': True,
            'seed_min_distinct_bases': 3, 'min_seed_separation_nt': a.k,
            'extension': 'Biopython PairwiseAligner, mode=local',
            'match_score': 2, 'mismatch_score': -3,
            'open_gap_score': -5, 'extend_gap_score': -2,
            'extension_flank_nt': 40, 'max_seed_pairs': 40,
            'min_alignment_nt': a.min_aln,
            'alignment_length_definition': 'alignment score / match_score',
            'interval_overlap_handling': 'seed pairs must be >= k apart; the main '
                                         'self-diagonal is excluded by construction '
                                         'since only distinct seed occurrences are paired',
            'identity_reported': 'not thresholded; extension score is the criterion',
        },
        'synthetic_panel': {'unit_lengths': UNIT_LENS, 'divergence': DIVERGENCE,
                            'indel_types': INDELS, 'copies_per_repeat': COPIES,
                            'spacer_nt': SPACER, 'replicates_per_cell': a.replicates,
                            'background': 'shuffled real genomic windows',
                            'seed': a.seed},
        'sensitivity_by_unit_length': agg(by_unit),
        'sensitivity_by_indel_type': agg(by_indel),
        'sensitivity_by_divergence': agg(by_div),
        'interpretation_rule': (
            'A negative from this probe means: no qualifying local similarity was '
            'detected under these parameters. It does NOT mean the region contains '
            'no repeat. Sensitivity is strongly dependent on unit length and '
            'divergence, so a negative is only as informative as the matching cell '
            'of this table.'),
    }
    (a.out / 'validation.json').write_text(json.dumps(out, indent=2) + '\n')
    print('sensitivity by unit length :', {k: v['sensitivity'] for k, v in out['sensitivity_by_unit_length'].items()})
    print('sensitivity by divergence  :', {k: v['sensitivity'] for k, v in out['sensitivity_by_divergence'].items()})
    print('sensitivity by indel type  :', {k: v['sensitivity'] for k, v in out['sensitivity_by_indel_type'].items()})


if __name__ == '__main__':
    main()
