#!/usr/bin/env python3
"""Copy-dependence diagnostic (exact scoring at unchanged target positions only).

DIAGNOSTIC ONLY: its unit definitions and labels are kept apart from the
discovery scorer and the benchmark ground truth.

For each locus, target units (interior copies) and their local flank are left
unchanged. Two perturbation arms, per fixed seed:
  copy    : every other copy of the unit outside the preserved zone is mutated
  noncopy : the same number of equal-length segments in non-copy sequence, at
            matched distance from the target, is mutated identically
Mutation: every base in a segment is replaced by a uniformly chosen different
base. The target positions are then re-scored with exact single-position masking
(both orientations, mapped back to forward coordinates).

  build : make sequences + manifest (CPU)        score : exact scoring (GPU)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_score as rs  # noqa: E402

SEEDS = [11, 23, 37, 41, 53]
FLANK = 20
MAX_TARGETS = 3
MATCH_TOLERANCE = 150   # nt; widened stepwise if no valid placement exists


# ----------------------------------------------------------------- units
def read_fa(path):
    op = gzip.open if str(path).endswith('.gz') else open
    out, name = {}, None
    for line in op(path, 'rt'):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v).lower() for k, v in out.items()}


def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_copies(seq, unit, max_mm):
    unit = unit.lower()
    hits, i = [], 0
    while i <= len(seq) - len(unit):
        if mismatches(seq[i:i + len(unit)], unit) <= max_mm:
            hits.append((i, i + len(unit)))
            i += len(unit)
        else:
            i += 1
    return hits


def locus_units(bench: Path, seqs: dict):
    """Diagnostic unit definitions: {locus: (record_id, flank, [(s,e)...], note)}."""
    man = {r['record_id']: r for r in csv.DictReader(
        (bench / 'data/seed19_manifest.tsv').open(), delimiter='\t')}
    subs = list(csv.DictReader((bench / 'data/annotated_repeat_subregions.tsv').open(),
                               delimiter='\t'))
    loci = {}
    for rid, name in [('w0004', 'SA1_ART'), ('w0005', 'MarsHill_ART')]:
        w0 = int(man[rid]['start0'])
        cores = sorted((int(s['start0']) - w0, int(s['end0']) - w0) for s in subs
                       if s['accession_version'] == man[rid]['accession_version']
                       and 'ART_core' in s['feature'])
        loci[name] = (rid, FLANK, cores, '14-nt ART cores (annotated sub-regions)')
    pos = {r['locus_id']: r for r in csv.DictReader(
        (bench / 'data/positive_loci.tsv').open(), delimiter='\t')}
    unit = pos['LMD9_CRISPR_1377229']['repeat_unit']
    loci['LMD9_CRISPR_1377229'] = ('w0002', FLANK, find_copies(seqs['w0002'], unit, 3),
                                   f'CRISPR repeat {len(unit)} nt, <=3 mismatches')
    # VIPR: 12 contiguous 5-nt YNNGG units of the mapped vrRNA, window [1487,1547)
    loci['SUSP1_VIPR_YNNGG'] = ('w0003', 10, [(1487 + 5 * k, 1492 + 5 * k) for k in range(12)],
                                '12 contiguous 5-nt units; flank reduced to 10 nt')
    # tRNA cluster: GenBank 1-based inclusive -> 0-based half-open window coords
    w0 = 86620
    trna = [(87462, 87536), (87547, 87619), (87629, 87701), (87711, 87782),
            (87793, 87865), (87965, 88036), (88042, 88116), (88298, 88371),
            (88960, 89032), (89038, 89113), (89124, 89196), (89316, 89388)]
    loci['SUSP1_tRNA_cluster'] = ('w0011', FLANK, [(s - 1 - w0, e - w0) for s, e in trna],
                                  '12 GenBank tRNAs (structurally related, not identical)')
    return loci


def pick_targets(units):
    interior = list(range(1, len(units) - 1)) or list(range(len(units)))
    if len(interior) <= MAX_TARGETS:
        return interior
    return sorted({interior[int(round(x))] for x in
                   np.linspace(0, len(interior) - 1, MAX_TARGETS)})


def mutate(seq, spans, rng):
    s = list(seq)
    for a, b in spans:
        for i in range(a, b):
            s[i] = rng.choice([c for c in 'acgt' if c != s[i]])
    return ''.join(s)


def place_noncopy(L, target, others, forbidden, rng):
    """Equal-length segments at matched |distance to target|, avoiding forbidden."""
    taken = forbidden.copy()
    placed, mism = [], []
    tc = (target[0] + target[1]) / 2
    for (a, b) in others:
        ln, d = b - a, abs((a + b) / 2 - tc)
        tol = MATCH_TOLERANCE
        while True:
            cands = [s for s in range(0, L - ln + 1)
                     if not taken[s:s + ln].any() and abs(abs(s + ln / 2 - tc) - d) <= tol]
            if cands or tol > L:
                break
            tol *= 2
        if not cands:
            raise RuntimeError('no valid non-copy placement')
        s = int(rng.choice(cands))
        taken[s:s + ln] = True
        placed.append((s, s + ln))
        mism.append(abs(abs(s + ln / 2 - tc) - d))
    return placed, mism


def build(a):
    bench = a.bench
    seqs = read_fa(bench / 'data/seed19_3kb.fa')
    loci = locus_units(bench, seqs)
    a.out.mkdir(parents=True, exist_ok=True)
    fa = (a.out / 'sequences.fa').open('w')
    rows = []
    for locus, (rid, flank, units, note) in loci.items():
        seq = seqs[rid]
        L = len(seq)
        copy_mask = np.zeros(L, bool)
        for s, e in units:
            copy_mask[s:e] = True
        for ti in pick_targets(units):
            t = units[ti]
            keep = (max(0, t[0] - flank), min(L, t[1] + flank))
            others = [u for j, u in enumerate(units) if j != ti
                      and (u[1] <= keep[0] or u[0] >= keep[1])]
            tpos = list(range(t[0], t[1]))
            base = dict(locus=locus, record_id=rid, target_index=ti,
                        target_start=t[0], target_end=t[1], keep_start=keep[0],
                        keep_end=keep[1], n_other_copies=len(others), note=note)
            sid = f'{locus}|t{ti}|original|0'
            fa.write(f'>{sid}\n{seq}\n')
            rows.append({**base, 'seq_id': sid, 'arm': 'original', 'seed': 0,
                         'mutated_spans': '', 'distance_mismatch_max': 0,
                         'target_positions': tpos})
            for seed in SEEDS:
                rng = np.random.default_rng([seed, ti, int(hashlib.sha256(locus.encode()).hexdigest()[:8], 16)])
                forbidden = copy_mask.copy()
                forbidden[keep[0]:keep[1]] = True
                nc, mism = place_noncopy(L, t, others, forbidden, rng)
                for arm, spans, mm in [('copy', others, [0]), ('noncopy', nc, mism)]:
                    mseq = mutate(seq, spans, rng)
                    # invariants: preserved zone identical; every mutated base changed
                    assert mseq[keep[0]:keep[1]] == seq[keep[0]:keep[1]]
                    assert sum(x != y for x, y in zip(mseq, seq)) == sum(e - s for s, e in spans)
                    sid = f'{locus}|t{ti}|{arm}|{seed}'
                    fa.write(f'>{sid}\n{mseq}\n')
                    rows.append({**base, 'seq_id': sid, 'arm': arm, 'seed': seed,
                                 'mutated_spans': ';'.join(f'{s}-{e}' for s, e in spans),
                                 'distance_mismatch_max': float(max(mm)),
                                 'target_positions': tpos})
    fa.close()
    (a.out / 'design.json').write_text(json.dumps(
        {'seeds': SEEDS, 'flank_default': FLANK, 'max_targets': MAX_TARGETS,
         'match_tolerance_nt': MATCH_TOLERANCE, 'loci': {k: {'record_id': v[0], 'flank': v[1],
         'units': v