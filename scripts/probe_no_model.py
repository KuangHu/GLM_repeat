#!/usr/bin/env python3
"""Indel-aware repeat probe for candidates the explainer could not model.

The v0.3 explainer aligns UNGAPPED, so a `no_model` candidate may still contain a
real repeat whose copies differ by insertions or deletions. This probe seeds on
short exact k-mers and extends each seed pair with a GAPPED local alignment
(Biopython PairwiseAligner, affine gaps), which tolerates indels.

Diagnostic only. It does not modify the detector, the explainer, or any frozen
output, and its result is never fed back into candidate selection. Its purpose is
to separate:

  explainer coverage gap : repeated seeds exist and extend into a gapped alignment
                           that the ungapped explainer cannot represent
  no repeat evidence     : no repeated seed, or seeds that do not extend

METHOD VALIDATION IS MANDATORY AND BUILT IN. `--positive-control` runs the same
probe over candidates the explainer DID model and reports their detection rate.
An earlier version of this probe used minimap2 -x asm20 and reported "no internal
similarity" for every no-model candidate; the positive control showed it also
found nothing for candidates with 11-34 verified copies, because minimizer seeding
cannot seed on 8-33 nt repeat units. Never trust a negative from this probe
without its positive-control line.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
from pathlib import Path

from Bio import Align


def read_fasta(path):
    out, name = {}, None
    for line in Path(path).open():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v).upper() for k, v in out.items()}


def make_aligner():
    a = Align.PairwiseAligner(mode='local', match_score=2, mismatch_score=-3,
                              open_gap_score=-5, extend_gap_score=-2)
    return a


def probe(seq, aligner, k=10, flank=40, max_pairs=40, min_aln=20):
    pos = collections.defaultdict(list)
    for i in range(len(seq) - k + 1):
        pos[seq[i:i + k]].append(i)
    rep = {m: p for m, p in pos.items() if len(p) > 1 and len(set(m)) > 2}
    if not rep:
        return {'repeated_seeds': 0, 'best_aln_nt': 0, 'best_has_gap': False,
                'verdict': 'no_repeated_seed'}
    pairs = []
    for m, p in rep.items():
        for i in range(len(p) - 1):
            if p[i + 1] - p[i] >= k:
                pairs.append((p[i], p[i + 1]))
    pairs = sorted(set(pairs))[:max_pairs]
    best_len, best_gap = 0, False
    for a, b in pairs:
        s1 = seq[max(0, a - flank):a + flank]
        s2 = seq[max(0, b - flank):b + flank]
        try:
            aln = aligner.align(s1, s2)[0]
        except Exception:
            continue
        ln = aln.score / 2
        if ln > best_len:
            try:
                rows = str(aln).splitlines()
                gap = any('-' in r for r in rows[::2][:2])
            except Exception:
                gap = False
            best_len, best_gap = ln, gap
    if best_len < min_aln:
        verdict = 'seeds_without_extension'
    elif best_gap:
        verdict = 'gapped_similarity_explainer_cannot_represent'
    else:
        verdict = 'ungapped_similarity'
    return {'repeated_seeds': len(rep), 'best_aln_nt': round(best_len, 1),
            'best_has_gap': best_gap, 'verdict': verdict}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates', type=Path, required=True)
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--genome', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--flank', type=int, default=500)
    p.add_argument('--k', type=int, default=10)
    p.add_argument('--min-aln', type=int, default=20)
    p.add_argument('--positive-control', type=int, default=20)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(a.genome)
    models = list(csv.DictReader(a.models.open(), delimiter='\t'))
    modelled = {m['candidate_id'] for m in models}
    cand = list(csv.DictReader(a.candidates.open(), delimiter='\t'))
    aligner = make_aligner()

    def region(contig, s0, e0):
        g = seqs[contig]
        return g[max(0, int(s0) - a.flank):min(len(g), int(e0) + a.flank)]

    rows = []
    for c in cand:
        cid = f"{c['contig']}:{c['start0']}-{c['end0']}"
        if cid in modelled:
            continue
        r = probe(region(c['contig'], c['start0'], c['end0']), aligner, a.k,
                  min_aln=a.min_aln)
        rows.append({'candidate_id': cid, 'contig': c['contig'],
                     'start0': c['start0'], 'end0': c['end0'], 'score': c['score'], **r})
    with (a.out / 'no_model_probe.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # ---- mandatory positive control on candidates the explainer DID model
    pc = [m for m in models if m.get('status') == 'candidate_associated'
          and int(m.get('n_copies', 0) or 0) >= 5]
    pc.sort(key=lambda m: -int(m['n_copies']))
    pc = pc[:a.positive_control]
    pc_rows = []
    for m in pc:
        r = probe(region(m['contig'], m['cand_start0'], m['cand_end0']), aligner, a.k,
                  min_aln=a.min_aln)
        pc_rows.append({'candidate_id': m['candidate_id'], 'n_copies': m['n_copies'],
                        'unit_len': m['unit_len'], **r})
    if pc_rows:
        with (a.out / 'positive_control.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, list(pc_rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(pc_rows)
    detected = sum(1 for r in pc_rows if r['verdict'] not in
                   ('no_repeated_seed', 'seeds_without_extension'))
    v = collections.Counter(r['verdict'] for r in rows)
    sens = detected / len(pc_rows) if pc_rows else None
    (a.out / 'summary.json').write_text(json.dumps({
        'n_no_model_candidates': len(rows), 'verdicts': dict(v),
        'positive_control_n': len(pc_rows), 'positive_control_detected': detected,
        'positive_control_sensitivity': sens,
        'method': f'exact {a.k}-mer seeds + Biopython local affine-gap extension',
        'diagnostic_only': 'does not modify detector, explainer or frozen outputs',
        'interpretation_rule': 'a negative result is only meaningful if '
                               'positive_control_sensitivity is high'}, indent=2) + '\n')
    print(f'{len(rows)} no-model candidates probed')
    for k, n in v.most_common():
        print(f'  {k:46s} {n}')
    if pc_rows:
        print(f'POSITIVE CONTROL: {detected}/{len(pc_rows)} modelled candidates detected '
              f'(sensitivity {sens:.0%})')


if __name__ == '__main__':
    main()
