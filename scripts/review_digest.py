#!/usr/bin/env python3
"""Compact per-candidate evidence digest for manual review.

Collects, for each browse-book candidate, the features a reviewer actually needs
to decide whether a repeat is real, whether the explainer's model covers it, and
whether it is worth following up -- including the low-complexity measures that
distinguish a structured repeat from a homopolymer or a 2-letter tract.

Computes nothing new about the detector; it reads frozen outputs and measures
sequence composition. No judgement is made here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def tsv(path):
    p = Path(path)
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


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


def entropy(seq):
    if not seq:
        return 0.0
    c = Counter(seq)
    n = len(seq)
    return -sum(v / n * math.log2(v / n) for v in c.values())


def max_homopolymer(seq):
    best = cur = 1
    for i in range(1, len(seq)):
        cur = cur + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, cur)
    return best if seq else 0


def dust_like(seq, k=3):
    """Fraction of k-mer mass in the most frequent k-mers: high = low complexity."""
    if len(seq) < k + 1:
        return 1.0
    c = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    tot = sum(c.values())
    return sum(v for _, v in c.most_common(3)) / tot


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='append', required=True, metavar='LABEL=DIR')
    ap.add_argument('--selection', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    runs = dict(s.partition('=')[::2] for s in a.run)

    data = {}
    for lab, d in runs.items():
        D = Path(d)
        models = defaultdict(list)
        for m in tsv(D / 'explain_full/hypotheses.tsv'):
            models[m['candidate_id']].append(m)
        copies = defaultdict(list)
        for x in tsv(D / 'explain_full/copies.tsv'):
            copies[x['candidate_id']].append(x)
        peaks = defaultdict(set)
        for x in tsv(D / 'calls/peaks.tsv'):
            peaks[x['candidate_id']].add((int(x['peak_start0']), int(x['peak_end0'])))
        probe = {x['candidate_id']: x for x in tsv(D / 'no_model_probe/no_model_probe.tsv')}
        genes = defaultdict(list)
        for g in tsv(D / 'neighborhoods/genes.tsv'):
            genes[g['candidate_id']].append(g)
        prim = {m['candidate_id']: m for m in tsv(D / 'explain_full/models.tsv')
                if m.get('primary', '').lower() in ('true', '1', 'yes')}
        data[lab] = dict(genome=read_fasta(D / 'input.fa'), models=models, copies=copies,
                         peaks=peaks, probe=probe, genes=genes, prim=prim,
                         cand={c['candidate_id']: c for c in tsv(D / 'calls/candidates.tsv')})

    rows = []
    for r in tsv(a.selection):
        lab, cid = r['run'], r['candidate_id']
        D = data[lab]
        c = D['cand'][cid]
        s0, e0 = int(c['start0']), int(c['end0'])
        seq = D['genome'][c['contig']][s0:e0]
        pks = sorted(D['peaks'].get(cid, ()))
        pk_nt = sum(min(e0, b) - max(s0, aa) for aa, b in pks if b > s0 and aa < e0)
        cps = D['copies'].get(cid, [])
        prim = D['prim'].get(cid, {})
        pm = prim.get('model_id', '')
        pcop = [x for x in cps if x['model_id'] == pm]
        inside = [x for x in pcop if x.get('contained_in_candidate', '').lower() == 'true']
        ids = [float(x['loo_consensus_identity']) for x in pcop
               if x.get('loo_consensus_identity') not in (None, '', '.')]
        starts = sorted(int(x['start0']) for x in pcop)
        gaps = [b - a_ for a_, b in zip(starts, starts[1:])]
        cons = prim.get('consensus', '')
        gn = D['genes'].get(cid, [])
        rows.append({
            'run': lab, 'candidate_id': cid, 'strata': r['rules'],
            'in_random_sample': r['in_random_sample'],
            'span_nt': e0 - s0, 'score': round(float(c['score']), 3),
            'n_peaks_dedup': len(pks),
            'peak_coverage_frac': round(pk_nt / max(1, e0 - s0), 3),
            'status': prim.get('status', 'no_model'),
            'architecture': prim.get('architecture', ''),
            'unit_len': prim.get('unit_len', ''),
            'n_models_with_copies': len({x['model_id'] for x in cps}),
            'n_hypotheses': len(D['models'].get(cid, [])),
            'primary_n_copies': len(pcop),
            'primary_copies_inside': len(inside),
            'median_loo_identity': round(sorted(ids)[len(ids) // 2], 3) if ids else '',
            'min_loo_identity': round(min(ids), 3) if ids else '',
            'spacing_cv': prim.get('spacing_cv', ''),
            'median_start_to_start': prim.get('median_start_to_start', ''),
            'start_gaps': ','.join(str(x) for x in gaps[:8]),
            'consensus': cons,
            'consensus_entropy': round(entropy(cons), 2) if cons else '',
            'cand_entropy_bits': round(entropy(seq), 2),
            'cand_max_homopolymer': max_homopolymer(seq),
            'cand_top3_3mer_frac': round(dust_like(seq), 3),
            'gc': round(100 * (seq.count('G') + seq.count('C')) / max(1, len(seq)), 1),
            'probe_verdict': D['probe'].get(cid, {}).get('verdict', ''),
            'probe_best_aln_nt': D['probe'].get(cid, {}).get('best_aln_nt', ''),
            'n_CDS_overlapping': prim.get('n_CDS_overlapping', ''),
            'CDS_overlap_fraction': prim.get('CDS_overlap_fraction', ''),
            'n_neighbor_CDS': len(gn),
            'seq_head': seq[:90],
        })
    rows.sort(key=lambda x: (x['in_random_sample'].lower() != 'true', x['run'],
                             -float(x['score'])))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} candidates digested -> {a.out}')


if __name__ == '__main__':
    main()
