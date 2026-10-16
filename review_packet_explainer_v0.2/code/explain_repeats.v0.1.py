#!/usr/bin/env python3
"""Repeat explainer: de novo sequence models for detector candidates.

SEPARATE from the detector. Reads detector outputs read-only and never modifies
them; no gLM2 inference is involved. Discovery is from native sequence alone: no
known motif, family name or annotation is supplied, and search is NOT restricted to
the gLM2 peak intervals (the candidate span plus fixed flanks is searched in full).

Not every candidate gets an explanation. Candidates with no supported model are
emitted as `unresolved`, which is a reported outcome, not a failure to hide.

LIMITATION - indels: all alignment here is UNGAPPED. Copies are matched by
substitution distance only. An array whose copies differ by insertions or deletions
will be under-reported: its copies may be split, truncated, or missed entirely, and
reported spacer lengths will absorb the indel. Period estimates from the self-match
profile are likewise substitution-based. This is recorded per model as
`indel_handling: ungapped_substitutions_only`.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

BASES = 'ACGT'


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


def self_match_profile(seq, min_p=2, max_p=600):
    """Fraction of identical bases at offset d, for each candidate period d."""
    a = np.frombuffer(seq.encode(), dtype='S1')
    out = {}
    for d in range(min_p, min(max_p, len(a) - 20) + 1):
        out[d] = float((a[:-d] == a[d:]).mean())
    return out


def low_complexity(seq):
    """Dominant mono/di-nucleotide fraction, for the homopolymer alternative."""
    if not seq:
        return 0.0, ''
    mono = Counter(seq).most_common(1)[0]
    di = Counter(seq[i:i + 2] for i in range(len(seq) - 1)).most_common(1)[0]
    return max(mono[1] / len(seq), di[1] / max(1, len(seq) - 1)), \
        mono[0] if mono[1] / len(seq) >= di[1] / max(1, len(seq) - 1) else di[0]


def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_copies(seq, motif, max_mm):
    """Non-overlapping ungapped occurrences within max_mm substitutions."""
    k, hits, i = len(motif), [], 0
    while i <= len(seq) - k:
        mm = mismatches(seq[i:i + k], motif)
        if mm <= max_mm:
            hits.append((i, i + k, mm))
            i += k
        else:
            i += 1
    return hits


def pwm_from(copies):
    if not copies:
        return None
    k = len(copies[0])
    m = np.zeros((k, 4))
    for c in copies:
        for j, ch in enumerate(c):
            if ch in BASES:
                m[j, BASES.index(ch)] += 1
    m += 0.25                      # Laplace, keeps probabilities finite
    return m / m.sum(1, keepdims=True)


def consensus(pwm):
    return ''.join(BASES[i] for i in pwm.argmax(1))


def information(pwm):
    return float((2 + (pwm * np.log2(pwm)).sum(1)).mean())


def extend(seq, hits, max_ext=60, min_conserve=0.75):
    """Grow the motif outward while columns stay conserved across copies."""
    s0 = [h[0] for h in hits]
    k = hits[0][1] - hits[0][0]
    left = 0
    while left < max_ext and all(h - left - 1 >= 0 for h in s0):
        col = [seq[h - left - 1] for h in s0]
        if Counter(col).most_common(1)[0][1] / len(col) < min_conserve:
            break
        left += 1
    right = 0
    while right < max_ext and all(h + k + right + 1 <= len(seq) for h in s0):
        col = [seq[h + k + right] for h in s0]
        if Counter(col).most_common(1)[0][1] / len(col) < min_conserve:
            break
        right += 1
    return left, right


def build_model(seq, min_copies=2, seed_k=12):
    """Seed on an over-represented k-mer, extend, then rescan with the PWM."""
    counts = Counter(seq[i:i + seed_k] for i in range(len(seq) - seed_k + 1))
    cands = [(kmer, n) for kmer, n in counts.items()
             if n >= min_copies and len(set(kmer)) > 2]
    if not cands:
        return None
    models = []
    for kmer, _ in sorted(cands, key=lambda kv: -kv[1])[:6]:
        hits = find_copies(seq, kmer, max_mm=0)
        if len(hits) < min_copies:
            continue
        left, right = extend(seq, hits)
        s = hits[0][0] - left
        motif = seq[s:s + seed_k + left + right]
        if len(motif) < 5:
            continue
        mm = max(1, int(round(0.2 * len(motif))))
        copies = find_copies(seq, motif, max_mm=mm)
        if len(copies) < min_copies:
            continue
        strs = [seq[a:b] for a, b, _ in copies]
        pwm = pwm_from(strs)
        models.append({'motif': consensus(pwm), 'unit_len': len(motif),
                       'copies': copies, 'copy_strings': strs, 'pwm': pwm,
                       'information_bits': information(pwm),
                       'score': len(copies) * len(motif) * information(pwm)})
    if not models:
        return None
    models.sort(key=lambda m: -m['score'])
    return models


def classify(model, seq_len):
    c = model['copies']
    starts = [a for a, _, _ in c]
    s2s = [starts[i + 1] - starts[i] for i in range(len(starts) - 1)]
    edges = [c[i + 1][0] - c[i][1] for i in range(len(c) - 1)]
    k = model['unit_len']
    if not s2s:
        return 'single_copy', s2s, edges
    cv = float(np.std(s2s) / np.mean(s2s)) if np.mean(s2s) else 1.0
    med_edge = float(np.median(edges))
    if k <= 10 and med_edge <= 2:
        kind = 'short_period_tandem'
    elif med_edge <= 2:
        kind = 'tandem_repeat'
    elif cv <= 0.15:
        kind = 'regular_spaced_array'
    elif cv <= 0.5:
        kind = 'variable_spacer_array'
    else:
        kind = 'dispersed_repeat'
    return kind, s2s, edges


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates', type=Path, required=True)
    p.add_argument('--genomes', type=Path, nargs='+', required=True,
                   help='FASTA dirs or files providing native sequence')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--flank', type=int, default=500)
    p.add_argument('--peaks', type=Path, help='peaks TSV for peak-to-copy correspondence')
    p.add_argument('--min-copies', type=int, default=2)
    p.add_argument('--label', default='candidates')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    seqs = {}
    for g in a.genomes:
        for f in ([g] if g.is_file() else sorted(g.glob('*.fna'))):
            seqs.update(read_fasta(f))
    cands = list(csv.DictReader(a.candidates.open(), delimiter='\t'))
    peaks = {}
    if a.peaks and a.peaks.exists():
        for r in csv.DictReader(a.peaks.open(), delimiter='\t'):
            peaks.setdefault(r['contig'], []).append((int(r['start0']), int(r['end0'])))

    models_rows, copy_rows, aln_rows, unresolved = [], [], [], []
    for i, c in enumerate(cands, 1):
        contig = c['contig']
        if contig not in seqs:
            unresolved.append({**c, 'reason': 'no sequence for contig'})
            continue
        g = seqs[contig]
        cs, ce = int(c['start0']), int(c['end0'])
        s = max(0, cs - a.flank)
        e = min(len(g), ce + a.flank)
        sub = g[s:e]
        cid = f'{contig}:{cs}-{ce}'
        if set(sub) - set(BASES):
            sub = ''.join(ch if ch in BASES else 'A' for ch in sub)
        found = build_model(sub, a.min_copies)
        span = sub[cs - s:ce - s]
        lc_frac, lc_unit = low_complexity(span)
        alts = []
        if lc_frac >= 0.6:
            alts.append({'type': 'low_complexity', 'unit': lc_unit,
                         'dominant_fraction': round(lc_frac, 3)})
        # Periodicity is evaluated on the candidate span as well as the flanked
        # region: a short candidate's period is swamped if only the flanked
        # region is profiled.
        for scope, s_ in (('candidate_span', span), ('span_plus_flanks', sub)):
            prof = self_match_profile(s_, max_p=min(600, max(2, len(s_) // 3)))
            if not prof:
                continue
            bp = max(prof, key=prof.get)
            if prof[bp] >= 0.5:
                alts.append({'type': 'periodic_self_match', 'scope': scope,
                             'period': bp, 'identity_at_period': round(prof[bp], 3)})
        if not found:
            unresolved.append({'candidate_id': cid, 'contig': contig, 'start0': cs,
                               'end0': ce, 'score': c.get('score', ''),
                               'reason': 'no motif with >=%d ungapped copies' % a.min_copies,
                               'alternatives': json.dumps(alts)})
            continue
        m = found[0]
        kind, s2s, edges = classify(m, len(sub))
        gstart = s
        cov = sum(b - aa for aa, b, _ in m['copies'])
        inside = sum(1 for aa, _, _ in m['copies'] if cs <= gstart + aa < ce)
        pk = peaks.get(contig, [])
        per_copy_peaks = []
        for aa, bb, mm in m['copies']:
            g0, g1 = gstart + aa, gstart + bb
            n = sum(1 for ps, pe in pk if ps < g1 and g0 < pe)
            per_copy_peaks.append(n)
            copy_rows.append({'candidate_id': cid, 'contig': contig, 'start0': g0,
                              'end0': g1, 'mismatches_to_consensus': mm,
                              'inside_candidate_span': cs <= g0 < ce,
                              'overlapping_peaks': n})
        for j, cstr in enumerate(m['copy_strings'], 1):
            aln_rows.append({'candidate_id': cid, 'copy_index': j, 'alignment': cstr})
        models_rows.append({
            'candidate_id': cid, 'contig': contig, 'cand_start0': cs, 'cand_end0': ce,
            'detector_score': c.get('score', ''),
            'annotation_category': c.get('annotation_category', ''),
            'architecture': kind, 'consensus': m['motif'], 'unit_len': m['unit_len'],
            'n_copies': len(m['copies']),
            'n_copies_inside_candidate': inside,
            'n_copies_in_flank_only': len(m['copies']) - inside,
            'copy_coverage_nt': cov,
            'information_bits_per_position': round(m['information_bits'], 3),
            'median_start_to_start': int(np.median(s2s)) if s2s else '',
            'start_to_start_range': f'{min(s2s)}-{max(s2s)}' if s2s else '',
            'median_spacer_edge_gap': int(np.median(edges)) if edges else '',
            'spacer_cv': round(float(np.std(s2s) / np.mean(s2s)), 3) if s2s and np.mean(s2s) else '',
            'copies_with_>=1_peak': sum(1 for n in per_copy_peaks if n),
            'alternatives': json.dumps(alts),
            'n_alternative_models': max(0, len(found) - 1),
            'indel_handling': 'ungapped_substitutions_only'})
        np.save(a.out / f'pwm_{i:04d}.npy', m['pwm'])
        models_rows[-1]['pwm_file'] = f'pwm_{i:04d}.npy'

    for name, rows in [('models.tsv', models_rows), ('copies.tsv', copy_rows),
                       ('alignments.tsv', aln_rows), ('unresolved.tsv', unresolved)]:
        if rows:
            with (a.out / name).open('w', newline='') as fh:
                w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
                w.writeheader()
                w.writerows(rows)
    arch = Counter(r['architecture'] for r in models_rows)
    (a.out / 'summary.json').write_text(json.dumps({
        'label': a.label, 'n_candidates': len(cands),
        'n_explained': len(models_rows), 'n_unresolved': len(unresolved),
        'architectures': dict(arch), 'flank_nt': a.flank,
        'search_scope': 'candidate span plus flanks, not restricted to peak intervals',
        'indel_limitation': 'ungapped substitution alignment only; indel-containing '
                            'arrays are under-reported, copies may be split or missed '
                            'and spacer lengths absorb the indel',
        'detector_outputs': 'read-only; not modified'}, indent=2) + '\n')
    print(f'{len(cands)} candidates: {len(models_rows)} explained, '
          f'{len(unresolved)} unresolved')
    for k, v in arch.most_common():
        print(f'  {k:26s} {v}')


if __name__ == '__main__':
    main()
