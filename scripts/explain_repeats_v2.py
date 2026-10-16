#!/usr/bin/env python3
"""Repeat explainer v0.2. v0.1 preserved at provenance/explain_repeats.v0.1.py.

Separate from the detector: detector outputs are read-only and no gLM2 inference is
used. Discovery is from native sequence only, over the candidate span plus fixed
flanks, never restricted to gLM2 peak intervals. Candidates with no supported model
remain `unresolved`; nothing is forced.

Two changes from v0.1.

1. PERIODIC branch is a first-class model source, emitting the same schema as the
   motif branch. It needs no exact k-mer seed: it infers interval, period and phase
   from the self-match profile, then emits per-unit coordinates, sequences,
   alignment, position probabilities, per-unit support and explicit boundary/phase
   uncertainty. No period, motif or locus boundary is hard-coded; the period is
   searched over a range and chosen by the data.

2. CORE vs EXTENDED repeat evidence is kept separate. The seed-supported core
   motif and its copies are preserved even when flank extension fails to match, and
   core-copy and extended-copy counts are reported separately. Any core copy that
   the extended motif rejects is recorded with the evidence for its rejection. No
   motif sequence or core length is hard-coded.

LIMITATION - indels: all alignment is UNGAPPED, substitution distance only. Arrays
whose copies differ by insertions or deletions are under-reported; copies may be
split, truncated or missed, and spacer lengths absorb the indel. Period estimates
are likewise substitution-based. Recorded per model as
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
MIN_PERIOD, MAX_PERIOD = 2, 300
PERIODIC_MIN_IDENTITY = 0.55     # local identity at the period to call periodicity
PERIODIC_MIN_UNITS = 3


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


def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b))


def low_complexity(seq):
    if not seq:
        return 0.0, ''
    mono = Counter(seq).most_common(1)[0]
    di = Counter(seq[i:i + 2] for i in range(len(seq) - 1)).most_common(1)[0]
    f_mono = mono[1] / len(seq)
    f_di = di[1] / max(1, len(seq) - 1)
    return (f_mono, mono[0]) if f_mono >= f_di else (f_di, di[0])


def pwm_from(copies):
    if not copies:
        return None
    k = max(len(c) for c in copies)
    m = np.zeros((k, 4))
    for c in copies:
        for j, ch in enumerate(c):
            if ch in BASES:
                m[j, BASES.index(ch)] += 1
    m += 0.25
    return m / m.sum(1, keepdims=True)


def consensus(pwm):
    return ''.join(BASES[i] for i in pwm.argmax(1))


def information(pwm):
    return float((2 + (pwm * np.log2(pwm)).sum(1)).mean())


# ----------------------------------------------------------------- periodic branch
def period_identity(seq, d):
    a = np.frombuffer(seq.encode(), dtype='S1')
    if len(a) <= d:
        return 0.0, np.zeros(0, bool)
    eq = a[:-d] == a[d:]
    return float(eq.mean()), eq


def best_period(seq):
    """Period maximising identity at offset d, searched over a range (not fixed)."""
    best = (0.0, None)
    hi = min(MAX_PERIOD, max(MIN_PERIOD, len(seq) // 3))
    for d in range(MIN_PERIOD, hi + 1):
        ident, _ = period_identity(seq, d)
        if ident > best[0]:
            best = (ident, d)
    return best[1], best[0]


def periodic_interval(seq, d, win=None):
    """Maximal run where local identity at offset d stays above threshold."""
    ident, eq = period_identity(seq, d)
    if eq.size == 0:
        return None
    win = win or max(2 * d, 20)
    k = min(win, eq.size)
    c = np.cumsum(np.insert(eq.astype(float), 0, 0))
    local = (c[k:] - c[:-k]) / k
    good = local >= PERIODIC_MIN_IDENTITY
    if not good.any():
        return None
    runs, s = [], None
    for i, b in enumerate(good):
        if b and s is None:
            s = i
        elif not b and s is not None:
            runs.append((s, i + k))
            s = None
    if s is not None:
        runs.append((s, len(good) + k))
    s, e = max(runs, key=lambda r: r[1] - r[0])
    return max(0, s), min(len(seq), e + d)


def periodic_model(seq, lo, hi, d):
    """Choose phase by PWM information, emit units, report phase/boundary margin."""
    region = seq[lo:hi]
    n = (hi - lo) // d
    if n < PERIODIC_MIN_UNITS:
        return None
    scored = []
    for phase in range(d):
        units = [region[phase + i * d: phase + (i + 1) * d]
                 for i in range((len(region) - phase) // d)]
        units = [u for u in units if len(u) == d]
        if len(units) < PERIODIC_MIN_UNITS:
            continue
        pwm = pwm_from(units)
        scored.append((information(pwm), phase, units, pwm))
    if not scored:
        return None
    scored.sort(key=lambda x: -x[0])
    info, phase, units, pwm = scored[0]
    phase_margin = round(info - scored[1][0], 4) if len(scored) > 1 else None
    cons = consensus(pwm)
    starts = [lo + phase + i * d for i in range(len(units))]
    support = [round(1 - mismatches(u, cons) / d, 3) for u in units]
    # boundary uncertainty: identity at the period just outside the called interval
    out_left = out_right = None
    if lo - d >= 0:
        out_left = round(float(np.mean([seq[lo - d + i] == seq[lo + i]
                                        for i in range(min(d, len(seq) - lo))])), 3)
    if hi + d <= len(seq):
        out_right = round(float(np.mean([seq[hi - d + i] == seq[hi + i]
                                         for i in range(min(d, len(seq) - hi))])), 3)
    return {'source': 'periodic', 'motif': cons, 'unit_len': d,
            'copies': [(s, s + d, mismatches(seq[s:s + d], cons)) for s in starts],
            'copy_strings': units, 'pwm': pwm, 'information_bits': info,
            'per_unit_support': support,
            'interval': (lo, hi), 'phase': phase, 'phase_info_margin': phase_margin,
            'boundary_identity_left': out_left, 'boundary_identity_right': out_right,
            'core_copies': len(units), 'extended_copies': None,
            'rejected_cores': []}


# -------------------------------------------------------------------- motif branch
def find_copies(seq, motif, max_mm):
    k, hits, i = len(motif), [], 0
    while i <= len(seq) - k:
        mm = mismatches(seq[i:i + k], motif)
        if mm <= max_mm:
            hits.append((i, i + k, mm))
            i += k
        else:
            i += 1
    return hits


def extension_offsets(seq, starts, k, max_ext=60, min_conserve=0.75):
    left = 0
    while left < max_ext and all(h - left - 1 >= 0 for h in starts):
        col = [seq[h - left - 1] for h in starts]
        if Counter(col).most_common(1)[0][1] / len(col) < min_conserve:
            break
        left += 1
    right = 0
    while right < max_ext and all(h + k + right + 1 <= len(seq) for h in starts):
        col = [seq[h + k + right] for h in starts]
        if Counter(col).most_common(1)[0][1] / len(col) < min_conserve:
            break
        right += 1
    return left, right


def motif_models(seq, min_copies=2, seed_k=12):
    """Core evidence first; extension is optional and never discards core copies."""
    counts = Counter(seq[i:i + seed_k] for i in range(len(seq) - seed_k + 1))
    cands = [(kmer, n) for kmer, n in counts.items()
             if n >= min_copies and len(set(kmer)) > 2]
    if not cands:
        return []
    models = []
    for kmer, _ in sorted(cands, key=lambda kv: -kv[1])[:6]:
        seeds = find_copies(seq, kmer, max_mm=0)
        if len(seeds) < min_copies:
            continue
        # --- CORE: the seed motif, with a mismatch budget scaled to its length
        core_mm = max(1, int(round(0.2 * seed_k)))
        core = find_copies(seq, kmer, max_mm=core_mm)
        if len(core) < min_copies:
            continue
        core_pwm = pwm_from([seq[a:b] for a, b, _ in core])
        core_cons = consensus(core_pwm)
        core = find_copies(seq, core_cons, max_mm=core_mm)
        # --- EXTENDED: optional, evaluated separately, never replaces the core set
        left, right = extension_offsets(seq, [a for a, _, _ in core], seed_k)
        ext_copies, rejected, ext_cons, ext_len = [], [], '', seed_k
        if left or right:
            s0 = core[0][0] - left
            ext_motif = seq[s0:s0 + seed_k + left + right]
            ext_len = len(ext_motif)
            ext_mm = max(1, int(round(0.2 * ext_len)))
            ext_copies = find_copies(seq, ext_motif, max_mm=ext_mm)
            ext_pwm = pwm_from([seq[a:b] for a, b, _ in ext_copies]) if ext_copies else None
            ext_cons = consensus(ext_pwm) if ext_pwm is not None else ''
            ext_starts = {a for a, _, _ in ext_copies}
            for a, b, mm in core:
                if any(abs(a - left - x) <= 2 for x in ext_starts):
                    continue
                cand = seq[max(0, a - left): a - left + ext_len]
                rejected.append({'core_start': a, 'core_mismatches': mm,
                                 'extended_mismatches': (mismatches(cand, ext_cons)
                                                         if ext_cons and len(cand) == ext_len
                                                         else None),
                                 'extended_budget': ext_mm,
                                 'reason': 'flank did not match extended consensus'
                                           if len(cand) == ext_len else
                                           'insufficient sequence for extension'})
        use_ext = len(ext_copies) >= len(core)
        chosen = ext_copies if use_ext else core
        strs = [seq[a:b] for a, b, _ in chosen]
        pwm = pwm_from(strs)
        cons = consensus(pwm)
        models.append({
            'source': 'motif_extended' if use_ext else 'motif_core',
            'motif': cons, 'unit_len': len(strs[0]), 'copies': chosen,
            'copy_strings': strs, 'pwm': pwm, 'information_bits': information(pwm),
            'per_unit_support': [round(1 - mm / len(strs[0]), 3) for _, _, mm in chosen],
            'core_copies': len(core), 'extended_copies': len(ext_copies),
            'rejected_cores': rejected, 'core_consensus': core_cons,
            'extended_consensus': ext_cons,
            'interval': None, 'phase': None, 'phase_info_margin': None,
            'boundary_identity_left': None, 'boundary_identity_right': None,
            'score': len(chosen) * len(strs[0]) * information(pwm)})
    models.sort(key=lambda m: -m['score'])
    return models


def classify(model):
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
    p.add_argument('--genomes', type=Path, nargs='+', required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--flank', type=int, default=500)
    p.add_argument('--peaks', type=Path)
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

    models_rows, copy_rows, aln_rows, unresolved, rej_rows = [], [], [], [], []
    for i, c in enumerate(cands, 1):
        contig = c['contig']
        if contig not in seqs:
            unresolved.append({'candidate_id': f"{contig}:{c['start0']}-{c['end0']}",
                               'contig': contig, 'start0': c['start0'],
                               'end0': c['end0'], 'score': c.get('score', ''),
                               'reason': 'no sequence for contig', 'alternatives': '[]'})
            continue
        g = seqs[contig]
        cs, ce = int(c['start0']), int(c['end0'])
        s, e = max(0, cs - a.flank), min(len(g), ce + a.flank)
        sub = ''.join(ch if ch in BASES else 'A' for ch in g[s:e])
        cid = f'{contig}:{cs}-{ce}'
        span = sub[cs - s:ce - s]

        found = motif_models(sub, a.min_copies)
        # periodic branch: independent, no seed required
        per = None
        d, ident = best_period(span)
        if d and ident >= PERIODIC_MIN_IDENTITY:
            iv = periodic_interval(sub, d)
            if iv:
                per = periodic_model(sub, iv[0], iv[1], d)
        alts = []
        lc_frac, lc_unit = low_complexity(span)
        if lc_frac >= 0.6:
            alts.append({'type': 'low_complexity', 'unit': lc_unit,
                         'dominant_fraction': round(lc_frac, 3)})
        if d:
            alts.append({'type': 'periodic_self_match', 'scope': 'candidate_span',
                         'period': d, 'identity_at_period': round(ident, 3)})

        pool = [m for m in [found[0] if found else None, per] if m]
        if not pool:
            unresolved.append({'candidate_id': cid, 'contig': contig, 'start0': cs,
                               'end0': ce, 'score': c.get('score', ''),
                               'reason': 'no motif copies and no periodic interval',
                               'alternatives': json.dumps(alts)})
            continue
        # prefer whichever model explains more sequence
        pool.sort(key=lambda m: -(len(m['copies']) * m['unit_len'] * m['information_bits']))
        m = pool[0]
        kind, s2s, edges = classify(m)
        gstart = s
        inside = sum(1 for aa, _, _ in m['copies'] if cs <= gstart + aa < ce)
        pk = peaks.get(contig, [])
        per_copy_peaks = []
        for idx, (aa, bb, mm) in enumerate(m['copies']):
            g0, g1 = gstart + aa, gstart + bb
            n = sum(1 for ps, pe in pk if ps < g1 and g0 < pe)
            per_copy_peaks.append(n)
            copy_rows.append({'candidate_id': cid, 'contig': contig, 'start0': g0,
                              'end0': g1, 'sequence': m['copy_strings'][idx],
                              'mismatches_to_consensus': mm,
                              'unit_support': m['per_unit_support'][idx],
                              'inside_candidate_span': cs <= g0 < ce,
                              'overlapping_peaks': n})
        for j, cstr in enumerate(m['copy_strings'], 1):
            aln_rows.append({'candidate_id': cid, 'copy_index': j, 'alignment': cstr})
        for r in m['rejected_cores']:
            rej_rows.append({'candidate_id': cid,
                             'core_start0': gstart + r['core_start'], **{
                                 k2: v for k2, v in r.items() if k2 != 'core_start'}})
        iv = m['interval']
        models_rows.append({
            'candidate_id': cid, 'contig': contig, 'cand_start0': cs, 'cand_end0': ce,
            'detector_score': c.get('score', ''),
            'annotation_category': c.get('annotation_category', ''),
            'model_source': m['source'], 'architecture': kind,
            'consensus': m['motif'], 'unit_len': m['unit_len'],
            'n_copies': len(m['copies']), 'n_copies_inside_candidate': inside,
            'core_copies': m['core_copies'],
            'extended_copies': m['extended_copies'] if m['extended_copies'] is not None else '',
            'n_rejected_cores': len(m['rejected_cores']),
            'information_bits_per_position': round(m['information_bits'], 3),
            'median_unit_support': round(float(np.median(m['per_unit_support'])), 3),
            'median_start_to_start': int(np.median(s2s)) if s2s else '',
            'start_to_start_range': f'{min(s2s)}-{max(s2s)}' if s2s else '',
            'median_spacer_edge_gap': int(np.median(edges)) if edges else '',
            'spacer_cv': round(float(np.std(s2s) / np.mean(s2s)), 3) if s2s and np.mean(s2s) else '',
            'periodic_interval_start0': gstart + iv[0] if iv else '',
            'periodic_interval_end0': gstart + iv[1] if iv else '',
            'periodic_phase': m['phase'] if m['phase'] is not None else '',
            'phase_info_margin_bits': m['phase_info_margin'] if m['phase_info_margin'] is not None else '',
            'boundary_identity_left': m['boundary_identity_left'] if m['boundary_identity_left'] is not None else '',
            'boundary_identity_right': m['boundary_identity_right'] if m['boundary_identity_right'] is not None else '',
            'copies_with_>=1_peak': sum(1 for n in per_copy_peaks if n),
            'alternatives': json.dumps(alts),
            'n_alternative_models': len(pool) - 1,
            'indel_handling': 'ungapped_substitutions_only',
            'pwm_file': f'pwm_{i:04d}.npy'})
        np.save(a.out / f'pwm_{i:04d}.npy', m['pwm'])

    for name, rows in [('models.tsv', models_rows), ('copies.tsv', copy_rows),
                       ('alignments.tsv', aln_rows), ('unresolved.tsv', unresolved),
                       ('rejected_cores.tsv', rej_rows)]:
        if rows:
            with (a.out / name).open('w', newline='') as fh:
                w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
                w.writeheader()
                w.writerows(rows)
    arch = Counter(r['architecture'] for r in models_rows)
    src = Counter(r['model_source'] for r in models_rows)
    (a.out / 'summary.json').write_text(json.dumps({
        'version': 'v0.2', 'label': a.label, 'n_candidates': len(cands),
        'n_explained': len(models_rows), 'n_unresolved': len(unresolved),
        'model_sources': dict(src), 'architectures': dict(arch),
        'n_rejected_core_records': len(rej_rows), 'flank_nt': a.flank,
        'periodic_branch': {'period_search_range': [MIN_PERIOD, MAX_PERIOD],
                            'min_identity': PERIODIC_MIN_IDENTITY,
                            'min_units': PERIODIC_MIN_UNITS,
                            'no_seed_required': True,
                            'nothing_hard_coded': 'period, phase, interval and '
                                                  'consensus are inferred per candidate'},
        'search_scope': 'candidate span plus flanks, not restricted to peak intervals',
        'indel_limitation': 'ungapped substitution alignment only; indel-containing '
                            'arrays are under-reported',
        'detector_outputs': 'read-only; not modified'}, indent=2) + '\n')
    print(f'{len(cands)} candidates: {len(models_rows)} explained, '
          f'{len(unresolved)} unresolved')
    for k, v in src.most_common():
        print(f'  source {k:18s} {v}')
    for k, v in arch.most_common():
        print(f'  arch   {k:18s} {v}')


if __name__ == '__main__':
    main()
