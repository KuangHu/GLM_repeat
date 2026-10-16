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
# --- everything from here on completes the truncated upload (186 lines). ---
# --- The upload ended mid-expression in the design.json write just below. ---
    (a.out / 'design.json').write_text(json.dumps(
        {'seeds': SEEDS, 'flank_default': FLANK, 'max_targets': MAX_TARGETS,
         'match_tolerance_nt': MATCH_TOLERANCE, 'loci': {k: {'record_id': v[0], 'flank': v[1],
         'units': [list(u) for u in v[2]], 'note': v[3]} for k, v in loci.items()},
         'mutation': 'every base in a mutated span replaced by a uniformly chosen '
                     'different base',
         'diagnostic_only': 'unit definitions are annotation-derived and are never '
                            'used by the discovery scorer'}, indent=2) + '\n')
    fields = [k for k in rows[0] if k != 'target_positions'] + ['target_positions']
    with (a.out / 'manifest.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({**r, 'target_positions': ','.join(map(str, r['target_positions']))})
    n_t = len({(r['locus'], r['target_index']) for r in rows})
    print(f"{len(rows)} sequences, {n_t} targets, {len(loci)} loci -> {a.out}")
    for locus, (rid, flank, units, note) in loci.items():
        print(f'  {locus} ({rid}): {len(units)} units, targets '
              f'{pick_targets(units)}, flank {flank} - {note}')


# ----------------------------------------------------------------- score
def score(a):
    rows = list(csv.DictReader((a.out / 'manifest.tsv').open(), delimiter='\t'))
    seqs = read_fa(a.out / 'sequences.fa')
    model, tokenizer = rs.load_model(a.snapshot.resolve(), a.device)
    manifest = json.loads((a.snapshot / 'snapshot_manifest.json').read_text())
    out_rows = []
    for i, r in enumerate(rows, 1):
        seq = seqs[r['seq_id']]
        pos = np.array([int(x) for x in r['target_positions'].split(',')])
        # the target must be byte-identical to the original in every arm
        orig = seqs[f"{r['locus']}|t{r['target_index']}|original|0"]
        if seq[pos[0]:pos[-1] + 1] != orig[pos[0]:pos[-1] + 1]:
            raise AssertionError(f"{r['seq_id']}: target positions were mutated")
        rec = {k: r[k] for k in ('locus', 'record_id', 'target_index', 'target_start',
                                 'target_end', 'arm', 'seed', 'n_other_copies',
                                 'mutated_spans', 'distance_mismatch_max')}
        for orient in ('+', '-'):
            sq = seq if orient == '+' else rs.reverse_complement(seq)
            # reverse position i maps to len-1-i; scores are stored in forward coords
            p = pos if orient == '+' else (len(seq) - 1 - pos)[::-1]
            t0 = time.monotonic()
            v = score_target(model, tokenizer, sq, p, device=a.device,
                             batch_size=a.batch_size,
                             attention_mask_policy=a.attention_mask_policy)
            if orient == '-':
                v = v[::-1]
            out_rows.append({**rec, 'orientation': orient, 'n_positions': len(pos),
                             'target_mean': float(v.mean()), 'target_min': float(v.min()),
                             'target_max': float(v.max()),
                             'seconds': time.monotonic() - t0,
                             'scores': ','.join(f'{x:.4f}' for x in v)})
        print(f"[{i}/{len(rows)}] {r['seq_id']} "
              f"{out_rows[-2]['target_mean']:+.3f} / {out_rows[-1]['target_mean']:+.3f}",
              flush=True)
    fields = [k for k in out_rows[0] if k != 'scores'] + ['scores']
    with (a.out / 'scores.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t')
        w.writeheader()
        w.writerows(out_rows)
    (a.out / 'score_run.json').write_text(json.dumps({
        'model': manifest['model_id'], 'revision': manifest['revision'],
        'dtype': str(next(model.parameters()).dtype),
        'attention_mask_policy': a.attention_mask_policy,
        'mode': 'exact single-position masking, target positions only',
        'both_orientations': True, 'n_rows': len(out_rows),
        'code_sha256': {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                        for f in (Path(__file__).resolve(), Path(rs.__file__).resolve())},
    }, indent=2) + '\n')
    print(a.out / 'scores.tsv')


def score_target(model, tokenizer, sequence, positions, *, device, batch_size,
                 attention_mask_policy):
    """Exact single-nucleotide masked score at the given positions only."""
    import torch
    original, attention, tok_pos, base_ids = rs.tokenize_dna(tokenizer, sequence)
    truth = np.array([rs.DNA.index(sequence[i]) for i in positions])
    out = np.zeros(len(positions), np.float64)
    with torch.inference_mode():
        for off in range(0, len(positions), batch_size):
            chunk = positions[off:off + batch_size]
            ids = original.repeat(len(chunk), 1).to(device)
            attn = rs.resolve_attention_mask(
                attention.repeat(len(chunk), 1).to(device), attention_mask_policy)
            for row, p in enumerate(chunk):
                ids[row, tok_pos[p]] = tokenizer.mask_token_id
            logits = model(input_ids=ids, attention_mask=attn, return_dict=True).logits
            sel = logits[torch.arange(len(chunk), device=device),
                         torch.as_tensor(tok_pos[chunk], device=device)][:, base_ids].float()
            lp = torch.log_softmax(sel, -1).gather(
                1, torch.as_tensor(truth[off:off + len(chunk)], device=device)[:, None])
            v = lp.squeeze(1).cpu().numpy()
            if not np.isfinite(v).all():
                raise FloatingPointError('Non-finite masked probability.')
            out[off:off + len(chunk)] = v
    return out + np.log(4)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build', help='write perturbed sequences + manifest (CPU)')
    b.add_argument('--bench', type=Path, required=True)
    b.add_argument('--out', type=Path, required=True)
    s = sub.add_parser('score', help='exact scoring of target positions (GPU)')
    s.add_argument('--out', type=Path, required=True, help='directory written by build')
    s.add_argument('--snapshot', type=Path, required=True)
    s.add_argument('--device', default='cuda:0')
    s.add_argument('--batch-size', type=int, default=8)
    s.add_argument('--attention-mask-policy', choices=rs.ATTENTION_MASK_POLICIES,
                   default=rs.DEFAULT_ATTENTION_MASK_POLICY)
    s.add_argument('--trust-reviewed-code', action='store_true')
    a = p.parse_args()
    if a.cmd == 'build':
        if a.out.exists() and any(a.out.iterdir()):
            p.error('Output directory is nonempty; use a new directory.')
        build(a)
    else:
        if not a.trust_reviewed_code:
            p.error('Review the downloaded model code, then pass --trust-reviewed-code.')
        if (a.out / 'scores.tsv').exists():
            p.error('scores.tsv already exists; use a new directory.')
        score(a)


if __name__ == '__main__':
    main()