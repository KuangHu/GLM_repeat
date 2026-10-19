#!/usr/bin/env python3
"""Near-duplicate detection for a candidate panel, against prior panels and itself.

This exists because an earlier inline version passed minimap2 both `-X` and
`--secondary=no`, which that tool refuses:

    [ERROR] -X/-P and --secondary=no can't be applied at the same time

stderr was discarded, the PAF came back empty, and the pipeline reported "0
internal near-duplicate pairs" for round-3 and db10 when the comparison had never
run. Two real duplicate pairs were scanned as a result.

So this script refuses to report a negative it has not earned:
  * minimap2's exit status is checked and stderr is kept, never discarded
  * an empty alignment file from a non-trivial input is an ERROR, not a zero
  * a positive control is run: the panel is compared to itself WITHOUT -X, which
    must find every sequence matching itself, or the comparison is not trusted

Criteria: identity >= 95% and both coverages >= 85%, aggregated over all
alignments between a pair.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def seq_lengths(fa):
    L, name, n = {}, None, 0
    with Path(fa).open() as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    L[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name:
        L[name] = n
    return L


def run_minimap(mm, ref, qry, out, err, extra):
    cmd = [mm, '-x', 'asm20', '-t', '8', *extra, str(ref), str(qry)]
    with Path(out).open('w') as fo, Path(err).open('w') as fe:
        rc = subprocess.call(cmd, stdout=fo, stderr=fe)
    msg = Path(err).read_text()
    if rc != 0:
        raise RuntimeError(f'minimap2 failed (exit {rc}) for {" ".join(cmd)}\n{msg}')
    if '[ERROR]' in msg:
        raise RuntimeError(f'minimap2 reported an error for {" ".join(cmd)}\n{msg}')
    return Path(out).stat().st_size


def pairs_from(paf, La, Lb, same):
    agg = defaultdict(lambda: defaultdict(int))
    with Path(paf).open() as fh:
        for line in fh:
            f = line.split('\t')
            q, t = f[0], f[5]
            if same and q == t:
                continue
            if int(f[10]) == 0:
                continue
            key = (q, t) if not same or q < t else (t, q)
            d = agg[key]
            d['m'] += int(f[9])
            d['a'] += int(f[10])
            if key[0] == q:
                d['qc'] += int(f[3]) - int(f[2])
                d['tc'] += int(f[8]) - int(f[7])
            else:
                d['tc'] += int(f[3]) - int(f[2])
                d['qc'] += int(f[8]) - int(f[7])
    out = []
    for (a, b), d in agg.items():
        ident = d['m'] / d['a']
        ca = min(1.0, d['qc'] / La.get(a, Lb.get(a, 1)))
        cb = min(1.0, d['tc'] / Lb.get(b, La.get(b, 1)))
        if ident >= 0.95 and ca >= 0.85 and cb >= 0.85:
            out.append((a, b, round(ident, 4), round(ca, 3), round(cb, 3)))
    return sorted(out), len(agg)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--panel', type=Path, required=True)
    ap.add_argument('--prior', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--minimap2', default='minimap2')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    Lp = seq_lengths(a.panel)
    print(f'panel: {len(Lp):,} sequences, {sum(Lp.values())/1e6:.2f} Mb')

    # ---- positive control: self-vs-self WITHOUT -X must find self matches
    size = run_minimap(a.minimap2, a.panel, a.panel, a.out / 'control.paf',
                       a.out / 'control.err', [])
    if size == 0:
        raise RuntimeError('positive control produced no alignments at all; the '
                           'comparison is not working, so no negative result from '
                           'this run can be trusted')
    selfhits = sum(1 for line in (a.out / 'control.paf').open()
                   if line.split('\t')[0] == line.split('\t')[5])
    if selfhits == 0:
        raise RuntimeError('positive control found no sequence matching itself')
    print(f'positive control: {selfhits:,} self-matches -> comparison works')

    drop, report = set(), []
    if a.prior:
        Lr = seq_lengths(a.prior)
        run_minimap(a.minimap2, a.prior, a.panel, a.out / 'vs_prior.paf',
                    a.out / 'vs_prior.err', ['--secondary=no'])
        hits, naln = pairs_from(a.out / 'vs_prior.paf', Lp, Lr, same=False)
        print(f'vs prior : {naln:,} aligned pairs -> {len(hits)} near-duplicates')
        for h in hits:
            drop.add(h[0])
            report.append(('vs_prior', *h))

    # -X skips self and dual mappings; it CANNOT be combined with --secondary=no
    run_minimap(a.minimap2, a.panel, a.panel, a.out / 'self.paf',
                a.out / 'self.err', ['-X'])
    hits, naln = pairs_from(a.out / 'self.paf', Lp, Lp, same=True)
    print(f'internal : {naln:,} aligned pairs -> {len(hits)} near-duplicates')
    for a_, b_, i, ca, cb in hits:
        keep, lose = (a_, b_) if Lp[a_] >= Lp[b_] else (b_, a_)
        if keep not in drop:
            drop.add(lose)
        report.append(('internal', a_, b_, i, ca, cb))

    with (a.out / 'dedup_pairs.tsv').open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['comparison', 'seq_a', 'seq_b', 'identity', 'coverage_a', 'coverage_b'])
        w.writerows(report)
    (a.out / 'drop.txt').write_text('\n'.join(sorted(drop)) + ('\n' if drop else ''))
    print(f'\n{len(drop)} sequences to drop -> {a.out}/drop.txt')


if __name__ == '__main__':
    main()
