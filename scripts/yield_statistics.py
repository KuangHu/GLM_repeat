#!/usr/bin/env python3
"""Authoritative candidate-yield statistics for every scan round.

Earlier rounds reported a single panel figure, total candidates divided by total
scanned Mb, and compared those figures across rounds. That was wrong twice over:

  * Candidates cluster inside genomes. The sampling unit is the GENOME, not the
    window and not the nucleotide, so any interval must resample genomes.
  * The pooled rate is a ratio estimator weighted by genome length, so a single
    long or repeat-rich genome moves it a long way. With 3-4 genomes it is not
    stable enough to compare across rounds, which is exactly the mistake that
    produced the retracted "host ~37/Mb" and "phage is 2.7x host" claims.

This module therefore reports, for every run:
  pooled_per_Mb        total candidates / total unique scanned Mb (ratio estimator)
  pooled_ci95          genome-level bootstrap CI on that ratio
  per_genome_*         the distribution that must accompany any headline number
  n_genomes            because the CI width is driven by this, not by Mb scanned

It computes nothing about the detector and changes no frozen output.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import statistics as st
from collections import defaultdict
from pathlib import Path

MIN_GENOME_NT = 50_000       # below this a per-genome rate is too noisy to plot


def pctl(v, p):
    if not v:
        return float('nan')
    k = (len(v) - 1) * p
    f = int(k)
    c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)


def unique_scanned_per_contig(run: Path):
    """Union of window spans per contig; never the summed window length."""
    spans = defaultdict(list)
    with (run / 'prepared/windows.tsv').open() as fh:
        for w in csv.DictReader(fh, delimiter='\t'):
            spans[w['contig']].append((int(w['start0']), int(w['end0'])))
    out = {}
    for c, sp in spans.items():
        sp.sort()
        tot, cur_end = 0, -1
        for s, e in sp:
            if e > cur_end:
                tot += e - max(s, cur_end)
                cur_end = max(cur_end, e)
        out[c] = tot
    return out


def load(run: Path):
    uniq = unique_scanned_per_contig(run)
    cnt = defaultdict(int)
    with (run / 'calls/candidates.tsv').open() as fh:
        for c in csv.DictReader(fh, delimiter='\t'):
            cnt[c['contig']] += 1
    trna = defaultdict(int)
    tp = run / 'calls/trna_overlap.tsv'
    if tp.exists():
        with tp.open() as fh:
            for c in csv.DictReader(fh, delimiter='\t'):
                if int(c['n_tRNA_overlapping'] or 0):
                    trna[c['contig']] += 1
    return [{'contig': c, 'nt': uniq[c], 'cands': cnt[c], 'trna': trna[c]}
            for c in sorted(uniq)]


MIN_GENOMES_FOR_CI = 5       # below this a genome bootstrap is not meaningful


def bootstrap_ratio(genomes, n_boot, rng):
    """Resample GENOMES with replacement; the ratio is recomputed each draw.

    With very few genomes the bootstrap has almost no distinct resamples to draw
    on (n=3 gives 10 distinct multisets), so the interval it returns is an
    artefact of the resampling, not a measure of uncertainty. Such intervals are
    withheld rather than printed, because a narrow-looking interval from n=3 is
    more misleading than no interval at all.
    """
    if len(genomes) < MIN_GENOMES_FOR_CI:
        return (float('nan'), float('nan'))
    out = []
    k = len(genomes)
    for _ in range(n_boot):
        draw = [genomes[rng.randrange(k)] for _ in range(k)]
        nt = sum(g['nt'] for g in draw)
        if nt:
            out.append(sum(g['cands'] for g in draw) / (nt / 1e6))
    out.sort()
    return (round(pctl(out, 0.025), 1), round(pctl(out, 0.975), 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='append', required=True, metavar='LABEL=DIR')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--bootstrap', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=20261017)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    a.out.mkdir(parents=True, exist_ok=True)

    rows, per_genome_rows = [], []
    for spec in a.run:
        lab, _, d = spec.partition('=')
        genomes = load(Path(d))
        big = [g for g in genomes if g['nt'] >= MIN_GENOME_NT]
        dens = sorted(g['cands'] / (g['nt'] / 1e6) for g in big)
        tot_c = sum(g['cands'] for g in genomes)
        tot_nt = sum(g['nt'] for g in genomes)
        lo, hi = bootstrap_ratio(genomes, a.bootstrap, rng)
        n_tr = sum(g['trna'] for g in genomes)
        ci_ok = len(genomes) >= MIN_GENOMES_FOR_CI
        rows.append({
            'run': lab, 'n_genomes': len(genomes),
            'n_genomes_rated_for_distribution': len(big),
            'ci_reportable': ci_ok,
            'ci_withheld_reason': '' if ci_ok else
                f'only {len(genomes)} genomes; a genome bootstrap below '
                f'{MIN_GENOMES_FOR_CI} genomes has too few distinct resamples to '
                f'mean anything, so no interval is reported',
            'unique_scanned_Mb': round(tot_nt / 1e6, 3),
            'n_candidates': tot_c,
            'pooled_per_Mb': round(tot_c / (tot_nt / 1e6), 1),
            'pooled_ci95_lo': lo, 'pooled_ci95_hi': hi,
            'ci95_width': round(hi - lo, 1) if hi == hi else float('nan'),
            'per_genome_median': round(st.median(dens), 1) if dens else '',
            'per_genome_p25': round(pctl(dens, 0.25), 1) if dens else '',
            'per_genome_p75': round(pctl(dens, 0.75), 1) if dens else '',
            'per_genome_min': round(min(dens), 1) if dens else '',
            'per_genome_max': round(max(dens), 1) if dens else '',
            'per_genome_fold_range': round(max(dens) / min(dens), 1) if dens and min(dens) else '',
            'pct_candidates_tRNA_overlapping': round(100 * n_tr / tot_c, 1) if tot_c else '',
        })
        for g in big:
            per_genome_rows.append({
                'run': lab, 'contig': g['contig'],
                'unique_scanned_nt': g['nt'], 'n_candidates': g['cands'],
                'per_Mb': round(g['cands'] / (g['nt'] / 1e6), 1),
                'n_tRNA_overlapping_candidates': g['trna'],
                'pct_tRNA': round(100 * g['trna'] / g['cands'], 1) if g['cands'] else '',
            })

    with (a.out / 'yield_by_run.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    per_genome_rows.sort(key=lambda r: (r['run'], -r['per_Mb']))
    with (a.out / 'yield_by_genome.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(per_genome_rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(per_genome_rows)

    (a.out / 'method.json').write_text(json.dumps({
        'sampling_unit': 'genome',
        'pooled_per_Mb': 'total candidates / total unique scanned Mb; a ratio '
                         'estimator weighted by genome length',
        'unique_scanned_Mb': 'union of window spans per contig, NOT summed window '
                             'length, which would double-count the 50% stride overlap',
        'interval': f'{a.bootstrap} genome-level bootstrap resamples, seed {a.seed}; '
                    'genomes are resampled with replacement and the ratio recomputed, '
                    'because candidates cluster within genomes',
        'why_not_window_level': 'window- or nucleotide-level intervals would assume '
                                'candidates are independent across a genome, which they '
                                'are not, and would understate the interval several-fold',
        'ci_withheld_below_n': MIN_GENOMES_FOR_CI,
        'min_genome_nt_for_per_genome_rate': MIN_GENOME_NT,
        'n_genomes_vs_rated': 'n_genomes counts every scanned contig; the per-genome '
                              f'distribution uses only those with >={MIN_GENOME_NT} nt '
                              'scanned, since a rate over a few kb is dominated by noise',
        'reporting_rule': 'a pooled rate must never be quoted without its per-genome '
                          'distribution and n_genomes. Panels below ~20 genomes have '
                          'intervals too wide to support cross-round comparison.',
        'not_an_FDR': 'these are candidate yields, not FDR, precision or biological '
                      'probabilities; no candidate here is human-confirmed',
    }, indent=2) + '\n')

    hdr = (f"{'run':10s} {'n_gen':>5s} {'rated':>5s} {'Mb':>6s} {'cands':>6s} "
           f"{'pooled/Mb':>10s} {'95% CI (genome bootstrap)':>26s} {'median':>7s} "
           f"{'range':>13s} {'tRNA%':>6s}")
    print(hdr)
    print('-' * len(hdr))
    for r in rows:
        ci = (f"{r['pooled_ci95_lo']}-{r['pooled_ci95_hi']}" if r['ci_reportable']
              else f'withheld (n={r["n_genomes"]})')
        rg = f"{r['per_genome_min']}-{r['per_genome_max']}"
        print(f"{r['run']:10s} {r['n_genomes']:5d} {r['n_genomes_rated_for_distribution']:5d} "
              f"{r['unique_scanned_Mb']:6.2f} {r['n_candidates']:6d} "
              f"{r['pooled_per_Mb']:10.1f} {ci:>26s} {r['per_genome_median']:7.1f} "
              f"{rg:>13s} {r['pct_candidates_tRNA_overlapping']:5.1f}%")
    print(f'\n-> {a.out}')


if __name__ == '__main__':
    main()
