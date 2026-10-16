#!/usr/bin/env python3
"""Compare threshold policies over one window panel.

Takes several peak_cluster output directories (label=dir) and tabulates window
scores, non-zero rates, ranking and unreachable-threshold tracks. Any window
singled out with --highlight is reported in its own section.
"""
import argparse
import csv
import json
from pathlib import Path


def load(d):
    run = json.loads((d / 'run.json').read_text())
    ws = {r['record_id']: float(r['peak_cluster_score'])
          for r in csv.DictReader((d / 'window_scores.tsv').open(), delimiter='\t')}
    cands = (list(csv.DictReader((d / 'candidates_merged.tsv').open(), delimiter='\t'))
             if (d / 'candidates_merged.tsv').exists() else [])
    return run, ws, cands


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs', nargs='+', required=True, help='label=path pairs')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--highlight', nargs='*', default=[])
    p.add_argument('--title', default='Threshold policy comparison')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = {r['record_id']: r for r in
           csv.DictReader(a.manifest.open(), delimiter='\t')}
    runs = {}
    for spec in a.runs:
        label, path = spec.split('=', 1)
        runs[label] = load(Path(path))
    ids = sorted(next(iter(runs.values()))[1])
    pos = [i for i in ids if man[i]['family'] not in ('background', 'calibration_background')]
    bg = [i for i in ids if i not in pos]

    L = [f'# {a.title}', '',
         'Thresholds were calibrated on background only and frozen before evaluation '
         '(results/calibration/thresholds/FROZEN.md). The original `window_local` '
         'scorer and its report are preserved unchanged.', '',
         '## Policy summary', '',
         '| policy | threshold | merged candidates | positives non-zero | '
         'backgrounds non-zero | tracks with unreachable threshold |',
         '|---|---:|---:|---:|---:|---:|']
    for label, (run, ws, cands) in runs.items():
        unr = run.get('unreachable_threshold_tracks', [])
        tau = run.get('threshold')
        L.append(f"| {label} | {'per-window' if tau is None else tau} | {len(cands)} | "
                 f"{sum(ws[i] > 0 for i in pos)}/{len(pos)} | "
                 f"{sum(ws[i] > 0 for i in bg)}/{len(bg)} | {len(unr)} |")
    unreach = {lab: run.get('unreachable_threshold_tracks', [])
               for lab, (run, _, _) in runs.items()}
    if any(unreach.values()):
        L += ['', 'Unreachable thresholds (tau >= ln(4) = 1.3863). Such a window yields '
              'zero peaks **by construction**; this is never evidence that no '
              'repetitive organization exists:', '']
        for lab, u in unreach.items():
            for t in u:
                L.append(f"- {lab}: {t['window_id']} ({t['orientation']}) tau = "
                         f"{t['threshold']}")
    L += ['', '## Window scores (all windows)', '',
          '| record | family | ' + ' | '.join(runs) + ' |',
          '|---|---|' + '---:|' * len(runs)]
    order = sorted(ids, key=lambda i: -max(runs[l][1][i] for l in runs))
    for i in order:
        mark = ' **<-**' if i in a.highlight else ''
        L.append(f"| {i}{mark} | {man[i]['family']} | " +
                 ' | '.join(f'{runs[l][1][i]:.3f}' for l in runs) + ' |')
    L += ['', '## Rank of each positive (of %d windows)' % len(ids), '',
          '| record | family | ' + ' | '.join(runs) + ' |',
          '|---|---|' + '---:|' * len(runs)]
    ranks = {l: {r: k + 1 for k, r in
                 enumerate(sorted(ids, key=lambda x: -runs[l][1][x]))} for l in runs}
    for i in pos:
        L.append(f"| {i} | {man[i]['family']} | " +
                 ' | '.join(str(ranks[l][i]) for l in runs) + ' |')
    for h in a.highlight:
        if h not in man:
            continue
        L += ['', f'## {h} reported separately', '',
              f"Panel label: **{man[h]['family']}** (unchanged). "
              f"Source {man[h]['source_id']}.", '',
              '| policy | score | rank | non-zero |', '|---|---:|---:|---|']
        for l in runs:
            s = runs[l][1][h]
            L.append(f'| {l} | {s:.3f} | {ranks[l][h]} | {"yes" if s > 0 else "no"} |')
    L += ['', 'Backgrounds are unlabeled real DNA, not confirmed negatives. '
          'No threshold was tuned against positives and none was re-tuned after '
          'this evaluation.']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
