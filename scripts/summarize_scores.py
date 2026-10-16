#!/usr/bin/env python3
"""Summarize scored windows: timing, per-window block features, label vs paired background.

Features follow WORKFLOW_SETUP.md section 7 (raw quantiles, block maxima with
quarter-block steps, top-20% block statistic). Window score = max over the two
scoring orientations. These are development diagnostics, not calibrated calls.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

SCALES = [64, 128, 256, 512, 1024, 2048]
TOP_FRACTION = 0.20


def running_mean(v, k=16):
    if len(v) < k:
        return np.array([v.mean()])
    c = np.cumsum(np.insert(v.astype(np.float64), 0, 0))
    return (c[k:] - c[:-k]) / k


def block_features(v):
    f = {'mean': float(v.mean()), 'sd': float(v.std()),
         'q50': float(np.quantile(v, .5)), 'q90': float(np.quantile(v, .9)),
         'max16': float(running_mean(v).max())}
    for size in SCALES:
        size_eff = min(size, len(v))
        step = max(1, size_eff // 4)
        starts = list(range(0, len(v) - size_eff + 1, step))
        if starts[-1] != len(v) - size_eff:
            starts.append(len(v) - size_eff)
        k = max(1, int(round(TOP_FRACTION * size_eff)))
        means, tops = [], []
        for s in starts:
            b = v[s:s + size_eff]
            means.append(b.mean())
            tops.append(np.partition(b, -k)[-k:].mean())
        f[f'blk{size}_mean_max'] = float(max(means))
        f[f'blk{size}_top20_max'] = float(max(tops))
        f[f'blk{size}_n'] = len(starts)
    return f


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, required=True, help='dir with run.json + NPZ tracks')
    p.add_argument('--manifest', type=Path, required=True, help='window manifest TSV')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    run = json.loads((a.scores / 'run.json').read_text())
    man = {r['record_id']: r for r in csv.DictReader(a.manifest.open(), delimiter='\t')}

    rows, tracks = [], {}
    for rec in run['records']:
        v = np.load(a.scores / rec['track'])['above_uniform'].astype(np.float64)
        m = man[rec['id']]
        tracks[(rec['id'], rec['scoring_orientation'])] = v
        row = {'record_id': rec['id'], 'family': m['family'] or 'background',
               'paired_positive_id': m.get('paired_positive_id', ''),
               'orientation': rec['scoring_orientation'], 'length': rec['length'],
               'gc': float(m['gc_fraction']), 'cds_fraction': float(m['cds_fraction']),
               'inference_s': rec.get('inference_seconds', rec.get('seconds')),
               'variants_per_s': rec.get('variants_per_second'),
               'peak_gpu_mb': rec.get('peak_gpu_mb')}
        if m['label_start_in_window0']:
            ls, le = int(m['label_start_in_window0']), int(m['label_end_in_window0'])
            inside = np.zeros(len(v), bool)
            inside[ls:le] = True
            row.update(label_mean=float(v[inside].mean()),
                       outside_mean=float(v[~inside].mean()),
                       label_max16=float(running_mean(v[ls:le]).max()))
        row.update(block_features(v))
        rows.append(row)

    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (a.out / 'window_orientation_features.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # window score = max over orientations, per feature
    feats = [k for k in rows[0] if k == 'mean' or k.startswith('blk') and not k.endswith('_n')
             or k in ('q90', 'max16')]
    by_rec = {}
    for r in rows:
        by_rec.setdefault(r['record_id'], []).append(r)
    win = {rid: {f: max(x[f] for x in rs) for f in feats} for rid, rs in by_rec.items()}

    lines = ['# Exact-score diagnostics (development only, not calibrated)', '',
             f"model {run['model']}@{run['revision'][:10]}  mode {run['mode']}  "
             f"batch {run.get('batch_size', 1)}  dtype {run.get('dtype', '?')}  "
             f"device {run.get('device_name', '?')}",
             f"model load {run.get('model_load_seconds', float('nan')):.1f}s", '',
             '## Timing per orientation', '',
             '| record | family | orient | variants | inference s | variants/s | peak GPU MB |',
             '|---|---|---|---:|---:|---:|---:|']
    for r, rec in zip(rows, run['records']):
        lines.append(f"| {r['record_id']} | {r['family']} | {r['orientation']} | "
                     f"{rec['n_variants']} | {r['inference_s']:.1f} | "
                     f"{(r['variants_per_s'] or float('nan')):.1f} | "
                     f"{(r['peak_gpu_mb'] or float('nan')):.0f} |")
    tot = sum(r['inference_s'] for r in rows)
    nt = sum(r['length'] for r in rows)
    lines += ['', f'Total inference {tot:.0f}s for {nt} scored nt-orientations '
              f'({nt / tot:.1f} nt/s).', '',
              '## Label region vs rest of window (above_uniform, nats)', '',
              '| record | family | orient | label mean | outside mean | label max16 |',
              '|---|---|---|---:|---:|---:|']
    for r in rows:
        if 'label_mean' in r:
            lines.append(f"| {r['record_id']} | {r['family']} | {r['orientation']} | "
                         f"{r['label_mean']:.3f} | {r['outside_mean']:.3f} | "
                         f"{r['label_max16']:.3f} |")
    show = ['mean', 'max16', 'blk64_mean_max', 'blk256_mean_max', 'blk1024_mean_max',
            'blk64_top20_max', 'blk256_top20_max', 'blk1024_top20_max']
    lines += ['', '## Positive vs paired background (window score = max over orientations)', '',
              '| positive | background | ' + ' | '.join(show) + ' |',
              '|---|---|' + '---:|' * len(show)]
    for rid, rs in by_rec.items():
        if rs[0]['family'] == 'background':
            continue
        for bid, brs in by_rec.items():
            if brs[0]['paired_positive_id'] != rid:
                continue
            cells = [f"{win[rid][f]:.3f} vs {win[bid][f]:.3f}" for f in show]
            lines.append(f"| {rid} ({rs[0]['family']}) | {bid} | " + ' | '.join(cells) + ' |')
    lines += ['', 'n=3 pairs: descriptive only. Backgrounds are unlabeled real DNA, '
              'not verified negatives.']
    (a.out / 'report.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
