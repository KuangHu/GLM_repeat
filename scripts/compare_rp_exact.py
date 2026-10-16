#!/usr/bin/env python3
"""Random-partition (approximate) vs exact masked scoring on the seed panel.

Evaluates (1) per-base agreement and bias, (2) local signal preservation:
16-nt smoothed tracks, top-5% position overlap, label and annotated sub-region
means, best-block location, (3) label-blind window ranking with the fixed
RANK_FEATURES, (4) seed-to-seed variability, (5) ART loci separately.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_tracks import load_runs  # noqa: E402
from report_seed_panel import RANK_FEATURES, auc  # noqa: E402
from summarize_scores import block_features, running_mean  # noqa: E402


def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def top_jaccard(a, b, frac=0.05):
    k = max(1, int(len(a) * frac))
    ta, tb = set(np.argsort(a)[-k:]), set(np.argsort(b)[-k:])
    return len(ta & tb) / len(ta | tb)


def best_block(v, size=256, step=64):
    return max(range(0, len(v) - size + 1, step), key=lambda i: v[i:i + size].mean())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exact', type=Path, nargs='+', required=True)
    p.add_argument('--rp', type=Path, nargs='+', required=True,
                   help='one run dir per seed; the first is the primary seed')
    p.add_argument('--bench', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = {r['record_id']: r for r in
           csv.DictReader((a.bench / 'data/seed19_manifest.tsv').open(), delimiter='\t')}
    subs = list(csv.DictReader((a.bench / 'data/annotated_repeat_subregions.tsv').open(),
                               delimiter='\t'))
    ex = load_runs(a.exact)
    rps = [load_runs([d]) for d in a.rp]
    rp = rps[0]
    keys = sorted(set(ex) & set(rp))
    ids = sorted({k[0] for k in keys})
    L = ['# Random partition (approximate) vs exact scoring: seed panel', '',
         f'exact: {", ".join(map(str, a.exact))}',
         f'random partition: {", ".join(map(str, a.rp))} (first = primary seed)',
         f'tracks compared: {len(keys)} ({len(ids)} windows)', '']
    rmeta = rp[keys[0]][1]
    L.append(f"rp settings: groups {rmeta['groups']}, repeats {rmeta['repeats']}, "
             f"variants/orientation {rp[keys[0]][2]['n_variants']} vs exact "
             f"{ex[keys[0]][2]['n_variants']}; policy {rmeta.get('attention_mask_policy')}")
    tsec = sum(r[2]['inference_seconds'] for r in rp.values())
    esec = sum(ex[k][2]['inference_seconds'] for k in keys)
    L += [f'measured inference time, these tracks: rp {tsec:.0f} s vs exact {esec:.0f} s '
          f'({esec / tsec:.0f}x)', '']

    # ---- per-base and local agreement per track
    rows = []
    for k in keys:
        e, r = ex[k][0], rp[k][0]
        es, rs_ = running_mean(e), running_mean(r)
        rows.append({'record_id': k[0], 'family': man[k[0]]['family'], 'orient': k[1],
                     'bias': float((r - e).mean()), 'mae': float(np.abs(r - e).mean()),
                     'p99': float(np.quantile(np.abs(r - e), .99)),
                     'pearson': float(np.corrcoef(e, r)[0, 1]),
                     'pearson16': float(np.corrcoef(es, rs_)[0, 1]),
                     'top5_jaccard16': top_jaccard(es, rs_),
                     'best256_exact': best_block(e), 'best256_rp': best_block(r)})
    with (a.out / 'per_track.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    agg = lambda f, rs: f"{np.median([x[f] for x in rs]):.3f} [{min(x[f] for x in rs):.3f}, {max(x[f] for x in rs):.3f}]"
    L += ['## Per-base and local agreement (median [min, max] over tracks)', '',
          '| group | n tracks | bias rp-exact | MAE | p99 abs | Pearson | Pearson 16-nt | top-5% Jaccard (16-nt) | best 256-nt block same |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for g in ['all', 'CRISPR', 'VIPR', 'ART', 'background']:
        rs = rows if g == 'all' else [x for x in rows if x['family'] == g]
        same = sum(x['best256_exact'] == x['best256_rp'] for x in rs)
        L.append(f"| {g} | {len(rs)} | {agg('bias', rs)} | {agg('mae', rs)} | {agg('p99', rs)} | "
                 f"{agg('pearson', rs)} | {agg('pearson16', rs)} | {agg('top5_jaccard16', rs)} | "
                 f"{same}/{len(rs)} |")

    # ---- seed variability
    if len(rps) > 1:
        sd_rows = []
        for k in keys:
            stack = np.stack([R[k][0] for R in rps if k in R])
            sd_rows.append((np.abs(stack[1:] - stack[0]).mean(), np.abs(rp[k][0] - ex[k][0]).mean()))
        L += ['', f'Seed variability ({len(rps)} seeds): mean |seed_i - primary| per base '
              f'{np.mean([s for s, _ in sd_rows]):.4f} vs mean |rp - exact| '
              f'{np.mean([d for _, d in sd_rows]):.4f}.']

    # ---- label and sub-region preservation
    L += ['', '## Label and annotated sub-region means (orientation-averaged)', '',
          '| record | family | region | exact | rp (seed mean +/- sd) | rp - exact |',
          '|---|---|---|---:|---:|---:|']
    for rid in ids:
        m = man[rid]
        if not m['label_start_in_window0']:
            continue
        w0 = int(m['start0'])
        ls, le = int(m['label_start_in_window0']), int(m['label_end_in_window0'])
        regions = [('label', [(ls, le)])]
        mine = [s for s in subs if s['accession_version'] == m['accession_version']
                and w0 + ls <= int(s['start0']) and int(s['end0']) <= w0 + le]
        if mine:
            regions.append((mine[0]['feature'].rsplit('_', 1)[0] if len(mine) > 1 else mine[0]['feature'],
                            [(int(s['start0']) - w0, int(s['end0']) - w0) for s in mine]))
        avg = lambda R: (R[(rid, '+')][0] + R[(rid, '-')][0]) / 2
        e = avg(ex)
        for name, spans in regions:
            msk = np.zeros(len(e), bool)
            for s, t in spans:
                msk[s:t] = True
            rv = [avg(R)[msk].mean() for R in rps if (rid, '+') in R]
            L.append(f"| {rid} | {m['family']} | {name} | {e[msk].mean():.3f} | "
                     f"{np.mean(rv):.3f} +/- {np.std(rv):.3f} | {np.mean(rv) - e[msk].mean():+.3f} |")

    # ---- ranking
    def windows(R):
        out = {}
        for rid in ids:
            fo = [block_features(R[(rid, o)][0]) for o in '+-']
            out[rid] = {f: max(x[f] for x in fo) for f in RANK_FEATURES}
        return out
    we, wr = windows(ex), windows(rp)
    pos = [r for r in ids if man[r]['family'] != 'background']
    bg = [r for r in ids if man[r]['family'] == 'background']
    L += ['', '## Label-blind window ranking (window score = max over orientations)', '',
          '| feature | Spearman exact vs rp (19 windows) | positive ranks exact | positive ranks rp | P(pos>bg) exact | P(pos>bg) rp |',
          '|---|---:|---|---|---:|---:|']
    for f in RANK_FEATURES:
        rank = lambda W: {r: i + 1 for i, r in enumerate(sorted(ids, key=lambda x: -W[x][f]))}
        rke, rkr = rank(we), rank(wr)
        L.append(f"| {f} | {spearman([we[r][f] for r in ids], [wr[r][f] for r in ids]):.3f} | "
                 + ', '.join(f"{man[r]['family']}:{rke[r]}" for r in pos) + ' | '
                 + ', '.join(f"{man[r]['family']}:{rkr[r]}" for r in pos) + ' | '
                 f"{auc([we[r][f] for r in pos], [we[r][f] for r in bg]):.2f} | "
                 f"{auc([wr[r][f] for r in pos], [wr[r][f] for r in bg]):.2f} |")

    # ---- ART separately
    L += ['', '## ART loci in detail', '']
    for rid in [r for r in pos if man[r]['family'] == 'ART']:
        m = man[rid]
        w0 = int(m['start0'])
        cores = [s for s in subs if s['accession_version'] == m['accession_version']
                 and 'ART_core' in s['feature']]
        avg = lambda R: (R[(rid, '+')][0] + R[(rid, '-')][0]) / 2
        e, r = avg(ex), avg(rp)
        L.append(f"- {rid} ({m['locus_id']}): " + '; '.join(
            f"core {c['feature'][-1]} exact {e[int(c['start0']) - w0:int(c['end0']) - w0].mean():.3f} "
            f"rp {r[int(c['start0']) - w0:int(c['end0']) - w0].mean():.3f}" for c in cores))
        for f in ['blk256_top20_max', 'blk1024_top20_max', 'blk1024_mean_max']:
            L.append(f"  - {f}: exact {we[rid][f]:.3f} vs rp {wr[rid][f]:.3f}; paired bg "
                     + ', '.join(f"{b} {we[b][f]:.3f}/{wr[b][f]:.3f}" for b in bg
                                 if man[b]['paired_positive_id'] == rid) + ' (exact/rp)')
    L += ['', 'Descriptive only (5 positives). The approximation needs its own calibration '
          'and the recall-loss gate (WORKFLOW_SETUP.md section 9) on a larger panel.']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
