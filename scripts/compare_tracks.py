#!/usr/bin/env python3
"""Compare candidate score tracks against an immutable reference run.

Reports per-base |diff| quantiles and maxima, block-feature differences, label-mean
differences and label-blind ranking agreement. With --criteria regression, also
applies the pre-declared attention-policy acceptance criteria
(results/regression/ACCEPTANCE_CRITERIA.md).
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_scores import block_features  # noqa: E402

QUANTILES = [0.5, 0.9, 0.99, 0.999]


def load_runs(dirs):
    """{(record, orientation): (track, run_meta, record_meta)} across run dirs."""
    out = {}
    for d in dirs:
        for run_json in sorted(Path(d).glob('**/run.json')):
            meta = json.loads(run_json.read_text())
            for rec in meta['records']:
                key = (rec['id'], rec['scoring_orientation'])
                if key in out:
                    raise ValueError(f'duplicate track {key} in {run_json}')
                v = np.load(run_json.parent / rec['track'])['above_uniform'].astype(np.float64)
                out[key] = (v, meta, rec)
    return out


def ranks(values):
    order = sorted(values, key=lambda k: -values[k])
    return {k: i + 1 for i, k in enumerate(order)}


def ranking_changes(ref, cand, tie=0.01):
    """Pairs whose order flips by more than the tie tolerance."""
    flips = []
    keys = list(ref)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            dr, dc = ref[a] - ref[b], cand[a] - cand[b]
            if dr * dc < 0 and abs(dr) > tie and abs(dc) > tie:
                flips.append((a, b))
    return flips


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ref', type=Path, nargs='+', required=True)
    p.add_argument('--cand', type=Path, nargs='+', required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--criteria', choices=['none', 'regression'], default='none')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = {r['record_id']: r for r in csv.DictReader(a.manifest.open(), delimiter='\t')}
    ref, cand = load_runs(a.ref), load_runs(a.cand)
    keys = sorted(set(ref) & set(cand))
    missing = sorted(set(ref) ^ set(cand))
    if not keys:
        raise SystemExit('No overlapping tracks.')

    rows, all_diff, feat_diffs, label_rows = [], [], {}, []
    win_ref, win_cand = {}, {}
    for key in keys:
        r, rm, rr = ref[key]
        c, cm, cr = cand[key]
        if len(r) != len(c) or rr['sequence_sha256'] != cr['sequence_sha256']:
            raise ValueError(f'{key}: sequence/length mismatch')
        d = np.abs(c - r)
        all_diff.append(d)
        row = {'record_id': key[0], 'orientation': key[1], 'n': len(d),
               'finite': bool(np.isfinite(c).all()),
               'n_variants_ref': rr['n_variants'], 'n_variants_cand': cr['n_variants'],
               'max_abs': float(d.max()), 'mean_abs': float(d.mean()),
               'argmax': int(d.argmax()), 'corr': float(np.corrcoef(r, c)[0, 1])}
        for q in QUANTILES:
            row[f'q{q}'] = float(np.quantile(d, q))
        rows.append(row)
        fr, fc = block_features(r), block_features(c)
        for f in fr:
            if f.endswith('_n'):
                continue
            feat_diffs.setdefault(f, []).append(abs(fc[f] - fr[f]))
            win_ref.setdefault(f, {}).setdefault(key[0], []).append(fr[f])
            win_cand.setdefault(f, {}).setdefault(key[0], []).append(fc[f])
        m = man[key[0]]
        if m['label_start_in_window0']:
            ls, le = int(m['label_start_in_window0']), int(m['label_end_in_window0'])
            label_rows.append({'record_id': key[0], 'orientation': key[1],
                               'family': m['family'], 'ref': float(r[ls:le].mean()),
                               'cand': float(c[ls:le].mean()),
                               'abs_diff': float(abs(c[ls:le].mean() - r[ls:le].mean()))})

    alld = np.concatenate(all_diff)
    overall = {'n_positions': int(alld.size), 'max_abs': float(alld.max()),
               'mean_abs': float(alld.mean()),
               **{f'q{q}': float(np.quantile(alld, q)) for q in QUANTILES}}
    feat_summary = {f: float(max(v)) for f, v in feat_diffs.items()}
    # window score = max over orientations (as in calibration), then rank
    flips = {}
    for f in win_ref:
        wr = {k: max(v) for k, v in win_ref[f].items()}
        wc = {k: max(v) for k, v in win_cand[f].items()}
        fl = ranking_changes(wr, wc)
        if fl:
            flips[f] = fl

    with (a.out / 'per_track.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    verdict = None
    checks = {}
    if a.criteria == 'regression':
        checks = {
            '1_perbase_p99<=0.03': overall['q0.99'] <= 0.03,
            '1_perbase_max<=0.15': overall['max_abs'] <= 0.15,
            '2_features<=0.01': max(feat_summary.values()) <= 0.01,
            '3_label_mean<=0.005': all(x['abs_diff'] <= 0.005 for x in label_rows),
            '4_ranking_unchanged': not flips,
            '5_finite_and_coverage': all(x['finite'] and x['n_variants_ref'] == x['n_variants_cand']
                                         for x in rows) and not missing,
        }
        verdict = 'PASS' if all(checks.values()) else 'FAIL'

    def policy(meta):
        return meta.get('attention_mask_policy', 'always (pre-policy code)')
    any_r, any_c = ref[keys[0]][1], cand[keys[0]][1]
    lines = ['# Track comparison', '',
             f"reference: {', '.join(map(str, a.ref))}  policy={policy(any_r)} "
             f"mode={any_r['mode']} batch={any_r.get('batch_size')} dtype={any_r.get('dtype')}",
             f"candidate: {', '.join(map(str, a.cand))}  policy={policy(any_c)} "
             f"mode={any_c['mode']} batch={any_c.get('batch_size')} dtype={any_c.get('dtype')}",
             f"model revision: ref {any_r['revision']} / cand {any_c['revision']}",
             f"candidate code: {json.dumps(any_c.get('code_sha256'))}",
             f"tracks compared: {len(keys)}; missing: {missing or 'none'}", '']
    if verdict:
        lines += [f'## Verdict: **{verdict}**', ''] + \
                 [f"- {k}: {'ok' if v else 'FAILED'}" for k, v in checks.items()] + ['']
    lines += ['## Per-base |diff| (nats), all positions', '',
              '| n | mean | p50 | p90 | p99 | p99.9 | max |', '|---:|---:|---:|---:|---:|---:|---:|',
              f"| {overall['n_positions']} | {overall['mean_abs']:.5f} | {overall['q0.5']:.5f} | "
              f"{overall['q0.9']:.5f} | {overall['q0.99']:.5f} | {overall['q0.999']:.5f} | "
              f"{overall['max_abs']:.5f} |", '',
              '## Per track', '',
              '| record | orient | mean | p99 | max | argmax | corr |',
              '|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['record_id']} | {r['orientation']} | {r['mean_abs']:.5f} | "
                     f"{r['q0.99']:.5f} | {r['max_abs']:.5f} | {r['argmax']} | {r['corr']:.5f} |")
    lines += ['', '## Max |diff| of window/block features (over tracks)', '',
              '| feature | max abs diff |', '|---|---:|']
    lines += [f'| {f} | {v:.5f} |' for f, v in feat_summary.items()]
    lines += ['', '## Label-region mean', '', '| record | family | orient | ref | cand | abs diff |',
              '|---|---|---|---:|---:|---:|']
    lines += [f"| {x['record_id']} | {x['family']} | {x['orientation']} | {x['ref']:.4f} | "
              f"{x['cand']:.4f} | {x['abs_diff']:.5f} |" for x in label_rows]
    lines += ['', '## Label-blind ranking (window score = max over orientations)', '',
              'Order flips beyond 0.01 tie tolerance: ' +
              ('none' if not flips else json.dumps({k: v for k, v in flips.items()}))]
    (a.out / 'report.md').write_text('\n'.join(lines) + '\n')
    (a.out / 'summary.json').write_text(json.dumps(
        {'verdict': verdict, 'checks': checks, 'overall': overall,
         'feature_max_abs_diff': feat_summary, 'label': label_rows,
         'ranking_flips': flips, 'missing': missing}, indent=2) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
