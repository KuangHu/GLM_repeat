#!/usr/bin/env python3
"""Consolidated seed-panel report: label-level signal and label-blind ranking.

Window score = max over the two scoring orientations. Ranking features are fixed
here (RANK_FEATURES) and were declared before the 19-window panel was scored; no
threshold is selected. Backgrounds are unlabeled real DNA, not confirmed negatives.
"""
import argparse
import csv
import json
import sys
from pathlib import Path
from statistics import median

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_tracks import load_runs  # noqa: E402
from summarize_scores import block_features  # noqa: E402

RANK_FEATURES = ['mean', 'q90', 'blk256_mean_max', 'blk1024_mean_max',
                 'blk256_top20_max', 'blk1024_top20_max']


def auc(pos, neg):
    """P(score_pos > score_bg), ties 0.5. Descriptive only at this n."""
    s = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return s / (len(pos) * len(neg))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, nargs='+', required=True)
    p.add_argument('--bench', type=Path, required=True, help='GLM_REPEAT_REAL_SEEDS_V1 dir')
    p.add_argument('--audit', type=Path, help='off-label candidate TSV (reported separately)')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = {r['record_id']: r for r in
           csv.DictReader((a.bench / 'data/seed19_manifest.tsv').open(), delimiter='\t')}
    subs = list(csv.DictReader((a.bench / 'data/annotated_repeat_subregions.tsv').open(),
                               delimiter='\t'))
    tracks = load_runs(a.scores)
    ids = sorted({k[0] for k in tracks})
    missing = [(r, o) for r in man for o in '+-' if (r, o) not in tracks]

    prov = sorted({(m['revision'], m.get('dtype'), m.get('attention_mask_policy', 'always'),
                    m['mode'], m.get('batch_size'),
                    json.dumps(m.get('code_sha256'), sort_keys=True) if m.get('code_sha256')
                    else 'not recorded: bundle v0.1 reference_score.py sha256 012c0e1a...')
                   for _, m, _ in tracks.values()})
    L = ['# Seed panel (19 windows): exact masked scoring', '',
         'Development/smoke panel only: not calibration or test data. Two CRISPR loci share',
         'the LMD9 genome; both ART loci share one Staphylococcus phage lineage. Backgrounds',
         'are unlabeled real DNA (not confirmed negatives). Scores are masked predictability',
         '(nats above uniform), not repeat probabilities. No thresholds were selected.', '',
         '## Provenance', '']
    for rev, dt, pol, mode, bs, code in prov:
        L.append(f'- model {rev}, dtype {dt}, attention_mask_policy {pol}, mode {mode}, '
                 f'batch {bs}, code {code}')
    L += [f'- tracks: {len(tracks)} ({len(ids)} windows x orientations); missing: '
          f'{missing or "none"}', '']

    inf = [rec.get('inference_seconds') for _, _, rec in tracks.values()]
    peak = [rec.get('peak_gpu_mb') for _, _, rec in tracks.values()]
    L += ['## Timing', '',
          f'- per 3,000-nt orientation: median {median(inf):.0f} s '
          f'(range {min(inf):.0f}-{max(inf):.0f}); peak GPU {max(peak):.0f} MB',
          f'- total GPU inference {sum(inf) / 3600:.2f} h for {len(tracks)} orientations', '']

    feats = {}
    for (rid, o), (v, _, _) in tracks.items():
        feats.setdefault(rid, {})[o] = block_features(v)
    win = {rid: {f: max(fo[o][f] for o in fo) for f in RANK_FEATURES}
           for rid, fo in feats.items()}

    # ---- label-level signal
    L += ['## Label-level signal (uses labels)', '',
          '| record | family | locus | orient | label nt | label mean | outside mean | difference |',
          '|---|---|---|---|---:|---:|---:|---:|']
    for rid in ids:
        m = man[rid]
        if not m['label_start_in_window0']:
            continue
        ls, le = int(m['label_start_in_window0']), int(m['label_end_in_window0'])
        for o in '+-':
            v = tracks[(rid, o)][0]
            inside = np.zeros(len(v), bool)
            inside[ls:le] = True
            L.append(f"| {rid} | {m['family']} | {m['locus_id']} | {o} | {le - ls} | "
                     f"{v[inside].mean():.3f} | {v[~inside].mean():.3f} | "
                     f"{v[inside].mean() - v[~inside].mean():+.3f} |")
    L += ['', 'Annotated sub-regions (orientation-averaged track):', '',
          '| record | feature | genomic [start0,end0) | sub-region mean | rest of label mean |',
          '|---|---|---|---:|---:|']
    for rid in ids:
        m = man[rid]
        if not m['label_start_in_window0']:
            continue
        w0, ls, le = int(m['start0']), int(m['label_start_in_window0']), int(m['label_end_in_window0'])
        v = (tracks[(rid, '+')][0] + tracks[(rid, '-')][0]) / 2
        mine = [s for s in subs if s['accession_version'] == m['accession_version']
                and w0 + ls <= int(s['start0']) and int(s['end0']) <= w0 + le]
        if not mine:
            continue
        sub_mask = np.zeros(len(v), bool)
        for s in mine:
            sub_mask[int(s['start0']) - w0:int(s['end0']) - w0] = True
        lab = np.zeros(len(v), bool)
        lab[ls:le] = True
        feature = mine[0]['feature'].rsplit('_', 1)[0] if len(mine) > 1 else mine[0]['feature']
        span = (f"{mine[0]['start0']}-..-{mine[-1]['end0']} ({len(mine)} parts)"
                if len(mine) > 1 else f"{mine[0]['start0']}-{mine[0]['end0']}")
        L.append(f"| {rid} | {feature} | {span} | {v[sub_mask].mean():.3f} | "
                 f"{v[lab & ~sub_mask].mean():.3f} |")

    # label-blind localisation: where is each window's best 256-nt block?
    L += ['', 'Localisation: position of the best 256-nt block (orientation-averaged track, '
          'chosen label-blind) relative to the label:', '']
    for rid in ids:
        m = man[rid]
        if not m['label_start_in_window0']:
            continue
        ls, le = int(m['label_start_in_window0']), int(m['label_end_in_window0'])
        v = (tracks[(rid, '+')][0] + tracks[(rid, '-')][0]) / 2
        s = max(range(0, len(v) - 255, 64), key=lambda i: v[i:i + 256].mean())
        ov = max(0, min(le, s + 256) - max(ls, s))
        L.append(f"- {rid} ({m['family']}): best block [{s},{s + 256}) mean "
                 f"{v[s:s + 256].mean():.3f}; overlaps label by {ov} nt")

    # ---- label-blind ranking
    pos = [r for r in ids if man[r]['family'] != 'background']
    bg = [r for r in ids if man[r]['family'] == 'background']
    L += ['', '## Label-blind window ranking (labels used only to read the result)', '',
          f'Features fixed in advance: {", ".join(RANK_FEATURES)}. '
          f'{len(pos)} positives, {len(bg)} backgrounds.', '',
          '| record | family | source | ' + ' | '.join(RANK_FEATURES) + ' |',
          '|---|---|---|' + '---:|' * len(RANK_FEATURES)]
    rank = {f: {r: i + 1 for i, r in enumerate(sorted(ids, key=lambda x: -win[x][f]))}
            for f in RANK_FEATURES}
    for rid in sorted(ids, key=lambda x: rank['blk1024_mean_max'][x]):
        m = man[rid]
        cells = [f"{win[rid][f]:.3f} (#{rank[f][rid]})" for f in RANK_FEATURES]
        L.append(f"| {rid} | {m['family']} | {m['source_id']} | " + ' | '.join(cells) + ' |')
    L += ['', '| feature | positive ranks (of 19) | P(pos > bg) | best-ranked background |',
          '|---|---|---:|---|']
    for f in RANK_FEATURES:
        pr = ', '.join(f"{man[r]['family']}:{rank[f][r]}" for r in sorted(pos, key=lambda x: rank[f][x]))
        top_bg = min(bg, key=lambda x: rank[f][x])
        L.append(f"| {f} | {pr} | {auc([win[r][f] for r in pos], [win[r][f] for r in bg]):.2f} | "
                 f"{top_bg} ({man[top_bg]['source_id']}, #{rank[f][top_bg]}) |")

    # ---- paired comparison
    L += ['', '## Positive vs its paired backgrounds (same genome, GC within 2 points)', '',
          '| positive | family | ' + ' | '.join(RANK_FEATURES) + ' |', '|---|---|' + '---:|' * len(RANK_FEATURES)]
    for rid in pos:
        pb = [b for b in bg if man[b]['paired_positive_id'] == rid]
        cells = [f"{win[rid][f]:.3f} vs " + '/'.join(f"{win[b][f]:.3f}" for b in pb)
                 for f in RANK_FEATURES]
        L.append(f"| {rid} vs {','.join(pb)} | {man[rid]['family']} | " + ' | '.join(cells) + ' |')

    if a.audit and a.audit.exists():
        L += ['', '## Separate audit (NOT ground truth, not used in any metric above)', '']
        for r in csv.DictReader(a.audit.open(), delimiter='\t'):
            L.append(f"- {r['candidate_id']}: {r['accession_version']}:[{r['start0']},{r['end0']}) "
                     f"{r['status']}; vrRNA identity {r['vrRNA_identity']}. {r['evidence']}")
    L += ['', 'n is tiny (5 positives, 3 families, 2 independent lineages for CRISPR/ART): '
          'all numbers are descriptive.']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    with (a.out / 'window_scores.tsv').open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['record_id', 'family', 'source_id'] + RANK_FEATURES)
        for rid in ids:
            w.writerow([rid, man[rid]['family'], man[rid]['source_id']] +
                       [f'{win[rid][f]:.5f}' for f in RANK_FEATURES])
    print('\n'.join(L))


if __name__ == '__main__':
    main()
