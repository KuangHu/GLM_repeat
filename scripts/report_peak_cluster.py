#!/usr/bin/env python3
"""Report peak-cluster results: candidates, peak support, label localization and
window ranking, each in its own section, plus optional RP-vs-exact comparison.

Labels are read only in the localization section and in the ranking read-out; they
never enter candidate formation or scoring.
"""
import argparse
import csv
import json
from pathlib import Path


def read(p):
    return list(csv.DictReader(p.open(), delimiter='\t')) if p.exists() else []


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def jaccard(a, b):
    """Reciprocal overlap of two merged candidate sets, by genomic span."""
    hit = 0
    for x in a:
        for y in b:
            if x['contig'] == y['contig'] and overlap(
                    int(x['start0']), int(x['end0']), int(y['start0']), int(y['end0'])) > 0:
                hit += 1
                break
    return hit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exact', type=Path, required=True)
    p.add_argument('--rp', type=Path, nargs='*', default=[])
    p.add_argument('--bench', type=Path)
    p.add_argument('--manifest', type=Path,
                   help='window manifest; defaults to <bench>/data/seed19_manifest.tsv')
    p.add_argument('--baselines', type=Path,
                   help='window_scores.tsv from report_seed_panel.py (block baselines)')
    p.add_argument('--annotations', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    mpath = a.manifest or (a.bench / 'data/seed19_manifest.tsv')
    man = {r['record_id']: r for r in csv.DictReader(mpath.open(), delimiter='\t')}
    run = json.loads((a.exact / 'run.json').read_text())
    cands = read(a.exact / 'candidates_merged.tsv')
    peaks = read(a.exact / 'peaks.tsv')
    percl = read(a.exact / 'clusters_per_window.tsv')
    wscore = {r['record_id']: float(r['peak_cluster_score'])
              for r in read(a.exact / 'window_scores.tsv')}
    base = {r['record_id']: r for r in read(a.baselines)} if a.baselines else {}
    extra = {r['record_id']: r for r in read(a.annotations)} if a.annotations else {}

    L = ['# Peak-cluster scorer on the seed panel', '',
         'Scorer design is agent-authored and was frozen before any track was scored: '
         'results/scorer/PEAK_CLUSTER_DESIGN.md.',
         f"Parameters: {json.dumps(run['params'])}", '',
         f"Tracks: {run['n_tracks']} ({run['track_mode']} mode, model "
         f"{run['model_revision'][:10]}, policy {run['attention_mask_policy']}).",
         'seed19 is development data: no FDR, no calibrated threshold, no recall claim.',
         '',
         '## 1. Candidate coordinates (merged across windows, genomic, 0-based half-open)',
         '', '| contig | start0 | end0 | span | score | peaks | median peak mean | orient | windows |',
         '|---|---:|---:|---:|---:|---:|---:|---|---|']
    for c in cands:
        L.append(f"| {c['contig']} | {c['start0']} | {c['end0']} | {c['span_nt']} | "
                 f"{c['score']} | {c['n_peaks']} | {c['median_peak_mean']} | "
                 f"{c['orientation']} | {c['source_windows']} |")
    L += ['', f'{len(percl)} per-window clusters merged into {len(cands)} candidates '
          f'(overlapping windows never counted as independent peaks).', '',
          '## 2. Peak support', '',
          '| candidate | peaks | peak lengths (nt) | gaps between peaks (nt) |',
          '|---|---:|---|---|']
    for c in cands:
        L.append(f"| {c['contig']}:{c['start0']}-{c['end0']} | {c['n_peaks']} | "
                 f"{c['peak_lengths']} | {c['gaps'] or '-'} |")
    npk = {}
    for pk in peaks:
        npk[pk['window_id']] = npk.get(pk['window_id'], 0) + 1
    L += ['', f'{len(peaks)} raw peaks over all windows and orientations; '
          f'median per window-orientation '
          f'{sorted(npk.values())[len(npk) // 2] if npk else 0}.', '']

    # ---- 3. label localization (labels enter only here)
    L += ['## 3. Label localization (labels used only to read the result)', '',
          '| record | family | label | best overlapping candidate | overlap nt | '
          'fraction of label | fraction of candidate |', '|---|---|---|---|---:|---:|---:|']
    recovered = 0
    labelled = [r for r in man.values() if r['label_start0']]
    for m in labelled:
        ls, le = int(m['label_start0']), int(m['label_end0'])
        best, bo = None, 0
        for c in cands:
            if c['contig'] != m['accession_version']:
                continue
            o = overlap(ls, le, int(c['start0']), int(c['end0']))
            if o > bo:
                best, bo = c, o
        recovered += bo > 0
        if best:
            L.append(f"| {m['record_id']} | {m['family']} | {ls}-{le} | "
                     f"{best['start0']}-{best['end0']} | {bo} | {bo / (le - ls):.2f} | "
                     f"{bo / (int(best['end0']) - int(best['start0'])):.2f} |")
        else:
            L.append(f"| {m['record_id']} | {m['family']} | {ls}-{le} | "
                     f"**none** | 0 | 0.00 | - |")
    L += ['', f'{recovered}/{len(labelled)} labels overlapped by a candidate. '
          'Overlap is loose recovery, not boundary reconstruction.', '',
          'Candidates in background windows (no label; not false positives by construction '
          '- backgrounds are unlabeled real DNA):', '']
    for c in cands:
        rid = c['best_window']
        if man[rid]['family'] == 'background':
            note = extra.get(rid, {}).get('annotation', '')
            L.append(f"- {c['contig']}:{c['start0']}-{c['end0']} score {c['score']} "
                     f"in {rid} ({man[rid]['source_id']})" + (f' - {note}' if note else ''))

    # ---- 4. window ranking
    L += ['', '## 4. Window ranking (peak-cluster vs unchanged block baselines)', '',
          '| rank | record | family | peak-cluster | blk256_mean_max | blk1024_mean_max | '
          'blk256_top20_max | blk1024_top20_max |', '|---:|---|---|---:|---:|---:|---:|---:|']
    order = sorted(wscore, key=lambda r: -wscore[r])
    for i, rid in enumerate(order, 1):
        b = base.get(rid, {})
        L.append(f"| {i} | {rid} | {man[rid]['family']} | {wscore[rid]:.3f} | "
                 + ' | '.join(b.get(k, '-') for k in
                              ['blk256_mean_max', 'blk1024_mean_max',
                               'blk256_top20_max', 'blk1024_top20_max']) + ' |')
    pos = [r for r in wscore if man[r]['family'] != 'background']
    bg = [r for r in wscore if man[r]['family'] == 'background']
    nz = lambda rs: sum(wscore[r] > 0 for r in rs)
    auc = sum((wscore[x] > wscore[y]) + 0.5 * (wscore[x] == wscore[y])
              for x in pos for y in bg) / (len(pos) * len(bg))
    L += ['', f'Non-zero peak-cluster score: {nz(pos)}/{len(pos)} positives, '
          f'{nz(bg)}/{len(bg)} backgrounds. P(positive > background) = {auc:.2f} '
          f'(ties at 0 counted as 0.5; descriptive only at n={len(pos)}).', '']

    # ---- RP comparison, against the NEW scorer's outputs
    if a.rp:
        L += ['## 5. Random partition vs exact, evaluated on this scorer', '',
              '| RP run | candidates | exact candidates recovered | RP candidates not in exact | '
              'window-score Spearman | positives with non-zero score |',
              '|---|---:|---|---:|---:|---:|']
        import numpy as np
        for rp_dir in a.rp:
            rc = read(rp_dir / 'candidates_merged.tsv')
            rw = {r['record_id']: float(r['peak_cluster_score'])
                  for r in read(rp_dir / 'window_scores.tsv')}
            ids = sorted(set(wscore) & set(rw))
            rank = lambda d: np.argsort(np.argsort([-d[i] for i in ids]))
            sp = float(np.corrcoef(rank(wscore), rank(rw))[0, 1])
            L.append(f"| {rp_dir.name} | {len(rc)} | {jaccard(cands, rc)}/{len(cands)} | "
                     f"{len(rc) - jaccard(rc, cands)} | {sp:.3f} | "
                     f"{sum(rw[r] > 0 for r in pos)}/{len(pos)} |")
        lost = []
        for rp_dir in a.rp:
            rc = read(rp_dir / 'candidates_merged.tsv')
            for c in cands:
                if not any(c['contig'] == y['contig'] and overlap(
                        int(c['start0']), int(c['end0']),
                        int(y['start0']), int(y['end0'])) > 0 for y in rc):
                    lost.append((rp_dir.name, c))
        L += ['', ('Exact-supported candidates missed by RP: ' +
                   ('none' if not lost else '')) ]
        for name, c in lost:
            L.append(f"- {name}: {c['contig']}:{c['start0']}-{c['end0']} "
                     f"(score {c['score']}, {c['n_peaks']} peaks, window {c['best_window']}, "
                     f"family {man[c['best_window']]['family']})")
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
