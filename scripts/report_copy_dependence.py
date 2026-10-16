#!/usr/bin/env python3
"""Summarize the copy-dependence diagnostic.

For each target: mean target score under `original`, `copy` (other copies mutated)
and `noncopy` (matched non-copy sequence mutated, same nt budget). The contrast of
interest is

    copy_dependence = (copy - original) - (noncopy - original) = copy - noncopy

i.e. how much more the target loses when its *copies* are destroyed than when the
same amount of unrelated sequence is destroyed. Seeds give the spread.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean, pstdev


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, required=True, help='dir containing scores.tsv')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((a.scores / 'scores.tsv').open(), delimiter='\t'))
    run = json.loads((a.scores / 'score_run.json').read_text())
    design = json.loads((a.scores / 'design.json').read_text())

    # key: (locus, target, orientation, arm) -> [means over seeds]
    acc = defaultdict(list)
    meta = {}
    for r in rows:
        k = (r['locus'], int(r['target_index']), r['orientation'], r['arm'])
        acc[k].append(float(r['target_mean']))
        meta[(r['locus'], int(r['target_index']))] = (
            int(r['target_start']), int(r['target_end']), int(r['n_other_copies']))

    L = ['# Copy-dependence diagnostic', '',
         'Target unit and its local flank are byte-identical in every arm; only the '
         'unchanged target positions are re-scored (exact single-position masking, '
         'both orientations).',
         f"Model {run['revision'][:10]}, dtype {run['dtype']}, policy "
         f"{run['attention_mask_policy']}. Seeds {design['seeds']}. "
         f"Mutation: {design['mutation']}.",
         '', 'DIAGNOSTIC ONLY: unit definitions are annotation-derived and are never '
         'inputs to the discovery scorer.', '',
         '## Per target (mean over seeds, nats above uniform)', '',
         '| locus | target | other copies | orient | original | copy | noncopy | '
         'copy-original | noncopy-original | copy-dependence |',
         '|---|---:|---:|---|---:|---:|---:|---:|---:|---:|']
    per_locus = defaultdict(list)
    for (locus, t), (ts, te, nother) in sorted(meta.items()):
        for o in '+-':
            orig = acc[(locus, t, o, 'original')]
            cp = acc[(locus, t, o, 'copy')]
            nc = acc[(locus, t, o, 'noncopy')]
            if not (orig and cp and nc):
                continue
            b, c, n = mean(orig), mean(cp), mean(nc)
            L.append(f'| {locus} | {t} | {nother} | {o} | {b:+.3f} | {c:+.3f} | {n:+.3f} | '
                     f'{c - b:+.3f} | {n - b:+.3f} | **{c - n:+.3f}** |')
            per_locus[locus].append((c - b, n - b, c - n,
                                     pstdev(cp) if len(cp) > 1 else 0.0,
                                     pstdev(nc) if len(nc) > 1 else 0.0))
    L += ['', '## Per locus (mean over targets and orientations)', '',
          '| locus | targets x orient | copy-original | noncopy-original | '
          'copy-dependence | seed SD (copy) | seed SD (noncopy) |',
          '|---|---:|---:|---:|---:|---:|---:|']
    for locus, v in sorted(per_locus.items(), key=lambda kv: mean(x[2] for x in kv[1])):
        L.append(f'| {locus} | {len(v)} | {mean(x[0] for x in v):+.3f} | '
                 f'{mean(x[1] for x in v):+.3f} | **{mean(x[2] for x in v):+.3f}** | '
                 f'{mean(x[3] for x in v):.3f} | {mean(x[4] for x in v):.3f} |')
    L += ['', 'A strongly negative copy-dependence means the target unit is predictable '
          'largely *because* its other copies are present. A value near zero means the '
          'target is predictable on its own, or equally sensitive to any perturbation.',
          '', 'Caveats: the non-copy arm matches perturbed-nt budget and approximate '
          'distance to the target, not local composition or annotation status. '
          'Distance matching widens when placement is tight (see '
          '`distance_mismatch_max` in scores.tsv). Target counts are small and these '
          'loci are not independent (SA1/MarsHill share a lineage).']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
