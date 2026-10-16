#!/usr/bin/env python3
"""Held-out pilot report.

Terminology used throughout:
  held-out loci         known loci in genomes never used for calibration, reported
                        with their genome / family grouping
  localization precision  fraction of a candidate interval lying inside the locus
  locus coverage          fraction of a locus interval covered by candidates
  candidate yield         merged candidates per unique scanned genomic Mb
                        (a yield, never a false-discovery rate)
CRISPR, ART and VIPR outcomes are reported separately and never pooled.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def union_len(spans):
    if not spans:
        return 0
    spans = sorted(spans)
    tot, cs, ce = 0, *spans[0]
    for s, e in spans[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs
            cs, ce = s, e
    return tot + ce - cs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scan', type=Path, required=True, help='pilot_scan scan dir')
    p.add_argument('--annotated', type=Path, help='candidates_annotated.tsv')
    p.add_argument('--loci', type=Path, required=True, help='known_loci.tsv')
    p.add_argument('--calibration-scan', type=Path, help='calibration scan dir')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--vipr-note', type=Path)
    p.add_argument('--geometry', type=Path,
                   help='geometry_<tier>.tsv from audit_candidates.py')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((a.scan / 'scan.json').read_text())
    cands = list(csv.DictReader((a.annotated or a.scan / 'candidates.tsv').open(),
                                delimiter='\t'))
    loci = list(csv.DictReader(a.loci.open(), delimiter='\t'))
    geom = {}
    if a.geometry and a.geometry.exists():
        for g in csv.DictReader(a.geometry.open(), delimiter='\t'):
            geom[(g['contig'], int(g['start0']), int(g['end0']))] = g

    L = ['# Held-out pilot scan', '',
         f"Configuration frozen before scanning: results/pilot/PILOT_CONFIG_FROZEN.md. "
         f"Mode {meta['mode']} (groups {meta['groups']}, repeats {meta['repeats']}, "
         f"seed {meta['seed']}), policy {meta['attention_mask_policy']}, model "
         f"{meta['model_revision'][:10]}.",
         f"Threshold {meta['threshold']} ({meta['threshold_label']}), calibrated on "
         'calibration material only and frozen before this scan.', '',
         '## Scan summary', '',
         f"- windows scored: {meta['n_windows_scored']}"
         + (f" (missing {meta['n_windows_missing']}, excluded from scanned nt, "
            'not counted as negatives)' if meta['n_windows_missing'] else ''),
         f"- unique scanned genomic span: {meta['unique_scanned_Mb']} Mb",
         f"- merged candidates: {meta['n_candidates']} "
         f"(from {meta['n_per_window_clusters']} per-window clusters)",
         f"- **candidate yield: {meta['candidates_per_Mb']} per unique scanned Mb** "
         '(a yield, not an FDR)', '']
    if a.calibration_scan and (a.calibration_scan / 'scan.json').exists():
        c = json.loads((a.calibration_scan / 'scan.json').read_text())
        L += [f"- calibration material at the same threshold: "
              f"{c['candidates_per_Mb']}/Mb over {c['unique_scanned_Mb']} Mb", '']
    L += ['### Yield per genome (densities are not forced to match)', '',
          '| contig | unique Mb | candidates/Mb |', '|---|---:|---:|']
    for contig, d in sorted(meta['per_contig_candidates_per_Mb'].items(),
                            key=lambda kv: -kv[1]):
        L.append(f"| {contig} | {meta['per_contig_unique_nt'][contig] / 1e6:.3f} | {d} |")

    # ---- held-out loci, grouped
    L += ['', '## Held-out loci (grouped by genome / family)', '',
          'Overlap alone is not evidence of precise localization, so peak support '
          'inside the best-overlapping candidate is reported alongside it.', '',
          '| locus | family | family group | locus nt | candidates overlapping | '
          'locus coverage | best localization precision | peak support of that '
          'candidate | its peak blocks | max internal peak gap |',
          '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    by_family = defaultdict(list)
    for lo in loci:
        ls, le = int(lo['start0']), int(lo['end0'])
        hit = [c for c in cands if c['contig'] == lo['contig']
               and ov(ls, le, int(c['start0']), int(c['end0'])) > 0]
        cov = union_len([(max(ls, int(c['start0'])), min(le, int(c['end0'])))
                         for c in hit]) / (le - ls) if hit else 0.0
        bestc, prec = None, 0.0
        for c in hit:
            pr = ov(ls, le, int(c['start0']), int(c['end0'])) / (
                int(c['end0']) - int(c['start0']))
            if pr >= prec:
                bestc, prec = c, pr
        g = geom.get((bestc['contig'], int(bestc['start0']), int(bestc['end0']))) if bestc else None
        ps = float(g['peak_support_fraction']) if g else float('nan')
        blocks = g['peak_blocks_in_span'] if g else '-'
        gap = g['max_internal_peak_gap_nt'] if g else '-'
        by_family[lo['family']].append((lo, bool(hit), cov, prec, ps))
        L.append(f"| {lo['locus_id']} | {lo['family']} | {lo['family_group']} | "
                 f"{le - ls} | {len(hit)} | {cov:.2f} | {prec:.2f} | "
                 f"{ps:.2f} | {blocks} | {gap} |")
    L += ['', '### Outcome by family (reported separately, never pooled)', '',
          '| family | family groups | loci | loci with >=1 overlapping candidate | '
          'median locus coverage | median localization precision | '
          'median peak support |', '|---|---:|---:|---:|---:|---:|---:|']
    for fam in ('CRISPR', 'ART', 'VIPR'):
        rows = by_family.get(fam, [])
        if not rows:
            L.append(f'| {fam} | 0 | 0 | not tested | - | - |')
            continue
        groups = len({r[0]['family_group'] for r in rows})
        med = lambda xs: sorted(xs)[len(xs) // 2]
        L.append(f"| {fam} | {groups} | {len(rows)} | {sum(r[1] for r in rows)} | "
                 f"{med([r[2] for r in rows]):.2f} | {med([r[3] for r in rows]):.2f} | "
                 f"{med([r[4] for r in rows]):.2f} |")
    L += ['', 'Loci within one family group are not independent: the Pyrococcus CRISPR '
          'arrays share a genome, and all ART loci share one Staphylococcus lineage. '
          'Counts of loci overstate the number of independent observations.', '']

    # ---- annotation categories
    if a.annotated:
        L += ['## Annotation of candidates (applied after scoring, never a filter)', '',
              '| category | candidates | median score | median span nt |',
              '|---|---:|---:|---:|']
        cats = defaultdict(list)
        for c in cands:
            cats[c['annotation_category']].append(c)
        for cat, rows in sorted(cats.items(), key=lambda kv: -len(kv[1])):
            sc = sorted(float(r['score']) for r in rows)
            sp = sorted(int(r['span_nt']) for r in rows)
            L.append(f'| {cat} | {len(rows)} | {sc[len(sc) // 2]:.3f} | '
                     f'{sp[len(sp) // 2]} |')
        unk = cats.get('unclassified', [])
        L += ['', f'`unclassified` candidates are retained and reported ({len(unk)} of '
              f'{len(cands)}); they are not discarded and are not implied to be '
              'artefacts. Highest-scoring unclassified candidates:', '']
        for c in sorted(unk, key=lambda r: -float(r['score']))[:10]:
            L.append(f"- {c['contig']}:{c['start0']}-{c['end0']} score {c['score']}, "
                     f"{c['n_peaks']} peaks, span {c['span_nt']} nt")

    # ---- VIPR
    L += ['', '## VIPR, reported separately', '',
          'No VIPR locus is present in the held-out material: no independent natural '
          'VIPR locus could be obtained (recorded `missing_reference`). **VIPR '
          'generalization is unresolved by this pilot and is not claimed.**', '',
          'The original VIPR 60-nt YNNGG region is development data (seed19 w0003, '
          'NC_028808.2:[29169,29229)) and is reported apart from any off-label hit. '
          'The off-label candidate NC_028808.2:[27682,27906) remains a '
          'sequence-supported repeat candidate with vrRNA identity unconfirmed; it is '
          'not counted as a VIPR recovery anywhere in this report.']
    if a.vipr_note and a.vipr_note.exists():
        L += ['', f'See {a.vipr_note} for the copy-dependence outcome on that region.']
    L += ['', '## Limits', '',
          '- Candidate yields are yields. No FDR, no calibrated biological probability.',
          '- Calibration genomes are ordinary DNA and may contain real repeats; they '
          'are not verified negatives.',
          '- Held-out loci are database or manuscript-derived annotations, not '
          'experimentally mapped RNA boundaries.',
          '- Copy dependence was not used for candidate eligibility.']
    (a.out / 'report.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
