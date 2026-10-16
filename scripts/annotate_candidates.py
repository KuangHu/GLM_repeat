#!/usr/bin/env python3
"""Annotate scan candidates AFTER scoring, for description only.

Annotation never gates, filters, penalises or re-ranks a candidate. Every candidate
is retained, including those matching nothing: `unclassified` is a reported
category, not a discard. Categories come from the frozen GenBank records of the
scanned genomes, plus the manuscript's 14-nt ART core motif (marked motif-derived).
"""
import argparse
import csv
import re
from pathlib import Path

ART_CORE = 'ATATGAATACGTAT'
FEATURES = ('tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'repeat_region', 'mobile_element')


def parse_features(path):
    """GenBank features of interest as (kind, start0, end0, qualifier blob)."""
    text = path.read_text()
    out = []
    pat = (r'^     (' + '|'.join(FEATURES) + r')\s+(complement\()?(?:join\()?'
           r'(\d+)\.\.(?:\d+\.\.)?(\d+)\)?\)?\s*\n((?:^(?!     \S).*\n)*)')
    for m in re.finditer(pat, text, re.M):
        out.append((m.group(1), int(m.group(3)) - 1, int(m.group(4)), m.group(5)))
    return out


def gb_sequence(path):
    out, on = [], False
    for line in path.open():
        if line.startswith('ORIGIN'):
            on = True
            continue
        if on:
            if line.startswith('//'):
                break
            out.append(re.sub(r'[^acgtACGTnN]', '', line))
    return ''.join(out).upper()


def art_core_spans(seq, max_mm=1):
    comp = ART_CORE.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    hits, i = [], 0
    while i <= len(seq) - len(ART_CORE):
        for q in (ART_CORE, comp):
            if sum(a != b for a, b in zip(seq[i:i + len(ART_CORE)], q)) <= max_mm:
                hits.append((i, i + len(ART_CORE)))
                i += len(ART_CORE) - 1
                break
        i += 1
    return hits


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def categorise(feats, cores, s, e):
    hit = [(k, fs, fe, q) for k, fs, fe, q in feats if overlap(s, e, fs, fe) > 0]
    crispr = [h for h in hit if h[0] == 'repeat_region' and 'CRISPR' in h[3]]
    trna = [h for h in hit if h[0] == 'tRNA']
    rrna = [h for h in hit if h[0] == 'rRNA']
    other_rpt = [h for h in hit if h[0] == 'repeat_region' and 'CRISPR' not in h[3]]
    other_rna = [h for h in hit if h[0] in ('ncRNA', 'tmRNA')]
    ncore = sum(1 for cs, ce in cores if overlap(s, e, cs, ce) > 0)
    if crispr:
        cat = 'CRISPR_array'
    elif ncore >= 3:
        cat = 'ART_core_array_motif'
    elif len(trna) >= 2:
        cat = 'tRNA_cluster'
    elif rrna:
        cat = 'rRNA'
    elif trna:
        cat = 'tRNA_single'
    elif other_rpt:
        cat = 'other_annotated_repeat'
    elif other_rna:
        cat = 'other_annotated_ncRNA'
    elif ncore:
        cat = 'ART_core_partial_motif'
    else:
        cat = 'unclassified'
    detail = ';'.join(f'{k}:{fs}-{fe}' for k, fs, fe, _ in hit[:6])
    return cat, len(hit), ncore, detail


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates', type=Path, required=True)
    p.add_argument('--gb-dirs', type=Path, nargs='+', required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cands = list(csv.DictReader(a.candidates.open(), delimiter='\t'))
    contigs = sorted({c['contig'] for c in cands})
    feats, cores = {}, {}
    for contig in contigs:
        gb = next((d / f'{contig}.gb' for d in a.gb_dirs
                   if (d / f'{contig}.gb').exists()), None)
        if gb is None:
            print(f'{contig}: no GenBank record; candidates stay unclassified')
            feats[contig], cores[contig] = [], []
            continue
        feats[contig] = parse_features(gb)
        cores[contig] = art_core_spans(gb_sequence(gb))
        print(f'{contig}: {len(feats[contig])} annotated features, '
              f'{len(cores[contig])} ART core motif hits')
    rows, counts = [], {}
    for c in cands:
        cat, nfeat, ncore, detail = categorise(
            feats[c['contig']], cores[c['contig']], int(c['start0']), int(c['end0']))
        counts[cat] = counts.get(cat, 0) + 1
        rows.append({**c, 'annotation_category': cat, 'n_overlapping_features': nfeat,
                     'n_overlapping_art_cores': ncore, 'overlapping_features': detail})
    with (a.out / 'candidates_annotated.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: -float(r['score'])))
    print('\ncategory counts (annotation is descriptive; nothing was filtered):')
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f'  {k:28s} {v}')


if __name__ == '__main__':
    main()
