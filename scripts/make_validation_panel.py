#!/usr/bin/env python3
"""Build an independent validation panel from freshly downloaded references.

CRISPR loci come from GenBank `repeat_region /rpt_family="CRISPR"` annotations in
hosts divergent from the seed panel. ART loci are the manuscript-listed
Staphylococcus phages, with coordinates defined here by the manuscript's 14-nt core
rule; the RT-start boundary is NOT applied because these records annotate almost
everything as hypothetical protein. Evidence tiers are recorded per locus and are
not promoted to experimentally mapped RNA boundaries.

Backgrounds are GC-matched windows from the same genome, excluding every locus with
a margin. No score, repeat finder or low-complexity filter selects anything.
"""
import argparse
import csv
import hashlib
import re
from pathlib import Path

import numpy as np

CORE = 'ATATGAATACGTAT'
ART = [('MW349129.1', 'Madawaska'), ('OR836606.1', 'LY01'), ('LC680885.1', 'S6'),
       ('MN091626.1', 'PALS_2'), ('MZ779063.1', 'UFV_DC4')]
CRISPR = [('NC_003413.1', 'Pyrococcus_furiosus_DSM3638', 'Archaea/Euryarchaeota'),
          ('NC_006461.1', 'Thermus_thermophilus_HB8', 'Deinococcus-Thermus')]
ART_GROUP = 'Staphylococcus_ART_type_II'   # same lineage as SA1/MarsHill


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


def gb_crispr_regions(path):
    """repeat_region features carrying /rpt_family="CRISPR" (0-based half-open)."""
    text = path.read_text()
    out = []
    for m in re.finditer(r'^     repeat_region\s+(complement\()?(\d+)\.\.(\d+)\)?\s*\n'
                         r'((?:^(?!     \S).*\n)*)', text, re.M):
        if 'CRISPR' in m.group(4):
            out.append((int(m.group(2)) - 1, int(m.group(3))))
    return sorted(out)


def core_hits(seq, max_mm=1):
    rc = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    del rc
    comp = CORE.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    hits, i = [], 0
    while i <= len(seq) - len(CORE):
        for strand, q in (('+', CORE), ('-', comp)):
            if sum(a != b for a, b in zip(seq[i:i + len(CORE)], q)) <= max_mm:
                hits.append((i, i + len(CORE), strand))
                i += len(CORE) - 1
                break
        i += 1
    return hits


def largest_array(hits, max_gap=600, min_cores=3):
    """Drop isolated singleton core matches; keep the biggest clustered array."""
    if not hits:
        return []
    groups, cur = [], [hits[0]]
    for h in hits[1:]:
        if h[0] - cur[-1][1] <= max_gap:
            cur.append(h)
        else:
            groups.append(cur)
            cur = [h]
    groups.append(cur)
    groups = [g for g in groups if len(g) >= min_cores]
    return max(groups, key=len) if groups else []


def pick_backgrounds(genome, blocked, target_gc, n, window, rng, tol=0.02):
    cands = []
    for s in range(0, len(genome) - window + 1, window):
        if blocked[s:s + window].any():
            continue
        w = genome[s:s + window]
        if set(w) - set('ACGT'):
            continue
        gc = (w.count('G') + w.count('C')) / window
        if abs(gc - target_gc) <= tol:
            cands.append(s)
    if not cands:
        return []
    return sorted(rng.choice(cands, size=min(n, len(cands)), replace=False).tolist())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, required=True, help='directory of .gb downloads')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--window', type=int, default=3000)
    p.add_argument('--bg-per-positive', type=int, default=1)
    p.add_argument('--margin', type=int, default=5000)
    p.add_argument('--seed', type=int, default=20261015)
    a = p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    loci, genomes, subregions = [], {}, []
    for acc, host, phylum in CRISPR:
        g = gb_sequence(a.raw / f'{acc}.gb')
        genomes[acc] = g
        for i, (s, e) in enumerate(gb_crispr_regions(a.raw / f'{acc}.gb'), 1):
            loci.append({'locus_id': f'{host}_CRISPR_{s + 1}', 'family': 'CRISPR',
                         'accession_version': acc, 'source_id': host,
                         'label_start0': s, 'label_end0': e,
                         'independence_group': host, 'taxon_note': phylum,
                         'evidence': 'GenBank repeat_region /rpt_family=CRISPR '
                                     '(computational annotation in a divergent host)'})
    for acc, name in ART:
        g = gb_sequence(a.raw / f'{acc}.gb')
        genomes[acc] = g
        arr = largest_array(core_hits(g))
        if not arr:
            print(f'{name}: no clustered core array; skipped')
            continue
        start = max(0, arr[0][0] - 60)      # manuscript: 60 nt before the first core
        end = arr[-1][1]
        loci.append({'locus_id': f'{name}_ART', 'family': 'ART',
                     'accession_version': acc, 'source_id': name,
                     'label_start0': start, 'label_end0': end,
                     'independence_group': ART_GROUP, 'taxon_note': 'Staphylococcus phage',
                     'evidence': f'manuscript-listed ART phage; extent defined here by '
                                 f'{len(arr)} 14-nt cores (<=1 mismatch), 60 nt before '
                                 f'first core to last core end; RT-start boundary NOT '
                                 f'applied (records annotate only hypothetical proteins)'})
        for k, (s, e, strand) in enumerate(arr, 1):
            subregions.append({'source_id': name, 'accession_version': acc, 'start0': s,
                               'end0': e, 'strand': strand, 'feature': f'ART_core_{k}',
                               'evidence': 'manuscript 14-nt core, <=1 mismatch'})

    rows = []
    with (a.out / 'validation.fa').open('w') as fa:
        for i, loc in enumerate(loci, 1):
            g = genomes[loc['accession_version']]
            mid = (loc['label_start0'] + loc['label_end0']) // 2
            s = max(0, min(mid - a.window // 2, len(g) - a.window))
            w = g[s:s + a.window]
            if set(w) - set('ACGT'):
                print(f"{loc['locus_id']}: window has non-ACGT; skipped")
                continue
            rid = f'v{i:04d}'
            fa.write(f'>{rid}\n{w.lower()}\n')
            ls, le = loc['label_start0'], loc['label_end0']
            rows.append({'record_id': rid, 'family': loc['family'],
                         'locus_id': loc['locus_id'], 'source_id': loc['source_id'],
                         'accession_version': loc['accession_version'],
                         'start0': s, 'end0': s + a.window, 'length_nt': a.window,
                         'label_start0': ls, 'label_end0': le,
                         'label_start_in_window0': max(0, ls - s),
                         'label_end_in_window0': min(a.window, le - s),
                         'label_nt_in_window': min(a.window, le - s) - max(0, ls - s),
                         'label_nt_total': le - ls,
                         'gc_fraction': f'{(w.count("G") + w.count("C")) / a.window:.4f}',
                         'independence_group': loc['independence_group'],
                         'taxon_note': loc['taxon_note'], 'paired_positive_id': '',
                         'split': 'independent_validation', 'evidence': loc['evidence'],
                         'sequence_sha256': hashlib.sha256(w.lower().encode()).hexdigest()})
        n = len(rows)
        for j, r in enumerate(list(rows), 1):
            g = genomes[r['accession_version']]
            blocked = np.zeros(len(g), bool)
            for loc in loci:
                if loc['accession_version'] == r['accession_version']:
                    blocked[max(0, loc['label_start0'] - a.margin):
                            loc['label_end0'] + a.margin] = True
            for prev in rows:
                if prev['accession_version'] == r['accession_version']:
                    blocked[int(prev['start0']):int(prev['end0'])] = True
            for s in pick_backgrounds(g, blocked, float(r['gc_fraction']),
                                      a.bg_per_positive, a.window, rng):
                n += 1
                rid = f'v{n:04d}'
                w = g[s:s + a.window]
                fa.write(f'>{rid}\n{w.lower()}\n')
                rows.append({**{k: '' for k in rows[0]}, 'record_id': rid,
                             'family': 'background', 'source_id': r['source_id'],
                             'accession_version': r['accession_version'],
                             'start0': s, 'end0': s + a.window, 'length_nt': a.window,
                             'gc_fraction': f'{(w.count("G") + w.count("C")) / a.window:.4f}',
                             'independence_group': r['independence_group'],
                             'taxon_note': r['taxon_note'],
                             'paired_positive_id': r['record_id'],
                             'split': 'independent_validation',
                             'evidence': 'unlabeled real DNA; GC-matched; '
                                         'not a confirmed negative',
                             'sequence_sha256': hashlib.sha256(w.lower().encode()).hexdigest()})
    with (a.out / 'validation_manifest.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    if subregions:
        with (a.out / 'validation_subregions.tsv').open('w', newline='') as fh:
            w = csv.DictWriter(fh, list(subregions[0]), delimiter='\t')
            w.writeheader()
            w.writerows(subregions)
    (a.out / 'VIPR_STATUS.txt').write_text(
        'VIPR: missing_reference. No independent natural VIPR locus was added. The '
        'public VIPR analysis archive (zenodo 20945944) states it does not include '
        'the frozen study datasets, and no confidently homologous vrRNA was located '
        'outside the SUSP1 seed (best non-self match 63% identity, inside the '
        'unconfirmed off-label candidate). VIPR is therefore NOT tested in this '
        'independent validation round.\n')
    pos = [r for r in rows if r['family'] != 'background']
    print(f"{len(rows)} windows: {len(pos)} positives "
          f"({sum(r['family'] == 'CRISPR' for r in pos)} CRISPR, "
          f"{sum(r['family'] == 'ART' for r in pos)} ART, 0 VIPR), "
          f"{len(rows) - len(pos)} backgrounds -> {a.out}")
    for r in pos:
        print(f"  {r['record_id']} {r['family']:7s} {r['locus_id']:34s} "
              f"label {r['label_nt_total']}nt, {r['label_nt_in_window']}nt in window")


if __name__ == '__main__':
    main()
