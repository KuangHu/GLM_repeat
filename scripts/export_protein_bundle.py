#!/usr/bin/env python3
"""Export every panel's predicted proteins with globally unique, traceable IDs.

pyrodigal numbers CDSs per run, so every panel restarts at cds00000001 and the
raw proteins.faa files share an ID space. Pooling them for clustering would merge
unrelated proteins under one identifier, silently and without error. Every protein
is therefore re-identified as `<panel>|<original_id>`, and a mapping table ties
each new ID back to its panel, contig, coordinates, strand and original ID.

The mapping also records, for each protein, whether it falls in the neighbourhood
of a repeat candidate and which candidate(s) -- that is the join a downstream
domain/function or clustering analysis needs in order to ask whether a protein
family sits near candidates more often than expected.

Nothing here re-predicts or re-annotates: it relabels and indexes frozen output.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ATTR_ID = re.compile(r'ID=([^;]+)')


def gff_coords(path):
    """{cds_id: (contig, start0, end0, strand, partial, conf)} from the CDS GFF."""
    out = {}
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            m = ATTR_ID.search(f[8])
            if not m:
                continue
            attrs = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            out[m.group(1)] = (f[0], int(f[3]) - 1, int(f[4]), f[6],
                               attrs.get('partial', ''), attrs.get('conf', ''))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='append', required=True, metavar='PANEL=DIR')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'proteins').mkdir(exist_ok=True)
    (a.out / 'provenance').mkdir(exist_ok=True)

    rows, totals, combined_parts = [], {}, []
    for spec in a.run:
        panel, _, d = spec.partition('=')
        D = Path(d)
        coords = gff_coords(D / 'cds/genes.gff3')
        # which candidate neighbourhoods each CDS belongs to
        near = defaultdict(list)
        gp = D / 'neighborhoods/genes.tsv'
        if gp.exists():
            with gp.open() as fh:
                for g in csv.DictReader(fh, delimiter='\t'):
                    near[g['gene_id']].append((g['candidate_id'], g['relation'],
                                               g.get('overlap_nt', '0')))
        out_buf = []
        n = 0
        cur_id, buf = None, []

        def flush():
            nonlocal n
            if cur_id is None:
                return
            seq = ''.join(buf)
            new = f'{panel}|{cur_id}'
            c, s0, e0, strand, partial, conf = coords.get(
                cur_id, ('', '', '', '', '', ''))
            hits = near.get(cur_id, [])
            cands = ','.join(h[0] for h in hits)
            rels = ','.join(h[1] for h in hits)
            ov = max([int(h[2] or 0) for h in hits], default=0)
            hdr = (f'>{new} panel={panel} contig={c} start0={s0} end0={e0} '
                   f'strand={strand} len_aa={len(seq)}')
            # Accumulate and write once per panel. Writing each 60-char line
            # straight to Lustre costs ~100 ms per protein here, which turns
            # 122k proteins into hours.
            out_buf.append(hdr + '\n'
                           + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60))
                           + '\n')
            rows.append({
                'protein_id': new, 'panel': panel, 'original_cds_id': cur_id,
                'contig': c, 'start0': s0, 'end0': e0, 'strand': strand,
                'length_aa': len(seq), 'partial': partial, 'pyrodigal_conf': conf,
                'in_candidate_neighborhood': bool(hits),
                'n_candidate_neighborhoods': len(hits),
                'candidate_ids': cands, 'relations': rels,
                'max_overlap_nt_with_candidate': ov,
            })
            n += 1

        with (D / 'cds/proteins.faa').open() as fh:
            for line in fh:
                if line.startswith('>'):
                    flush()
                    cur_id, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
            flush()
        text = ''.join(out_buf)
        (a.out / 'proteins' / f'{panel}.faa').write_text(text)
        combined_parts.append(text)
        totals[panel] = n
        src = D / 'cds/provenance.json'
        if src.exists():
            (a.out / 'provenance' / f'{panel}_cds_provenance.json').write_text(
                src.read_text())
        for extra, dest in [('prepared/provenance.json', 'prepared'),
                            ('calls/provenance.json', 'calls'),
                            ('calls/summary.json', 'calls_summary'),
                            ('neighborhoods/provenance.json', 'neighborhoods')]:
            sp = D / extra
            if sp.exists():
                (a.out / 'provenance' / f'{panel}_{dest}.json').write_text(sp.read_text())
        print(f'  {panel:15s} {n:7,} proteins', flush=True)
    (a.out / 'proteins' / 'all_panels.faa').write_text(''.join(combined_parts))

    with (a.out / 'protein_index.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    ids = {r['protein_id'] for r in rows}
    orig = defaultdict(set)
    for r in rows:
        orig[r['original_cds_id']].add(r['panel'])
    collide = sum(1 for k, v in orig.items() if len(v) > 1)
    meta = {
        'panels': totals, 'total_proteins': len(rows),
        'unique_protein_ids': len(ids),
        'ids_are_unique': len(ids) == len(rows),
        'original_ids_shared_between_panels': collide,
        'why_relabelled':
            'pyrodigal numbers CDSs per run, so every panel restarts at '
            'cds00000001. Pooling the raw proteins.faa files would place '
            f'{collide:,} original IDs in more than one panel and merge unrelated '
            'proteins under one identifier without raising an error. IDs here are '
            '<panel>|<original_id>.',
        'id_format': '<panel>|<original_cds_id>',
        'coordinates': '0-based half-open, on the contig named in the header',
        'annotation_status':
            'pyrodigal ab initio prediction only. There is NO functional '
            'annotation: every product in the source GFF is "hypothetical '
            'protein". Any domain or function assignment must come from the '
            'analysis being run on these sequences, not from this bundle.',
        'candidate_join':
            'in_candidate_neighborhood / candidate_ids come from the frozen '
            'neighbourhood export (candidate +/- flank). A protein being near a '
            'candidate is a positional fact, not evidence of a relationship.',
    }
    (a.out / 'bundle_manifest.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(f"\n  total {len(rows):,} proteins, unique ids: {meta['ids_are_unique']}")
    print(f"  original IDs that appear in >1 panel: {collide:,}")
    print(f"  in a candidate neighbourhood: "
          f"{sum(1 for r in rows if r['in_candidate_neighborhood']):,}")


if __name__ == '__main__':
    main()
