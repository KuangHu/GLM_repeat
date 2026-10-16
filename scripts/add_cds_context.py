#!/usr/bin/env python3
"""Append CDS context without filtering, ranking or changing detector coordinates.

Existing annotation_category is retained as a selected RNA/repeat category, not
as a claim that other GenBank annotations are absent. GenBank inputs are explicit.
"""
import argparse
import hashlib
import json
from pathlib import Path

from explain_repeats_v3 import read_table, write_table, union_length


def cds_context(record, start, end):
    if record is None:
        return dict(annotation_available=False,n_CDS_overlapping='',CDS_overlap_nt='',
                    CDS_overlap_fraction='',CDS_products='',CDS_locations='')
    hits, spans = [], []
    for feature in record.features:
        if feature.type != 'CDS' or feature.location is None:
            continue
        clipped = [(max(start,int(p.start)),min(end,int(p.end))) for p in feature.location.parts
                   if int(p.start)<end and start<int(p.end)]
        if clipped:
            hits.append(feature)
            spans.extend(clipped)
    nt = union_length(spans)
    return dict(annotation_available=True,n_CDS_overlapping=len(hits),CDS_overlap_nt=nt,
                CDS_overlap_fraction=nt/(end-start),
                CDS_products=';'.join('|'.join(f.qualifiers.get('product',[])) for f in hits),
                CDS_locations=';'.join(str(f.location) for f in hits))


def main():
    from Bio import SeqIO
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--genbank',type=Path,nargs='+',required=True)
    p.add_argument('--out',type=Path,required=True)
    a = p.parse_args()
    records = {}
    for path in a.genbank:
        for r in SeqIO.parse(str(path),'genbank'):
            if r.id in records:
                raise ValueError(f'Duplicate GenBank record: {r.id}')
            records[r.id] = r
    original = read_table(a.candidates)
    rows = [{**c,**cds_context(records.get(c['contig']),int(c['start0']),int(c['end0'])),
             'annotation_category_scope':'selected_RNA_repeat_categories_only',
             'CDS_annotation_source':';'.join(str(p) for p in a.genbank)} for c in original]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    fields = list(rows[0]) if rows else ['contig','start0','end0',*cds_context(None,0,1),
                                       'annotation_category_scope','CDS_annotation_source']
    write_table(a.out,rows,fields)
    sources = [a.candidates,*a.genbank]
    a.out.with_suffix('.provenance.json').write_text(json.dumps({
        'sources_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'policy':'annotation is descriptive; all candidates and original categories retained',
        'category_scope':'unclassified means no selected RNA/repeat category, not no annotation',
        'n_candidates':len(rows)},indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
