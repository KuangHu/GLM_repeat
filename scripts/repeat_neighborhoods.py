#!/usr/bin/env python3
"""Export candidate-centred gene neighborhoods and optional pyGenomeViz figures.

Overlapping neighborhoods remain separate. No proximity-based biological merge,
operon assignment, protein-family assignment or transcription strand inference.
"""
from __future__ import annotations

import argparse
import html
import json
from collections import defaultdict
from pathlib import Path

from region_io import (SpanIndex, fasta, interval, new_dir, provenance, read_gff,
                       stable_id, table, write_table)
from make_windows import read_fasta


def neighborhood(s,e,length,genes,flank,nearest):
    left = sorted((g for g in genes if g['end0']<=s), key=lambda g:(-g['end0'],g['gene_id']))[:nearest]
    right = sorted((g for g in genes if g['start0']>=e), key=lambda g:(g['start0'],g['gene_id']))[:nearest]
    chosen = left+right+[g for g in genes if SpanIndex(g['parts']).overlap(s,e)]
    start = max(0,min([s-flank]+[g['start0'] for g in chosen]))
    end = min(length,max([e+flank]+[g['end0'] for g in chosen]))
    return start,end


def make_record(key, contig, seq, start,end,cs,ce,genes,proteins,copies,members):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
    record = SeqRecord(Seq(seq[start:end]), id=key, name=key,
                       description=f'{contig}:{start}-{end}; original coordinates 0-based half-open')
    record.annotations.update(molecule_type='DNA', topology='linear')
    record.features.append(SeqFeature(SimpleLocation(cs-start,ce-start),type='misc_feature',
        qualifiers={'label':['candidate'], 'note':['Informative candidate envelope; function unresolved']}))
    for g in genes:
        parts = [(max(s,start),min(e,end)) for s,e in g['parts'] if s<end and e>start]
        if not parts:
            continue
        strand = 1 if g['strand']=='+' else -1
        locations = [SimpleLocation(s-start,e-start,strand=strand) for s,e in parts]
        if strand == -1:
            locations.reverse()
        loc = locations[0] if len(locations)==1 else CompoundLocation(locations)
        clipped = parts != g['parts']
        quals = dict(locus_tag=[g['gene_id']], product=[g['product'] or 'hypothetical protein'],
                     note=[f'original_contig={contig}; original_parts={g["parts"]}; clipped={clipped}'],
                     transl_table=[g['attrs'].get('transl_table','11')])
        if not clipped and g['gene_id'] in proteins:
            quals['translation'] = [proteins[g['gene_id']].rstrip('*')]
        record.features.append(SeqFeature(loc,type='CDS',qualifiers=quals))
    for kind,rows in [('repeat_region',copies),('misc_feature',members)]:
        for r in rows:
            s,e = max(start,int(r['start0'])),min(end,int(r['end0']))
            if s<e:
                record.features.append(SeqFeature(SimpleLocation(s-start,e-start),type=kind,
                    qualifiers={'label':[r.get('copy_id','detector_cluster')],
                                'note':[r.get('model_id',r.get('window_id',''))]}))
    return record


def draw(path,c,s,e,cs,ce,genes,copies,members):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pygenomeviz import GenomeViz
    gv = GenomeViz(fig_width=14,fig_track_height=1.0)
    gv.set_scale_xticks()
    cds = gv.add_feature_track('CDS', (s,e))
    cds.add_sublabel(f'{c}:{s}-{e} [0-based); axis relative to region start; gene IDs in genes.tsv')
    for g in genes:
        for i,(gs,ge) in enumerate(g['parts']):
            gs,ge = max(s,gs),min(e,ge)
            if gs<ge:
                cds.add_feature(gs,ge,1 if g['strand']=='+' else -1,fc='#6486a5',
                                label=g['gene_id'] if i==0 else '',text_kws={'size':6})
    cand = gv.add_feature_track('Candidate', (s,e))
    cand.add_feature(cs,ce,plotstyle='box',fc='#e4a23b')
    if members:
        mt = gv.add_feature_track('Source clusters',(s,e))
        for r in members:
            x,y = max(s,int(r['start0'])),min(e,int(r['end0']))
            if x<y:
                mt.add_feature(x,y,plotstyle='box',fc='#b8a1d9',alpha=0.3)
    if copies:
        rt = gv.add_feature_track('Repeat hypotheses',(s,e))
        for x,y in sorted({(max(s,int(r['start0'])),min(e,int(r['end0']))) for r in copies}):
            if x<y:
                rt.add_feature(x,y,plotstyle='box',fc='#207b69',alpha=0.7)
    fig = gv.plotfig()
    fig.savefig(path, bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'),dpi=110,bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--gff',type=Path,required=True)
    p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--proteins',type=Path,help='Protein FASTA IDs must match GFF CDS ID')
    p.add_argument('--copies',type=Path,help='Explainer v0.3 copies.tsv')
    p.add_argument('--models',type=Path,help='Only display copy model IDs listed here (e.g. primary models.tsv)')
    p.add_argument('--members',type=Path,help='call_region_candidates members.tsv')
    p.add_argument('--flank',type=int,default=10000)
    p.add_argument('--nearest-genes',type=int,default=5)
    p.add_argument('--plot',action='store_true',help='Requires pygenomeviz; SVG + PNG per candidate')
    a = p.parse_args()
    if min(a.flank,a.nearest_genes)<0 or (a.models and not a.copies):
        p.error('Lengths/counts must be >=0; --models requires --copies')
    seqs = fasta(a.fasta)
    genes = read_gff(a.gff,seqs)
    by = defaultdict(list)
    for g in genes:
        by[g['contig']].append(g)
    proteins = dict(read_fasta(a.proteins)) if a.proteins else {}
    if a.proteins:
        ids = [g['gene_id'] for g in genes]
        if len(ids)!=len(set(ids)):
            raise ValueError('Protein association requires globally unique CDS IDs')
        if not set(ids)&set(proteins) and ids:
            raise ValueError('No protein FASTA IDs match GFF CDS IDs')
    candidates = table(a.candidates)
    seen = set()
    for r in candidates:
        c,s,e = interval(r,seqs)
        # The v0.3 explainer also derives IDs from these coordinates.
        cid = f'{c}:{s}-{e}'
        if r.get('candidate_id',cid) != cid or cid in seen:
            raise ValueError(f'Candidate ID must be unique and coordinate-derived: {cid}')
        r['candidate_id'] = cid
        seen.add(cid)
    model_ids = {r['model_id'] for r in table(a.models)} if a.models else None
    copy_by, member_by = defaultdict(list),defaultdict(list)
    for path,dest in [(a.copies,copy_by),(a.members,member_by)]:
        if path:
            for r in table(path):
                interval(r,seqs)
                if r['candidate_id'] not in seen:
                    raise ValueError(f'{path}: unknown candidate {r["candidate_id"]}')
                if path==a.copies and model_ids is not None and r['model_id'] not in model_ids:
                    continue
                dest[r['candidate_id']].append(r)
    from Bio import SeqIO
    out = new_dir(a.out)
    (out/'genbank').mkdir()
    if a.plot:
        (out/'plots').mkdir()
    loci,neighbors,items = [],[],[]
    with (out/'regions.fa').open('w',encoding='utf-8') as fh:
        for r in candidates:
            c,cs,ce = interval(r,seqs)
            cid,key = r['candidate_id'],stable_id('n',r['candidate_id'])
            s,e = neighborhood(cs,ce,len(seqs[c]),by[c],a.flank,a.nearest_genes)
            near = [g for g in by[c] if SpanIndex(g['parts']).overlap(s,e)]
            copies,members = copy_by[cid],member_by[cid]
            if any(x['contig']!=c for x in copies+members):
                raise ValueError(f'Evidence contig mismatch: {cid}')
            record = make_record(key,c,seqs[c],s,e,cs,ce,near,proteins,copies,members)
            SeqIO.write(record,out/'genbank'/f'{key}.gbk','genbank')
            fh.write(f'>{key} {c}:{s}-{e}\n{seqs[c][s:e]}\n')
            full_translations = sum('translation' in f.qualifiers for f in record.features if f.type=='CDS')
            loci.append(dict(candidate_id=cid,neighborhood_id=key,contig=c,start0=s,end0=e,
                candidate_start0=cs,candidate_end0=ce,n_CDS=len(near),n_translations=full_translations,
                edge_truncated=s==0 or e==len(seqs[c]),contains_ambiguous=any(x not in 'ACGT' for x in seqs[c][s:e]),
                genbank=f'genbank/{key}.gbk',plot=f'plots/{key}.svg' if a.plot else ''))
            for g in near:
                ov = SpanIndex(g['parts']).overlap(cs,ce)
                relation = 'overlap' if ov else 'left' if g['end0']<=cs else 'right' if g['start0']>=ce else 'between_CDS_parts'
                distance = cs-g['end0'] if relation=='left' else g['start0']-ce if relation=='right' else 0
                neighbors.append(dict(candidate_id=cid,gene_id=g['gene_id'],contig=c,
                    start0=g['start0'],end0=g['end0'],strand=g['strand'],parts=json.dumps(g['parts']),
                    relation=relation,distance_nt=distance,overlap_nt=ov,product=g['product'],
                    protein_available=g['gene_id'] in proteins,attributes=json.dumps(g['attrs'],ensure_ascii=False)))
            if a.plot:
                draw(out/'plots'/f'{key}.svg',c,s,e,cs,ce,near,copies,members)
            items.append(f'<li>{html.escape(cid)} — <a href="genbank/{key}.gbk">GenBank</a>'+
                         (f' | <a href="plots/{key}.svg">figure</a>' if a.plot else '')+'</li>')
    write_table(out/'neighborhoods.tsv',loci,['candidate_id','neighborhood_id','contig','start0','end0',
        'candidate_start0','candidate_end0','n_CDS','n_translations','edge_truncated','contains_ambiguous','genbank','plot'])
    write_table(out/'genes.tsv',neighbors,['candidate_id','gene_id','contig','start0','end0','strand','parts',
        'relation','distance_nt','overlap_nt','product','protein_available','attributes'])
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Repeat neighborhoods</title>'
        '<h1>Repeat neighborhoods</h1><p>Candidate envelopes and repeat hypotheses; gene proximity does not '
        'establish function or co-transcription. Overlapping neighborhoods remain separate.</p><ul>'+
        '\n'.join(items)+'</ul>',encoding='utf-8')
    provenance(out,a,[a.fasta,a.gff,a.candidates,a.proteins,a.copies,a.models,a.members],
        n_neighborhoods=len(loci),limitations=['no circular-origin wrap-around',
            'protein families and synteny links require external homology analysis (e.g. clinker)',
            'plot shows exact-interval deduplication only; original model hypotheses are retained upstream',
            'gene left/right is genomic direction, not inferred repeat transcription direction'])
    print(f'{len(loci)} neighborhoods, {len(neighbors)} gene associations -> {out}')


if __name__ == '__main__':
    main()
