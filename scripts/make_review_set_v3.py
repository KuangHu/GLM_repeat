#!/usr/bin/env python3
"""Review v0.3 hypotheses with explicit statuses and CDS context.

Exact strand-canonical consensus groups are navigation aids, not biological
families. Periodic consensus keys also allow cyclic column rotations.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from explain_repeats_v3 import read_table, write_table, read_fasta


def pattern_key(m):
    seq=m['consensus']; rc=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
    variants=[seq,rc]
    if m['model_source']=='periodic':
        variants=[s[i:]+s[:i] for s in variants for i in range(len(s))]
    return m['model_source']+':'+min(variants)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--genomes',type=Path,required=True)
    p.add_argument('--previous-review',type=Path)
    p.add_argument('--fixed-random',type=Path)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    models={m['candidate_id']:m for m in read_table(a.results/'models.tsv')}
    candidates={c['candidate_id']:c for c in read_table(a.results/'candidates.tsv')}
    groups=defaultdict(list)
    for m in models.values():
        if m['association']=='candidate_associated':
            groups[pattern_key(m)].append(m)
    group_rows=[]
    for i,(key,ms) in enumerate(sorted(groups.items()),1):
        ms.sort(key=lambda m:(-float(m['detector_score']),m['candidate_id']))
        m=ms[0]
        group_rows.append(dict(pattern_group_id=f'g{i:04d}',pattern_key=key,n_members=len(ms),
            representative_candidate=m['candidate_id'],architecture=m['architecture'],
            detector_score=m['detector_score'],member_candidates=';'.join(x['candidate_id'] for x in ms)))
    write_table(a.out/'pattern_groups.tsv',group_rows,['pattern_group_id','pattern_key','n_members',
                'representative_candidate','architecture','detector_score','member_candidates'])
    selected=defaultdict(list); per_arch=Counter()
    for group in sorted(group_rows,key=lambda r:-float(r['detector_score'])):
        m=models[group['representative_candidate']]
        if m['annotation_category']=='unclassified' and per_arch[m['architecture']]<3:
            per_arch[m['architecture']]+=1
            selected[m['candidate_id']].append('current_representative_no_selected_RNA_repeat_category')
    for status in ('unresolved','flank_only'):
        for c in sorted((c for c in candidates.values() if c['status']==status),
                        key=lambda c:(-float(c['detector_score']),c['candidate_id']))[:8]:
            selected[c['candidate_id']].append('high_scoring_'+status)
    if a.previous_review:
        for r in read_table(a.previous_review):
            selected[r['representative_candidate']].append('fixed_v02_review')
    if a.fixed_random:
        for r in read_table(a.fixed_random):
            selected[r['candidate_id']].append('fixed_random_seed_20261015')
    rows=[]
    for cid,reasons in selected.items():
        c=candidates[cid]; m=models.get(cid,{})
        rows.append(dict(candidate_id=cid,status=c['status'],review_reason=';'.join(reasons),
            primary_model_id=c['primary_model_id'],n_candidate_associated_models=c['n_candidate_associated_models'],
            architecture=m.get('architecture',''),consensus=m.get('consensus',''),unit_len=m.get('unit_len',''),
            n_primary_copies=m.get('n_copies',''),detector_score=c['detector_score'],
            annotation_category=c['annotation_category'],n_CDS_overlapping=c['n_CDS_overlapping'],
            CDS_overlap_fraction=c['CDS_overlap_fraction'],annotation_available=c['annotation_available']))
    fields=['candidate_id','status','review_reason','primary_model_id','n_candidate_associated_models','architecture',
            'consensus','unit_len','n_primary_copies','detector_score','annotation_category','n_CDS_overlapping',
            'CDS_overlap_fraction','annotation_available']
    write_table(a.out/'review_set.tsv',rows,fields)
    copies=read_table(a.results/'copies.tsv')
    ids={r['primary_model_id'] for r in rows if r['primary_model_id']}
    subset=[c for c in copies if c['model_id'] in ids]
    write_table(a.out/'primary_copies.tsv',subset,list(copies[0]) if copies else ['copy_id'])
    genomes=read_fasta(a.genomes)
    with (a.out/'regions.fa').open('w',encoding='utf-8') as f:
        for cid in selected:
            c=candidates[cid]; g=genomes[c['contig']]
            s=max(0,int(c['start0'])-500); e=min(len(g),int(c['end0'])+500)
            f.write(f'>{cid}|region={s}-{e}\n{g[s:e]}\n')
    summary=dict(n_groups=len(groups),n_unique_review_candidates=len(rows),
        current_representatives_by_architecture=dict(per_arch),
        grouping='exact strand-canonical consensus; cyclic rotations for periodic profiles; not biological systems',
        selection='current representatives + unresolved/flank-only + fixed previous review + unchanged random10')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
