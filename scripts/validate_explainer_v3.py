#!/usr/bin/env python3
"""Check saved v0.3 evidence against original DNA and annotation (CPU only).

Annotations are evaluation inputs only. Validation here means engineering
invariants; annotated-unit overlap is not a statement of biological validation.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
from explain_repeats_v3 import read_fasta, read_table, write_table, mismatches, find_copies, overlap, union_length


def check_outputs(folder, genomes):
    hypotheses = read_table(folder/'hypotheses.tsv')
    models = {m['model_id']:m for m in hypotheses}
    copies = read_table(folder/'copies.tsv')
    by_model, by_copy = defaultdict(list), {}
    failures = []
    def check(ok, detail):
        if not ok:
            failures.append(detail)
    for c in copies:
        check(c['copy_id'] not in by_copy,'duplicate copy id: '+c['copy_id'])
        by_copy[c['copy_id']]=c
        by_model[c['model_id']].append(c)
        m=models[c['model_id']]
        a,b,cs,ce=map(int,[c['start0'],c['end0'],m['cand_start0'],m['cand_end0']])
        check(genomes[c['contig']][a:b]==c['sequence'],'coordinate mismatch: '+c['copy_id'])
        mm=mismatches(c['sequence'],m['consensus'])
        check(mm==int(c['mismatches_to_consensus']),'consensus mismatch: '+c['copy_id'])
        check(abs(float(c['consensus_identity'])-(1-mm/(b-a)))<1e-10,'identity mismatch: '+c['copy_id'])
        check((c['overlaps_candidate']=='True')==(overlap(a,b,cs,ce)>0),'overlap mismatch: '+c['copy_id'])
        check((c['contained_in_candidate']=='True')==(cs<=a and b<=ce),'containment mismatch: '+c['copy_id'])
        if c['search_motif']:
            check(int(c['search_mismatches'])==mismatches(c['sequence'],c['search_motif']),'search distance mismatch: '+c['copy_id'])
            check(int(c['search_mismatches'])<=int(c['search_budget']),'unsupported core: '+c['copy_id'])
    for mid,m in models.items():
        cc=by_model[mid]
        check(len(cc)==int(m['n_copies']),'model count mismatch: '+mid)
        n=sum(c['overlaps_candidate']=='True' for c in cc)
        check(n==int(m['n_copies_overlapping_candidate']),'overlap count mismatch: '+mid)
        check((m['association']=='candidate_associated')==(n>0),'association mismatch: '+mid)
        cs,ce=int(m['cand_start0']),int(m['cand_end0'])
        covered=union_length([(max(cs,int(c['start0'])),min(ce,int(c['end0']))) for c in cc if c['overlaps_candidate']=='True'])
        check(covered==int(m['candidate_coverage_nt']),'coverage mismatch: '+mid)
    summary=json.loads((folder/'summary.json').read_text())
    flank=summary['parameters']['flank']
    # For these unambiguous source regions, reconstruct each saved search-motif
    # core set independently of extensions and verify every coordinate is emitted.
    core_models_checked=0
    for mid,m in models.items():
        if m['model_source']!='motif_core':
            continue
        cc=by_model[mid]; g=genomes[m['contig']]
        lo=max(0,int(m['cand_start0'])-flank); hi=min(len(g),int(m['cand_end0'])+flank)
        region=g[lo:hi]
        if set(region)-set('ACGT'):
            continue  # Ambiguous segments are covered by the synthetic regression.
        expected={(a+lo,b+lo) for a,b in find_copies(region,cc[0]['search_motif'],int(cc[0]['search_budget']))}
        actual={(int(c['start0']),int(c['end0'])) for c in cc}
        check(actual==expected,'core set lost/recruited: '+mid)
        core_models_checked+=1
    attempts=read_table(folder/'extensions.tsv')
    accepted=Counter()
    rejected=0
    for x in attempts:
        c=by_copy[x['core_copy_id']]; a,b=int(x['start0']),int(x['end0'])
        check(x['model_id']==c['model_id'],'wrong extension link: '+x['core_copy_id'])
        check(a<=int(c['start0'])<int(c['end0'])<=b,'extension does not contain core: '+x['core_copy_id'])
        check(genomes[c['contig']][a:b]==x['sequence'],'extension coordinate mismatch: '+x['core_copy_id'])
        mm=mismatches(x['sequence'],x['search_motif'])
        check(mm==int(x['search_mismatches']),'extension search mismatch: '+x['core_copy_id'])
        check((x['accepted']=='True')==(mm<=int(x['search_budget'])),'extension decision mismatch: '+x['core_copy_id'])
        if x['final_consensus']:
            check(mismatches(x['sequence'],x['final_consensus'])==int(x['mismatches_to_consensus']),
                  'extension final distance mismatch: '+x['core_copy_id'])
        if x['accepted']=='True':
            accepted[x['model_id']]+=1
        else:
            rejected+=1
    for mid,m in models.items():
        check(accepted[mid]==int(m['extended_copies']),'extension count mismatch: '+mid)
    cand=read_table(folder/'candidates.tsv')
    for c in cand:
        mm=[m for m in hypotheses if m['candidate_id']==c['candidate_id']]
        status='candidate_associated' if any(m['association']=='candidate_associated' for m in mm) else 'flank_only' if mm else 'unresolved'
        check(status==c['status'],'candidate status mismatch: '+c['candidate_id'])
        check(len(mm)==int(c['n_models']),'candidate hypothesis count mismatch: '+c['candidate_id'])
    return dict(n_candidates=len(cand),n_hypotheses=len(models),n_copies=len(copies),
                n_core_models_reconstructed=core_models_checked,n_extension_attempts=len(attempts),
                rejected_extensions_with_core_retained=rejected,failures=failures)


def control_overlap(folder, definitions, units):
    models=read_table(folder/'models.tsv')
    copies=defaultdict(list)
    for c in read_table(folder/'copies.tsv'):
        copies[c['model_id']].append(c)
    definitions={f"{d['contig']}:{d['start0']}-{d['end0']}":d for d in definitions}
    rows=[]
    for m in models:
        d=definitions[m['candidate_id']]
        gold=[u for u in units if u['locus_id']==d['locus']]
        cc=copies[m['model_id']]
        hit=lambda c,u: overlap(int(c['start0']),int(c['end0']),int(u['start0']),int(u['end0']))>0
        rows.append(dict(locus=d['locus'],candidate_id=m['candidate_id'],primary_source=m['model_source'],
            core_or_periodic_copies=len(cc),unit_len=m['unit_len'],extended_copies=m['extended_copies'],
            extension_unit_len=m['extension_unit_len'],period=m['period'],
            annotated_units=len(gold),annotated_units_overlapped=sum(any(hit(c,u) for c in cc) for u in gold) if gold else '',
            extra_copies=sum(not any(hit(c,u) for u in gold) for c in cc) if gold else '',
            evaluation_scope='interval_overlap_only' if gold else 'no_copy_level_annotation'))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,required=True)
    p.add_argument('--heldout',type=Path,required=True)
    p.add_argument('--genomes',type=Path,required=True)
    p.add_argument('--control-definitions',type=Path,required=True)
    p.add_argument('--annotated-units',type=Path,required=True)
    p.add_argument('--v2-heldout',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    genomes=read_fasta(a.genomes)
    results={name:check_outputs(folder,genomes) for name,folder in [('controls',a.controls),('heldout',a.heldout)]}
    overlaps=control_overlap(a.controls,read_table(a.control_definitions),read_table(a.annotated_units))
    write_table(a.out/'controls_overlap.tsv',overlaps,list(overlaps[0]))
    old={m['candidate_id']:m for m in read_table(a.v2_heldout/'models.tsv')}
    new={m['candidate_id']:m for m in read_table(a.heldout/'models.tsv')}
    changes=[]
    for c in read_table(a.heldout/'candidates.tsv'):
        o=old.get(c['candidate_id'],{}); n=new.get(c['candidate_id'],{})
        changes.append(dict(candidate_id=c['candidate_id'],v2_status='explained' if o else 'unresolved',
            v3_status=c['status'],v2_copies=o.get('n_copies',''),v3_primary_copies=n.get('n_copies',''),
            v2_source=o.get('model_source',''),v3_source=n.get('model_source',''),
            v3_associated_hypotheses=c['n_candidate_associated_models'],
            comparison_note='v3 primary motif copies are cores; extensions and alternative hypotheses are separate'))
    write_table(a.out/'candidate_changes.tsv',changes,list(changes[0]))
    results['status_transitions']=dict(Counter(f"{r['v2_status']} -> {r['v3_status']}" for r in changes))
    results['passed']=not any(r['failures'] for key,r in results.items() if key in ('controls','heldout'))
    (a.out/'validation.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2))
    if not results['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
