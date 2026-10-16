"""Read-only audit of frozen outputs. Writes evidence under independent_review only."""
import csv
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path
import random
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
PKT = ROOT / 'review_packet_explainer_v0.2'
spec = importlib.util.spec_from_file_location('ex', ROOT / 'scripts/explain_repeats_v2.py')
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)

def table(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def emit(name, rows):
    if not rows:
        return
    with (OUT / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def regions(path):
    result = {}
    for header, seq in ex.read_fasta(path).items():
        cid, region, *_ = header.split('|')
        s, e = map(int, region.removeprefix('region=').split('-'))
        assert len(seq) == e-s
        result[cid] = (s, e, seq)
    return result

def by_candidate(rows):
    out = defaultdict(list)
    for r in rows:
        out[r['candidate_id']].append(r)
    return out

summary = {}
genome_path = OUT/'source/genomes.fa'
genomes = ex.read_fasta(genome_path) if genome_path.exists() else {}
allregions = {}
for kind, fn in [('controls', 'controls_regions.fa'), ('heldout_primary', 'review_set_regions.fa')]:
    reg = regions(PKT / 'data/regions' / fn)
    allregions.update(reg)
    models = table(PKT / 'outputs' / kind / 'models.tsv')
    copies = table(PKT / 'outputs' / kind / 'copies.tsv')
    modelmap = {m['candidate_id']: m for m in models}
    cp = by_candidate(copies)
    wrong_mm, wrong_seq = [], []
    checked_seq = 0
    for c in copies:
        cid = c['candidate_id']
        actual = ex.mismatches(c['sequence'], modelmap[cid]['consensus'])
        if actual != int(c['mismatches_to_consensus']):
            wrong_mm.append({**c, 'actual_mismatches_to_consensus': actual,
                             'actual_support_to_consensus': round(1-actual/len(c['sequence']), 3)})
        if c['contig'] in genomes:
            checked_seq += 1
            if genomes[c['contig']][int(c['start0']):int(c['end0'])] != c['sequence']:
                wrong_seq.append(c)
        elif cid in reg:
            s,e,seq = reg[cid]
            checked_seq += 1
            if seq[int(c['start0'])-s:int(c['end0'])-s] != c['sequence']:
                wrong_seq.append(c)
    emit(f'{kind}_incorrect_mismatch_fields.tsv', wrong_mm)
    emit(f'{kind}_incorrect_coordinate_sequences.tsv', wrong_seq)
    outside = [m for m in models if int(m['n_copies_inside_candidate']) == 0]
    emit(f'{kind}_zero_copies_inside.tsv', outside)
    no_overlap = [m for m in models if not any(int(c['start0'])<int(m['cand_end0']) and int(m['cand_start0'])<int(c['end0']) for c in cp[m['candidate_id']])]
    emit(f'{kind}_zero_copies_overlapping.tsv',no_overlap)
    summary[kind] = dict(models=len(models), copies=len(copies),
        checked_coordinate_sequences=checked_seq, sequence_mismatches=len(wrong_seq),
        incorrect_mismatch_fields=len(wrong_mm), zero_copies_inside=len(outside),
        zero_copies_overlapping=len(no_overlap),
        low_copy_models=sum(int(m['n_copies'])<=3 for m in models),
        model_sources=dict(Counter(m['model_source'] for m in models)))
    if kind == 'controls':
        units = table(PKT/'data/annotated_units.tsv')
        defs = table(PKT/'data/controls_definition.tsv')
        metrics = []
        for d in defs:
            cid = f"{d['contig']}:{d['start0']}-{d['end0']}"
            gold = [u for u in units if u['locus_id']==d['locus']]
            matched = [u for u in gold if any(int(c['start0'])<int(u['end0']) and int(u['start0'])<int(c['end0']) for c in cp[cid])]
            extra = [c for c in cp[cid] if not any(int(c['start0'])<int(u['end0']) and int(u['start0'])<int(c['end0']) for u in gold)]
            metrics.append(dict(locus=d['locus'], candidate_id=cid, source=modelmap[cid]['model_source'],
                unit_len=modelmap[cid]['unit_len'], copies=len(cp[cid]), annotated_units=len(gold),
                annotated_units_overlapped=len(matched), copies_not_overlapping_annotation=len(extra),
                extra_starts=','.join(c['start0'] for c in extra),
                missed_starts=','.join(u['start0'] for u in gold if u not in matched)))
        emit('controls_annotation_overlap.tsv', metrics)
        summary['controls']['annotation_overlap']=metrics
        # Recompute frozen models on the EXACT saved input, without changing any rule.
        reproduction = []
        for m in models:
            cid=m['candidate_id']; s,e,seq=reg[cid]
            span=seq[int(m['cand_start0'])-s:int(m['cand_end0'])-s]
            found=ex.motif_models(seq)
            d,ident=ex.best_period(span)
            per=None
            if d and ident>=ex.PERIODIC_MIN_IDENTITY:
                iv=ex.periodic_interval(seq,d)
                if iv: per=ex.periodic_model(seq,*iv,d)
            pool=[v for v in [found[0] if found else None,per] if v]
            pool.sort(key=lambda v: -len(v['copies'])*v['unit_len']*v['information_bits'])
            v=pool[0]
            reproduction.append(dict(candidate_id=cid, copies_equal=len(v['copies'])==int(m['n_copies']),
                consensus_equal=v['motif']==m['consensus'], source_equal=v['source']==m['model_source']))
        emit('controls_reproduction.tsv',reproduction)
        summary['controls']['reproduction_all_equal']=all(all(v for k,v in r.items() if k!='candidate_id') for r in reproduction)
        vipr = next(m for m in models if m['annotation_category']=='control_VIPR_labelled_60nt')
        vipr_rows=[]
        for c in cp[vipr['candidate_id']]:
            others=[x['sequence'] for x in cp[vipr['candidate_id']] if x is not c]
            loo_cons=ex.consensus(ex.pwm_from(others))
            start,end=int(c['start0']),int(c['end0'])
            relation='contained' if 29169<=start and end<=29229 else ('boundary_crossing' if start<29229 and 29169<end else 'outside')
            vipr_rows.append({**c,'label_relation':relation,'loo_consensus':loo_cons,
                'loo_hamming_support':round(1-ex.mismatches(c['sequence'],loo_cons)/5,3),
                'rna_orientation_unit':c['sequence'].translate(str.maketrans('ACGT','TGCA'))[::-1]})
        emit('VIPR_units.tsv',vipr_rows)
        summary['VIPR_label_relation']=dict(Counter(c['label_relation'] for c in vipr_rows))

models=table(PKT/'outputs/heldout_primary/models.tsv')
copies=table(PKT/'outputs/heldout_primary/copies.tsv')
cp=by_candidate(copies)
lost_core_models=[m for m in models if m['model_source']=='motif_extended' and int(m['n_rejected_cores'])>0]
emit('extended_models_with_rejected_cores.tsv',lost_core_models)
summary['extended_models_with_rejected_cores']=len(lost_core_models)
summary['rejected_cores_in_selected_extended_models']=sum(int(m['n_rejected_cores']) for m in lost_core_models)
if genomes:
    core_rechecks=[]
    for m in lost_core_models:
        cs,ce=int(m['cand_start0']),int(m['cand_end0'])
        g=genomes[m['contig']]; s,e=max(0,cs-500),min(len(g),ce+500)
        seq=g[s:e]; v=ex.motif_models(seq)[0]
        core=ex.find_copies(seq,v['core_consensus'],2)
        absent=[(a,b) for a,b,mm in core if not any(x<=a and b<=y for x,y,z in v['copies'])]
        core_rechecks.append(dict(candidate_id=m['candidate_id'],core_copies=len(core),
            emitted_extended_copies=len(v['copies']),core_not_contained_in_any_emitted_copy=len(absent),
            absent_core_starts=','.join(str(s+a) for a,b in absent)))
    emit('lost_core_coordinates.tsv',core_rechecks)
    summary['selected_extended_models_with_physically_missing_cores']=sum(r['core_not_contained_in_any_emitted_copy']>0 for r in core_rechecks)
    summary['core_copies_not_contained_in_selected_extended_copies']=sum(r['core_not_contained_in_any_emitted_copy'] for r in core_rechecks)
random10=random.Random(20261015).sample(sorted(models,key=lambda x:x['candidate_id']),10)
emit('random10_fixed.tsv',random10)
emit('random10_copies.tsv',[c for m in random10 for c in cp[m['candidate_id']]])
summary['random10_seed']=20261015
summary['random10_original_regions_available']=sum(m['candidate_id'] in allregions for m in random10)
if genomes:
    to_inspect=table(PKT/'outputs/review_set/review_set.tsv')
    wanted=list(dict.fromkeys([r['representative_candidate'] for r in to_inspect]+[m['candidate_id'] for m in random10]))
    with (OUT/'review33_regions.fa').open('w',encoding='utf-8') as f:
        for cid in wanted:
            contig,interval=cid.split(':'); cs,ce=map(int,interval.split('-'))
            seq=genomes[contig]; s,e=max(0,cs-500),min(len(seq),ce+500)
            f.write(f'>{cid}|region={s}-{e}\n{seq[s:e]}\n')
    summary['genome_source']=dict(records=len(genomes),sha256=hashlib.sha256(genome_path.read_bytes()).hexdigest())
    summary['packet_regions_matching_NCBI']=sum(genomes[cid.split(':')[0]][s:e]==seq for cid,(s,e,seq) in allregions.items())
    summary['packet_regions_checked']=len(allregions)
    gbpath=OUT/'source/genomes.gb'
    if gbpath.exists():
        from Bio import SeqIO
        records={r.id:r for r in SeqIO.parse(str(gbpath),'genbank')}
        rows=[]
        for cid in wanted:
            contig,interval=cid.split(':'); cs,ce=map(int,interval.split('-'))
            cds=[ft for ft in records[contig].features if ft.type=='CDS' and any(int(p.start)<ce and cs<int(p.end) for p in ft.location.parts)]
            rows.append(dict(candidate_id=cid,n_CDS_overlapping=len(cds),
                products=';'.join(','.join(ft.qualifiers.get('product',[])) for ft in cds),
                cds_coordinates=';'.join(str(ft.location) for ft in cds),
                source='fresh NCBI GenBank; sequence matches frozen region'))
        emit('review_candidates_CDS.tsv',rows)
        summary['fresh_GenBank_sequence_equals_FASTA']=all(str(records[k].seq).upper()==v for k,v in genomes.items())
    # Compact casebook: every preselected representative, every unresolved case,
    # and the fixed random sample, retaining the raw alignments for human review.
    modelmap={m['candidate_id']:m for m in models}
    with (OUT/'casebook.txt').open('w',encoding='utf-8') as f:
        for cid in wanted:
            f.write('\n'+cid+'\n')
            if cid not in modelmap:
                contig,interval=cid.split(':'); cs,ce=map(int,interval.split('-'))
                seq=genomes[contig][cs:ce]
                base=Counter(seq)
                f.write('UNRESOLVED base composition '+str(base)+'\n')
                for i in range(0,len(seq),100): f.write(seq[i:i+100]+'\n')
                continue
            m=modelmap[cid]
            f.write(' '.join(f'{k}={m[k]}' for k in ['model_source','architecture','unit_len','n_copies','n_copies_inside_candidate','core_copies','extended_copies','n_rejected_cores'])+'\n')
            f.write('consensus '+m['consensus']+'\n')
            for c in cp[cid]:
                f.write(f"{c['start0']}-{c['end0']} inside={c['inside_candidate_span']} mm_actual={ex.mismatches(c['sequence'],m['consensus'])} {c['sequence']}\n")

manifest_bad=[]; manifest_count=0; manifest_git_bad=[]
for line in (PKT/'MANIFEST.sha256').read_text().splitlines():
    h,p=line.split(maxsplit=1); f=PKT/p.lstrip('*')
    actual=hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else 'MISSING'
    manifest_count+=1
    if actual!=h: manifest_bad.append(dict(path=p,expected=h,actual=actual))
    gitpath='review_packet_explainer_v0.2/'+p.removeprefix('./')
    blob=subprocess.run(['git','show','HEAD:'+gitpath],cwd=ROOT,capture_output=True)
    gitactual=hashlib.sha256(blob.stdout).hexdigest() if blob.returncode==0 else 'MISSING_IN_GIT'
    if gitactual!=h: manifest_git_bad.append(dict(path=p,expected=h,actual=gitactual))
emit('manifest_mismatches.tsv',manifest_bad)
emit('manifest_git_mismatches.tsv',manifest_git_bad)
summary['manifest']=dict(checked=manifest_count,working_tree_mismatches=len(manifest_bad),git_blob_mismatches=len(manifest_git_bad))
summary['code_packet_equal']=(ROOT/'scripts/explain_repeats_v2.py').read_bytes()==(PKT/'code/explain_repeats_v2.py').read_bytes()
reproduced=OUT/'reproduced_heldout/models.tsv'
if reproduced.exists():
    saved=table(PKT/'outputs/heldout_primary/models.tsv'); new=table(reproduced)
    ignore={'copies_with_>=1_peak'}
    discrepancies=[dict(candidate_id=a['candidate_id'],field=k,old=a[k],new=b[k]) for a,b in zip(saved,new) for k in a if k not in ignore and a[k]!=b[k]]
    emit('reproduction_differences.tsv',discrepancies)
    summary['heldout_reproduction']=dict(saved=len(saved),reproduced=len(new),model_field_differences=len(discrepancies),ignored_fields=sorted(ignore))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
