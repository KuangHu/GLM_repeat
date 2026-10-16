"""Behavioral regression tests for selection, audit, provenance and neighborhoods."""
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
from region_io import SpanIndex, read_gff, table
from prepare_region_scan import regions_for, choose_audit
from call_region_candidates import validate_tracks, decorate
from repeat_neighborhoods import neighborhood, make_record


def gene(s,e,conf=100,partial='00',gid='g1',parts=None,strand='+'):
    parts = parts or [(s,e)]
    return dict(contig='c',gene_id=gid,start0=s,end0=e,strand=strand,parts=parts,
                length_nt=sum(b-a for a,b in parts),product='test',source='test',
                attrs=dict(conf=str(conf),cscore='100',partial=partial))


def cli(script,*args):
    return subprocess.run([sys.executable,str(SCRIPTS/script),*map(str,args)],
                          check=True,capture_output=True,text=True)


def test_interval_union_index_matches_basewise_oracle():
    rng = random.Random(7)
    spans = [(a,a+rng.randrange(1,30)) for a in [rng.randrange(200) for _ in range(40)]]
    covered = {x for s,e in spans for x in range(s,e)}
    idx = SpanIndex(spans)
    for s in range(0,250,3):
        for length in (1,7,63):
            assert idx.overlap(s,s+length)==len(covered & set(range(s,s+length)))


def test_weak_orf_does_not_split_long_noncoding():
    targets,_,_=regions_for('A'*1000,[gene(300,700,conf=30)],'conservative',500,0,300,99)
    assert targets==[(0,1000)]
    # The explicitly aggressive policy illustrates the exact false-negative risk.
    assert regions_for('A'*1000,[gene(300,700,conf=30)],'all',500,0,300,99)[0]==[]


@pytest.mark.parametrize('g',[gene(0,1000,partial='10'),gene(0,1000,conf=float('nan')),
                             {**gene(0,1000),'attrs':{}}])
def test_uncertain_predictions_are_not_hard_exclusions(g):
    assert regions_for('A'*1000,[g],'conservative',40,0,300,99)[0]==[(0,1000)]


def test_uncertain_overlapping_cds_rescues_high_scoring_cds():
    targets,_,_=regions_for('A'*1000,[gene(0,1000),gene(400,600,conf=30,gid='weak')],
                            'conservative',40,20,300,99)
    assert targets==[(380,620)]


def test_minimum_length_is_applied_before_padding():
    targets,rows,_=regions_for('A'*1000,[gene(60,1000)],'all',60,75,300,99)
    assert targets==[(0,135)]
    assert regions_for('A'*1000,[gene(59,1000)],'all',60,75,300,99)[0]==[]


def test_ambiguous_bases_are_never_bridged_by_padding():
    targets,_,_=regions_for('A'*100+'N'*10+'C'*100,[],'conservative',40,75,300,99)
    assert targets==[(0,100),(110,210)]


def test_joined_gff_and_negative_strand_preserve_noncds_gaps(tmp_path):
    gff=tmp_path/'x.gff'
    gff.write_text('c\tt\tCDS\t11\t30\t.\t-\t0\tID=g;product=a%3Bb\n'
                   'c\tt\tCDS\t51\t70\t.\t-\t0\tID=g\n')
    g=read_gff(gff,{'c':'A'*100})[0]
    assert g['parts']==[(10,30),(50,70)] and g['product']=='a;b'
    assert regions_for('A'*100,[g],'all',0,0,300,99)[0]==[(0,10),(30,50),(70,100)]
    rec=make_record('n','c','A'*100,0,100,35,45,[g],{},[],[])
    feature=next(f for f in rec.features if f.type=='CDS')
    assert str(feature.extract(rec.seq))=='T'*40


def test_bad_gff_contig_or_out_of_bounds_rejected(tmp_path):
    gff=tmp_path/'x.gff'
    gff.write_text('c\tt\tCDS\t1\t101\t.\t+\t0\tID=g\n')
    with pytest.raises(ValueError): read_gff(gff,{'c':'A'*100})


def test_audit_selection_order_independent_and_explicit():
    assert choose_audit(['a','b','c'],.02,17,[])==choose_audit(['c','a','b'],.02,17,[])
    assert len(choose_audit(['a','b','c'],.02,17,[]))==1
    assert choose_audit(['a','b'],0,17,['b'])=={'b'}
    with pytest.raises(ValueError): choose_audit(['a'],0,17,['missing'])


def fixture_prepared(tmp_path):
    rng=random.Random(5)
    seq=''.join(rng.choices('ACGT',k=12000))
    fa=tmp_path/'genome.fa'; fa.write_text('>c\n'+seq+'\n')
    gff=tmp_path/'genes.gff3'
    # 600 nt noncoding target, then a large CDS. Audit keeps the full contig.
    gff.write_text('c\tt\tCDS\t601\t12000\t.\t+\t0\tID=g;conf=100;cscore=100;partial=00\n')
    prep=tmp_path/'prepared'
    cli('prepare_region_scan.py','--fasta',fa,'--gff',gff,'--out',prep,
        '--boundary-pad',0,'--audit-fraction',1)
    return fa,gff,prep,seq


def test_cli_keeps_original_windows_and_reports_actual_work(tmp_path):
    fa,gff,prep,seq=fixture_prepared(tmp_path)
    rows=table(prep/'windows.tsv')
    assert len(rows)==7
    assert sum(r['policy_selected']=='True' for r in rows)==1
    from make_windows import read_fasta
    for name,s in read_fasta(prep/'windows.fa'):
        r=next(x for x in rows if x['window_id']==name)
        assert s==seq[int(r['start0']):int(r['end0'])]
    summary=json.loads((prep/'summary.json').read_text())
    assert summary['window_nt_reduction_fraction']==0
    assert summary['policy_window_nt_reduction_fraction']==pytest.approx(6/7)


def fake_tracks(prep,path):
    rows=table(prep/'windows.tsv')
    path.mkdir()
    meta=dict(revision='test',mode='random_partition',groups=20,repeats=2,seed=20261014,
              attention_mask_policy='always',records=[])
    for i,r in enumerate(rows):
        length=int(r['length_nt']); start=int(r['start0'])
        for ori in ('+','-'):
            v=np.zeros(length)
            # Two loci: first is a target, second lies deep in a filtered CDS.
            for s,e in [(100,160),(300,360),(7100,7160),(7300,7360)]:
                x,y=max(0,s-start),min(length,e-start)
                if x<y: v[x:y]=1.2
            name=f'{i}{ori}.npz'; np.savez_compressed(path/name,above_uniform=v)
            meta['records'].append(dict(id=r['window_id'],scoring_orientation=ori,track=name,
                sequence_sha256=r['window_sha256'],coordinates='input_fasta_record_0_based'))
    (path/'run.json').write_text(json.dumps(meta))


def test_synthetic_end_to_end_audit_explainer_and_neighborhood(tmp_path):
    fa,gff,prep,seq=fixture_prepared(tmp_path)
    scores=tmp_path/'scores'; fake_tracks(prep,scores)
    called=tmp_path/'called'
    cli('call_region_candidates.py','--prepared',prep,'--scores',scores,'--threshold',.8,'--out',called)
    rows=table(called/'candidates.tsv')
    assert len(rows)==2
    assert {r['selection_class'] for r in rows}=={'target_associated','audit_only'}
    summary=json.loads((called/'summary.json').read_text())
    assert summary['audit_full_without_policy_overlap']==1
    ex=tmp_path/'explainer'
    cli('explain_repeats_v3.py','--candidates',called/'candidates.tsv','--genomes',fa,'--out',ex)
    out=tmp_path/'neighborhood'
    cli('repeat_neighborhoods.py','--fasta',fa,'--gff',gff,'--candidates',called/'candidates.tsv',
        '--copies',ex/'copies.tsv','--models',ex/'models.tsv','--members',called/'members.tsv',
        '--flank',500,'--out',out)
    # Overlapping neighborhoods remain two records; original sequence preserved.
    from Bio import SeqIO
    regions=table(out/'neighborhoods.tsv')
    assert len(regions)==2
    for r in regions:
        record=SeqIO.read(out/r['genbank'],'genbank')
        assert str(record.seq)==seq[int(r['start0']):int(r['end0'])]


@pytest.mark.parametrize('failure',['missing_strand','hash','length','mixed','coordinates','nan'])
def test_incomplete_or_incompatible_scores_fail(failure):
    windows={'w':dict(start0='0',end0='100',window_sha256='x')}
    meta=dict(revision='r',mode='exact',groups=20,repeats=2,seed=1,attention_mask_policy='always')
    rec=dict(sequence_sha256='x',coordinates='input_fasta_record_0_based')
    tracks={('w',ori):(np.zeros(100),dict(meta),dict(rec)) for ori in ('+','-')}
    if failure=='missing_strand': del tracks[('w','-')]
    elif failure=='hash': tracks[('w','-')][2]['sequence_sha256']='y'
    elif failure=='length': tracks[('w','-')]=(np.zeros(99),meta,rec)
    elif failure=='mixed': tracks[('w','-')][1]['revision']='other'
    elif failure=='coordinates': tracks[('w','-')][2]['coordinates']='reverse'
    elif failure=='nan': tracks[('w','-')][0][0]=np.nan
    with pytest.raises(ValueError): validate_tracks(tracks,windows)


def test_transitive_merge_is_flagged_and_members_retained():
    raw=[dict(contig='c',start0=s,end0=e,window_id='w') for s,e in [(0,100),(90,190),(180,280)]]
    merged=[dict(contig='c',start0=0,end0=280)]
    members=decorate(merged,raw,{'c':SpanIndex([(0,300)])},set())
    assert merged[0]['merge_geometry']=='transitive_chain' and len(members)==3
    assert merged[0]['n_windows']==1


def test_frozen_threshold_checks_model_and_window_context(tmp_path):
    _,_,prep,_=fixture_prepared(tmp_path)
    scores=tmp_path/'scores'; fake_tracks(prep,scores)
    from peak_cluster import PARAMS
    threshold=tmp_path/'thresholds.json'
    data=dict(model_revision='test',mode='random_partition',groups=20,repeats=2,
              seed=20261014,attention_mask_policy='always',peak_cluster_params=PARAMS,
              thresholds={'bg100_per_Mb':{'tau':.8}})
    threshold.write_text(json.dumps(data))
    cli('call_region_candidates.py','--prepared',prep,'--scores',scores,
        '--threshold-file',threshold,'--out',tmp_path/'ok')
    data['model_revision']='wrong-model'; threshold.write_text(json.dumps(data))
    with pytest.raises(subprocess.CalledProcessError) as ex:
        cli('call_region_candidates.py','--prepared',prep,'--scores',scores,
            '--threshold-file',threshold,'--out',tmp_path/'bad-model')
    assert 'do not match' in ex.value.stderr
    data['model_revision']='test'; threshold.write_text(json.dumps(data))
    prov=prep/'provenance.json'; values=json.loads(prov.read_text())
    values['parameters']['window']=1000; prov.write_text(json.dumps(values))
    with pytest.raises(subprocess.CalledProcessError) as ex:
        cli('call_region_candidates.py','--prepared',prep,'--scores',scores,
            '--threshold-file',threshold,'--out',tmp_path/'bad-context')
    assert 'original 3000/1500/40' in ex.value.stderr


def test_neighborhood_nearest_genes_and_clipped_cds_translation():
    genes=[gene(10,100,gid='left'),gene(900,1000,gid='right')]
    assert neighborhood(400,450,2000,genes,50,1)==(10,1000)
    rec=make_record('n','c','A'*2000,50,950,400,450,genes,{'left':'X','right':'Y'},[],[])
    assert all('translation' not in f.qualifiers for f in rec.features)


def test_empty_array_shard_does_not_load_model(tmp_path):
    prep=tmp_path/'prep'; prep.mkdir(); (prep/'records.txt').write_text('only-one\n')
    out=tmp_path/'gpu'
    result=cli('score_region_chunk.py','--prepared',prep,'--out',out,'--snapshot','does-not-exist',
               '--chunks',2,'--task',2,'--trust-reviewed-code')
    assert 'Empty shard' in result.stdout and not out.exists()


def test_zero_selected_windows_complete_without_scores(tmp_path):
    fa=tmp_path/'g.fa'; fa.write_text('>c\n'+'A'*1000+'\n')
    gff=tmp_path/'g.gff'; gff.write_text('c\tt\tCDS\t1\t1000\t.\t+\t0\tID=g\n')
    prep=tmp_path/'prep'
    cli('prepare_region_scan.py','--fasta',fa,'--gff',gff,'--out',prep,'--policy','all','--audit-fraction',0)
    out=tmp_path/'called'
    cli('call_region_candidates.py','--prepared',prep,'--out',out,'--threshold',.8)
    assert table(out/'candidates.tsv')==[]
    assert json.loads((out/'summary.json').read_text())['unique_scanned_nt']==0


def test_real_predictor_smoke(tmp_path):
    pytest.importorskip('pyrodigal')
    rng=random.Random(42)
    seq=''.join(rng.choices('ACGT',k=10000))+'NN'+'ATG'+'GCC'*100+'TAA'
    fa=tmp_path/'x.fa'; fa.write_text('>c\n'+seq+'\n')
    out=tmp_path/'cds'
    cli('predict_cds.py','--fasta',fa,'--out',out)
    for g in read_gff(out/'genes.gff3',{'c':seq}):
        assert 'N' not in seq[g['start0']:g['end0']]
    assert table(out/'contigs.tsv')[0]['length_nt']==str(len(seq))
