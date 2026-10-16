"""Semantic regressions from the independent review; no GPU or downloaded data."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
spec = importlib.util.spec_from_file_location('ex3',SCRIPTS/'explain_repeats_v3.py')
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


def test_distances_and_support_refer_to_final_consensus():
    strings = ['ACGTACGTACGA','ACGTACGTACGT','ACGTACGTACGT']
    p = ex.profile(strings,ex.background(''.join(strings)))
    assert p['consensus']==strings[1]
    assert p['mismatches']==[1,0,0]
    assert p['consensus_identity'][1]==1
    assert p['loo_identity'][0]==11/12


def test_unequal_lengths_are_not_silently_compared():
    with pytest.raises(ValueError):
        ex.mismatches('ACG','AC')


def test_cores_survive_rejected_optional_extension():
    # Four exact cores with a supported flank consensus and one discordant flank.
    core='ACGTCAGTGCAT'
    seq=('GGGG'+core+'CCCC'+'AT'*15)*3 + 'TTTT'+core+'AAAA'
    models=ex.motif_models(seq)
    m=next(m for m in models if m['search_motif']==core)
    assert len(m['copies'])==4
    assert len(m['extensions'])==4
    assert any(not x['accepted'] for x in m['extensions'])
    for (a,b),x in zip(m['copies'],m['extensions']):
        assert seq[a:b]==core
        assert x['start']<=a<b<=x['end']


def test_extension_never_recruits_unlinked_copies():
    seq='GGACGTCAGTGCATCC'+'T'*35+'GGACGTCAGTGCATCC'+'A'*35+'GGACGTCAGTGCATCC'
    for m in ex.motif_models(seq):
        assert len(m['extensions']) in (0,len(m['copies']))
        for (a,b),x in zip(m['copies'],m['extensions']):
            assert x['start']<=a and x['end']>=b


def test_flank_repeat_is_not_candidate_associated():
    seq='ACGTTGCACTGA'*4+''.join(np.random.default_rng(1).choice(list('ACGT'),200))
    models,diagnostics,status=ex.analyze(seq,100,150)
    assert models and status=='flank_only'
    assert all(m['n_overlapping']==0 for m in models)


def test_association_uses_overlap_and_full_containment_separately():
    seq='ACGTTGCACTGA'*4
    models,_,status=ex.analyze(seq,5,7)
    assert status=='candidate_associated'
    assert models[0]['n_overlapping']>0
    assert models[0]['n_contained']==0
    assert models[0]['candidate_coverage_nt']==2


def test_all_periodic_intervals_retained_instead_of_longest_flank_only():
    rng=np.random.default_rng(42)
    seq='ACGTT'*10+''.join(rng.choice(list('ACGT'),200))+'ACGTT'*30
    intervals=ex.periodic_intervals(seq,5)
    assert len(intervals)>=2
    assert any(a<50 and b>0 for a,b in intervals)
    assert any(b>250 for a,b in intervals)


def test_periodic_interval_no_extra_base_at_right_boundary():
    seq='ACGTT'*10
    assert ex.periodic_intervals(seq,5)==[(0,len(seq))]


def test_degenerate_period_and_harmonics_are_retained():
    rng=np.random.default_rng(3)
    seq=''.join('C'+''.join(rng.choice(list('ACGT'),2))+'GG' for _ in range(60))
    ds={d['period']:d for d in ex.period_diagnostics(seq)}
    assert ds[5]['passes_identity_floor']
    assert ds[10]['passes_identity_floor']
    assert '5' in ds[10]['supported_divisors'].split(',')
    models,_,_=ex.analyze(seq,0,len(seq))
    assert any(m['period']==5 for m in models)
    assert models[0]['period']==5
    assert all(m['supported_divisor']==5 for m in models if m['period'] in (10,15))


def test_noise_argmax_not_promoted_to_periodic_model():
    rng=np.random.default_rng(5)
    seq=''.join(rng.choice(list('ACGT'),1000))
    ds=ex.period_diagnostics(seq)
    assert ds and not any(d['passes_identity_floor'] for d in ds)
    models,_,_=ex.analyze(seq,0,len(seq))
    assert not any(m['source']=='periodic' for m in models)


def test_phase_is_a_display_convention_not_fitted_confidence(tmp_path):
    ex.run([dict(contig='test',start0=0,end0=100)],{'test':'ACGTT'*20},tmp_path,flank=0)
    rows=ex.read_table(tmp_path/'hypotheses.tsv')
    periodic=[r for r in rows if r['model_source']=='periodic']
    assert periodic
    assert all(r['phase_status']=='unanchored_cyclic' and r['display_phase']=='0' for r in periodic)
    assert all('phase_info_margin_bits' not in r for r in rows)


def test_ambiguity_splits_without_inventing_A_bases(tmp_path):
    seq='ACGTTGCACTGA'*2+'N'*30+'ACGTTGCACTGA'*2
    ex.run([dict(contig='x',start0=0,end0=len(seq))],{'x':seq},tmp_path,flank=0)
    for r in ex.read_table(tmp_path/'copies.tsv'):
        a,b=int(r['start0']),int(r['end0'])
        assert seq[a:b]==r['sequence'] and 'N' not in r['sequence']
    assert ex.read_table(tmp_path/'candidates.tsv')[0]['ambiguous_bases']=='30'


def test_two_copies_do_not_establish_regular_spacing():
    kind,pitch,cv=ex.classify([(0,39),(143,182)],39)
    assert kind=='repeat_pair' and pitch==143 and cv is None
    p=ex.profile(['ACGT'*9+'ACG']*2,np.full(4,0.25))
    assert all(x>0 for x in p['loo_enrichment'])


def test_empty_rerun_overwrites_tables_and_missing_peaks_are_blank(tmp_path):
    seq='ACGTTGCACTGA'*4
    ex.run([dict(contig='x',start0=0,end0=len(seq))],{'x':seq},tmp_path,flank=0)
    assert all(r['overlapping_peaks']=='' for r in ex.read_table(tmp_path/'copies.tsv'))
    ex.run([],{},tmp_path)
    assert ex.read_table(tmp_path/'copies.tsv')==[]
    assert ex.read_table(tmp_path/'models.tsv')==[]


def test_copy_ids_and_extension_links_round_trip(tmp_path):
    seq=('GGGGACGTCAGTGCATCCCC'+'AT'*15)*3+'TTTTACGTCAGTGCATAAAA'
    ex.run([dict(contig='x',start0=4,end0=len(seq)-4)],{'x':seq},tmp_path,flank=4)
    models={m['model_id']:m for m in ex.read_table(tmp_path/'hypotheses.tsv')}
    copies={c['copy_id']:c for c in ex.read_table(tmp_path/'copies.tsv')}
    for c in copies.values():
        a,b=int(c['start0']),int(c['end0'])
        assert seq[a:b]==c['sequence']
        assert int(c['mismatches_to_consensus'])==ex.mismatches(c['sequence'],models[c['model_id']]['consensus'])
    for x in ex.read_table(tmp_path/'extensions.tsv'):
        core=copies[x['core_copy_id']]
        assert int(x['start0'])<=int(core['start0'])<int(core['end0'])<=int(x['end0'])
        assert seq[int(x['start0']):int(x['end0'])]==x['sequence']
        assert (x['accepted']=='True')==(int(x['search_mismatches'])<=int(x['search_budget']))


def test_missing_sequence_is_explicit_unresolved(tmp_path):
    ex.run([dict(contig='missing',start0=10,end0=20)],{},tmp_path)
    row=ex.read_table(tmp_path/'candidates.tsv')[0]
    assert row['status']=='unresolved' and row['reason']=='sequence_unavailable'


def test_cds_joined_locations_do_not_cover_the_gap():
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature,FeatureLocation,CompoundLocation
    from Bio.Seq import Seq
    from add_cds_context import cds_context
    record=SeqRecord(Seq('A'*100))
    record.features=[SeqFeature(CompoundLocation([FeatureLocation(0,10),FeatureLocation(30,40)]),type='CDS')]
    assert cds_context(record,10,30)['n_CDS_overlapping']==0
    assert cds_context(record,5,35)['CDS_overlap_nt']==10
    assert cds_context(None,0,10)['annotation_available'] is False
