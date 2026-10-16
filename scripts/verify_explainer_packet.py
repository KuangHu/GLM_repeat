#!/usr/bin/env python3
"""Verify a standalone explainer v0.3 packet and optionally regenerate its evidence.

The manifest is checked with Python's standard library before dependencies load.
Rebuilding requires only the packet, numpy, Biopython and (for tests) pytest.
No repository checkout, external genomes, model, network or GPU is used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath


def verify_manifest(root):
    expected = {}
    for line in (root/'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
        digest,name = line.split('  ',1)
        posix = PurePosixPath(name)
        if (not re.fullmatch('[0-9a-f]{64}',digest) or posix.is_absolute()
                or '..' in posix.parts or '\\' in name or name in expected):
            raise ValueError('Invalid or duplicate manifest entry: '+name)
        path = root.joinpath(*posix.parts)
        path.resolve().relative_to(root.resolve())
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError('Manifest mismatch: '+name)
        expected[name]=digest
    if not expected:
        raise ValueError('Empty manifest')
    return expected


def rows(path):
    with path.open(encoding='utf-8',newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        return reader.fieldnames,list(reader)


def same_cell(a,b):
    if a==b:
        return True
    # Only decimal/scientific cells can use tolerance; counts/DNA/IDs stay exact.
    numeric=r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'
    if (not all(re.fullmatch(numeric,x) for x in (a,b))
            or not all(any(mark in value for mark in '.eE') for value in (a,b))):
        return False
    x,y=float(a),float(b)
    return math.isfinite(x) and math.isfinite(y) and math.isclose(x,y,rel_tol=1e-12,abs_tol=1e-12)


def compare_table(expected,actual,ignore=()):
    ef,er=rows(expected); af,ar=rows(actual)
    if ef!=af or len(er)!=len(ar):
        raise ValueError('Table schema/row count differs: '+str(actual))
    for i,(e,a) in enumerate(zip(er,ar),2):
        for field in ef:
            if field not in ignore and not same_cell(e[field],a[field]):
                raise ValueError(f'Table differs: {actual}:{i}, field={field}')
    return dict(rows=len(er),columns=len(ef),ignored_fields=list(ignore))


def compare_json(expected,actual):
    a=json.loads(expected.read_text(encoding='utf-8'))
    b=json.loads(actual.read_text(encoding='utf-8'))
    # Physical source-file bytes can differ after a Git checkout; this field is
    # provenance, while actual sequence/annotation sources are checked separately.
    if isinstance(a,dict):
        a.pop('code_sha256',None); b.pop('code_sha256',None)
    if a!=b:
        raise ValueError('JSON evidence differs: '+str(actual))


def check_sequences(root):
    from Bio import SeqIO
    records=list(SeqIO.parse(root/'data/genomes.fa','fasta'))
    seqs={r.id:str(r.seq).upper() for r in records}
    annotations=list(SeqIO.parse(root/'data/genomes.gb','genbank'))
    gb={r.id:str(r.seq).upper() for r in annotations}
    if len(seqs)!=len(records) or len(gb)!=len(annotations) or seqs!=gb:
        raise ValueError('FASTA/GenBank IDs or original DNA do not agree')
    sources=json.loads((root/'data/SOURCES.json').read_text(encoding='utf-8'))
    expected={r['accession']:r for r in sources['records']}
    if set(expected)!=set(seqs):
        raise ValueError('Source inventory does not match genomes')
    for name,seq in seqs.items():
        item=expected[name]
        if len(seq)!=item['length_nt'] or hashlib.sha256(seq.encode()).hexdigest()!=item['sequence_sha256']:
            raise ValueError('Sequence inventory mismatch: '+name)
    for name in ('controls','heldout'):
        _,rr=rows(root/f'outputs/inputs/{name}.tsv')
        for r in rr:
            if r['contig'] not in seqs or not 0<=int(r['start0'])<int(r['end0'])<=len(seqs[r['contig']]):
                raise ValueError('Candidate outside source sequence')
    return dict(n_records=len(seqs),total_nt=sum(map(len,seqs.values())),fasta_equals_genbank=True)


def execute(root,script,*args):
    result=subprocess.run([sys.executable,str(root/'scripts'/script),*map(str,args)],
                          cwd=root,check=True,capture_output=True,text=True,encoding='utf-8')
    print(f'{script}: completed',flush=True)
    return result.stdout


def rebuild(root,dest):
    if dest.exists() and any(dest.iterdir()):
        raise ValueError('Rebuild directory must be new or empty: '+str(dest))
    dest.mkdir(parents=True,exist_ok=True)
    comparisons={}
    for name in ('controls','heldout'):
        annotated=dest/'inputs'/f'{name}.tsv'
        execute(root,'add_cds_context.py','--candidates',root/f'data/candidates_pre_cds/{name}.tsv',
                '--genbank',root/'data/genomes.gb','--out',annotated)
        comparisons[f'inputs/{name}.tsv']=compare_table(root/f'outputs/inputs/{name}.tsv',annotated,
                                                        ignore=('CDS_annotation_source',))
        folder='controls' if name=='controls' else 'heldout_primary'
        execute(root,'explain_repeats_v3.py','--candidates',annotated,
                '--genomes',root/'data/genomes.fa','--out',dest/folder)
        expected=root/'outputs'/folder
        for path in sorted(expected.glob('*.tsv')):
            comparisons[f'{folder}/{path.name}']=compare_table(path,dest/folder/path.name)
        compare_json(expected/'summary.json',dest/folder/'summary.json')
        # Path keys differ; hashes of regenerated candidate bytes legitimately
        # differ because the descriptive annotation-source path was relocated.
        sources=json.loads((dest/folder/'inputs.sha256.json').read_text(encoding='utf-8'))
        if sources[str(root/'data/genomes.fa')]!=hashlib.sha256((root/'data/genomes.fa').read_bytes()).hexdigest():
            raise ValueError('Regenerated genome provenance hash differs')
    execute(root,'validate_explainer_v3.py','--controls',dest/'controls','--heldout',dest/'heldout_primary',
            '--genomes',root/'data/genomes.fa','--control-definitions',root/'data/controls_definition.tsv',
            '--annotated-units',root/'data/annotated_units.tsv','--v2-heldout',root/'data/v2_heldout',
            '--out',dest/'validation')
    for name in ('candidate_changes.tsv','controls_overlap.tsv'):
        comparisons['validation/'+name]=compare_table(root/'outputs/validation'/name,dest/'validation'/name)
    compare_json(root/'outputs/validation/validation.json',dest/'validation/validation.json')
    execute(root,'make_review_set_v3.py','--results',dest/'heldout_primary','--genomes',root/'data/genomes.fa',
            '--previous-review',root/'data/previous_review.tsv','--fixed-random',root/'data/fixed_random.tsv',
            '--out',dest/'review_set')
    for path in sorted((root/'outputs/review_set').glob('*.tsv')):
        comparisons['review_set/'+path.name]=compare_table(path,dest/'review_set'/path.name)
    compare_json(root/'outputs/review_set/summary.json',dest/'review_set/summary.json')
    if (root/'outputs/review_set/regions.fa').read_text()!= (dest/'review_set/regions.fa').read_text():
        raise ValueError('Review-region original DNA differs')
    return dict(compared_tables=comparisons,validation=json.loads((dest/'validation/validation.json').read_text()),
                ignored_provenance=['annotation source path','summary source-code byte hash',
                                    'path-keyed input/provenance JSONs','historical package_check.json'],
                numeric_comparison='same rows/columns/order; counts and text exact; decimal rel/abs tolerance 1e-12')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--packet',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--rerun',type=Path,help='Rebuild all evidence in a new directory, then compare to frozen outputs')
    p.add_argument('--run-tests',action='store_true')
    a=p.parse_args(); root=a.packet.resolve()
    manifest=verify_manifest(root)
    print(f'Manifest: {len(manifest)} exact-byte entries passed',flush=True)
    report=dict(manifest_entries=len(manifest),sequences=check_sequences(root),
                environment=dict(python=platform.python_version(),**{name:importlib.metadata.version(name)
                    for name in ('numpy','biopython')}),biological_validation=False)
    if a.run_tests:
        test=subprocess.run([sys.executable,'-m','pytest','-q','tests/'],cwd=root,
                            check=True,capture_output=True,text=True,encoding='utf-8')
        report['tests']=test.stdout.strip()
        print(report['tests'],flush=True)
    if a.rerun:
        dest=a.rerun.resolve()
        # Prevent overwriting the frozen evidence, sources or packaged code.
        if dest==root or dest in root.parents or any(dest==root/x or root/x in dest.parents
               for x in ('outputs','data','scripts','tests','environment')):
            raise ValueError('Unsafe rebuild destination')
        report['reproduction']=rebuild(root,dest)
        report['passed']=True
        (dest/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='reproduction'},indent=2))


if __name__=='__main__':
    main()
