#!/usr/bin/env python3
"""Build a portable CPU review ZIP from the completed v0.3 outputs."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import platform
import zipfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--genomes',type=Path,required=True)
    p.add_argument('--genbank',type=Path,required=True)
    p.add_argument('--fixed-random',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); root=a.root.resolve()
    if a.out.exists():
        p.error('Output archive already exists; choose a new path')
    payload={}
    def add(path,name):
        target=PurePosixPath(name)
        if target.is_absolute() or '..' in target.parts or name in payload:
            raise ValueError('Unsafe or duplicate archive path: '+name)
        payload[name]=Path(path).read_bytes()
    for name in ('explain_repeats_v3.py','add_cds_context.py','validate_explainer_v3.py',
                 'make_review_set_v3.py','verify_explainer_packet.py'):
        add(root/'scripts'/name,'scripts/'+name)
    add(root/'tests/test_explainer_v3.py','tests/test_explainer_v3.py')
    add(root/'requirements-review.txt','requirements-review.txt')
    add(a.genomes,'data/genomes.fa'); add(a.genbank,'data/genomes.gb')
    add(a.fixed_random,'data/fixed_random.tsv')
    add(root/'review_packet_explainer_v0.2/data/controls_definition.tsv','data/candidates_pre_cds/controls.tsv')
    add(root/'results/pilot/heldout/annotated_primary/candidates_annotated.tsv','data/candidates_pre_cds/heldout.tsv')
    from Bio import SeqIO
    fasta=list(SeqIO.parse(a.genomes,'fasta'))
    gb=list(SeqIO.parse(a.genbank,'genbank'))
    seqs={r.id:str(r.seq).upper() for r in fasta}
    annotations={r.id:str(r.seq).upper() for r in gb}
    if len(seqs)!=len(fasta) or len(annotations)!=len(gb) or seqs!=annotations:
        raise ValueError('Input FASTA and GenBank must contain the same unique records and original DNA')
    inventory=dict(provider='NCBI public accession-version records',
        annotation_policy='Archived GenBank snapshot; annotations may differ from newer live records',
        records=[dict(accession=name,length_nt=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                      url='https://www.ncbi.nlm.nih.gov/nuccore/'+name) for name,seq in seqs.items()])
    payload['data/SOURCES.json']=(json.dumps(inventory,indent=2)+'\n').encode('utf-8')
    versions={name:importlib.metadata.version(name) for name in ('numpy','biopython','pytest')}
    payload['environment/requirements-tested.txt']=''.join(f'{name}=={version}\n' for name,version in versions.items()).encode()
    payload['environment/build.json']=(json.dumps(dict(python=platform.python_version(),packages=versions,
        dependencies='CPU only; numpy, Biopython, pytest; no torch/model/network at verification time'),indent=2)+'\n').encode()
    packet=root/'review_packet_explainer_v0.2'
    for src,dst in [('data/controls_definition.tsv','data/controls_definition.tsv'),
                    ('data/annotated_units.tsv','data/annotated_units.tsv'),
                    ('outputs/heldout_primary/models.tsv','data/v2_heldout/models.tsv'),
                    ('outputs/review_set/review_set.tsv','data/previous_review.tsv')]:
        add(packet/src,dst)
    for name in ('controls','heldout_primary','inputs','validation','review_set'):
        for path in sorted((root/'results/explainer_v3'/name).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix!='.pyc':
                add(path,'outputs/'+name+'/'+path.relative_to(root/'results/explainer_v3'/name).as_posix())
    for name in ('EXPLAINER_V3.md','FIXES_zh.txt'):
        add(root/'results/explainer_v3'/name,name)
    payload['README.md']=('''# Explainer v0.3 CPU review packet

This packet includes code, 17 semantic tests, all eight original public DNA/GenBank
records, frozen candidates before and after CDS annotation, complete v0.3 evidence,
the unchanged random10 and previous review selection, and validation references.
After installing CPU dependencies, regeneration needs no repository checkout,
genome downloads, torch, model weights, GPU, or detector recalibration.
Start with FIXES_zh.txt or EXPLAINER_V3.md. The latter describes the original
repository layout; use the packet-relative commands below for this archive.

## Rebuild and compare everything

Extract into a fresh directory, cd into it, create a Python 3.10+ environment:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\\Scripts\\Activate.ps1
python -m pip install -r environment/requirements-tested.txt
python scripts/verify_explainer_packet.py --run-tests --rerun rerun
```

requirements-review.txt gives the supported dependency ranges; requirements-tested
records the exact CPU dependency versions used to build and verify this release.
Installing packages can require network access; the verification and rebuild do not.
`python scripts/verify_explainer_packet.py` checks the manifest, sequence inventory,
FASTA/GenBank identity and candidate bounds without regeneration.

The one-command rebuild regenerates CDS context from candidates_pre_cds and the
archived GenBank, then controls, all 467 heldout candidates, engineering validation,
and the review set. It compares every scientific TSV column in row order against
outputs/, along with summaries and original review DNA, and writes
rerun/reproduction.json. Decimal cells allow relative/absolute tolerance 1e-12;
counts, coordinates, DNA, labels and IDs must match exactly. The descriptive
CDS_annotation_source path is the only ignored annotation-table column. Relocated
path-keyed provenance JSONs, source-code byte hashes and historical package-check
metadata are not scientific output comparisons; fresh provenance is still emitted.
The script never rewrites outputs/, code or source data. Reuse a new/empty rerun
directory for each experiment.

Expected result: controls 15 candidates / 92 hypotheses / 1526 copy-evidence rows;
heldout 467 candidates / 1641 hypotheses / 9313 rows; candidate states
259 associated, 84 flank-only, 124 unresolved; 5991 extension attempts, 147 rejected
attempts retaining their cores; review set 48 candidates and 219 navigation groups.
These are engineering/evidence counts, not validated biological systems.

## Individual stages (use a separate new directory)

```bash
python scripts/explain_repeats_v3.py --candidates outputs/inputs/controls.tsv --genomes data/genomes.fa --out rerun/controls
python scripts/explain_repeats_v3.py --candidates outputs/inputs/heldout.tsv --genomes data/genomes.fa --out rerun/heldout_primary
python scripts/validate_explainer_v3.py --controls rerun/controls --heldout rerun/heldout_primary --genomes data/genomes.fa --control-definitions data/controls_definition.tsv --annotated-units data/annotated_units.tsv --v2-heldout data/v2_heldout --out rerun/validation
python scripts/make_review_set_v3.py --results rerun/heldout_primary --genomes data/genomes.fa --previous-review data/previous_review.tsv --fixed-random data/fixed_random.tsv --out rerun/review_set
```

MANIFEST.sha256 covers exact archived bytes, with no LF/CRLF normalization, and
excludes machine caches and the manifest itself. The archive intentionally retains
the frozen files' original line endings. Extracting the ZIP preserves these bytes;
opening/resaving a file with newline conversion can invalidate its checksum.
Existing frozen provenance retains historical source paths; data/SOURCES.json
documents accession versions, sequence lengths and sequence hashes of this snapshot.
Candidate association and annotation overlap do not establish biological function.
Alternative hypotheses share physical sequence; do not count evidence rows as
independent repeat copies. No independent natural VIPR validation is supplied.
''').encode('utf-8')
    manifest=''.join(hashlib.sha256(data).hexdigest()+'  '+name+'\n' for name,data in sorted(payload.items()))
    payload['MANIFEST.sha256']=manifest.encode('utf-8')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.out,'w',zipfile.ZIP_DEFLATED) as archive:
        for name,data in sorted(payload.items()):
            entry=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
            entry.compress_type=zipfile.ZIP_DEFLATED
            entry.create_system=3
            entry.external_attr=0o100644 << 16
            archive.writestr(entry,data)
    print(f'{a.out}: {len(payload)} files, {a.out.stat().st_size} bytes')


if __name__=='__main__':
    main()
