#!/usr/bin/env python3
"""Build a portable CPU review ZIP from the completed v0.3 outputs."""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import zipfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--genomes',type=Path,required=True)
    p.add_argument('--genbank',type=Path,required=True)
    p.add_argument('--fixed-random',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); root=a.root.resolve()
    payload={}
    def add(path,name):
        target=PurePosixPath(name)
        if target.is_absolute() or '..' in target.parts or name in payload:
            raise ValueError('Unsafe or duplicate archive path: '+name)
        payload[name]=Path(path).read_bytes()
    for name in ('explain_repeats_v3.py','add_cds_context.py','validate_explainer_v3.py','make_review_set_v3.py'):
        add(root/'scripts'/name,'scripts/'+name)
    add(root/'tests/test_explainer_v3.py','tests/test_explainer_v3.py')
    add(root/'requirements-review.txt','requirements-review.txt')
    add(a.genomes,'data/genomes.fa'); add(a.genbank,'data/genomes.gb')
    add(a.fixed_random,'data/fixed_random.tsv')
    packet=root/'review_packet_explainer_v0.2'
    for src,dst in [('data/controls_definition.tsv','data/controls_definition.tsv'),
                    ('data/annotated_units.tsv','data/annotated_units.tsv'),
                    ('outputs/heldout_primary/models.tsv','data/v2_heldout/models.tsv'),
                    ('outputs/review_set/review_set.tsv','data/previous_review.tsv')]:
        add(packet/src,dst)
    for name in ('controls','heldout_primary','inputs','validation','review_set'):
        for path in sorted((root/'results/explainer_v3'/name).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts:
                add(path,'outputs/'+name+'/'+path.relative_to(root/'results/explainer_v3'/name).as_posix())
    for name in ('EXPLAINER_V3.md','FIXES_zh.txt'):
        add(root/'results/explainer_v3'/name,name)
    payload['README.md']=('''# Explainer v0.3 CPU review packet

This packet includes code, 17 semantic tests, original public DNA/GenBank records,
the frozen candidate inputs, complete v0.3 evidence and coordinate validation.
No model weights, GPU inference or detector recalibration are needed.
Start with FIXES_zh.txt or EXPLAINER_V3.md. The latter describes the original
repository layout; use the packet-relative commands below for this archive.

```text
python -m pip install -r requirements-review.txt
python -m pytest -q tests/
python scripts/explain_repeats_v3.py --candidates outputs/inputs/controls.tsv --genomes data/genomes.fa --out rerun/controls
python scripts/explain_repeats_v3.py --candidates outputs/inputs/heldout.tsv --genomes data/genomes.fa --out rerun/heldout_primary
python scripts/validate_explainer_v3.py --controls rerun/controls --heldout rerun/heldout_primary --genomes data/genomes.fa --control-definitions data/controls_definition.tsv --annotated-units data/annotated_units.tsv --v2-heldout data/v2_heldout --out rerun/validation
python scripts/make_review_set_v3.py --results rerun/heldout_primary --genomes data/genomes.fa --previous-review data/previous_review.tsv --fixed-random data/fixed_random.tsv --out rerun/review_set
```

MANIFEST.sha256 covers exact archived bytes and excludes machine caches and the
manifest itself. Existing provenance files retain paths from their original run;
their content hashes can be compared to the corresponding archived inputs.
Candidate association and annotation overlap do not establish biological function.
Alternative hypotheses share physical sequence; do not count evidence rows as
independent repeat copies. No independent natural VIPR validation is supplied.
''').encode('utf-8')
    manifest=''.join(hashlib.sha256(data).hexdigest()+'  '+name+'\n' for name,data in sorted(payload.items()))
    payload['MANIFEST.sha256']=manifest.encode('utf-8')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.out,'w',zipfile.ZIP_DEFLATED) as archive:
        for name,data in sorted(payload.items()):
            archive.writestr(name,data)
    print(f'{a.out}: {len(payload)} files, {a.out.stat().st_size} bytes')


if __name__=='__main__':
    main()
