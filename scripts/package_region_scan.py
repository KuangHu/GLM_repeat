#!/usr/bin/env python3
"""Build a Linux-ready source bundle; excludes models, virtualenvs and raw genomes."""
import argparse
import hashlib
import io
import json
import tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=['region_io.py','predict_cds.py','prepare_region_scan.py','benchmark_region_selection.py',
         'score_region_chunk.py','call_region_candidates.py','repeat_neighborhoods.py',
         'make_windows.py','timing_score.py','reference_score.py','compare_tracks.py',
         'summarize_scores.py','peak_cluster.py','pilot_scan.py','explain_repeats_v3.py',
         'add_cds_context.py','validate_explainer_v3.py','package_region_scan.py']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'dist'/'GLM_repeat_region_scan_v0.1.tar.gz')
    a=p.parse_args()
    if a.out.exists():
        p.error('Archive already exists; choose a new path')
    files=[ROOT/'scripts'/s for s in SCRIPTS]
    files += sorted((ROOT/'tests').glob('test_*.py'))
    files += [ROOT/s for s in ['REGION_SCAN_zh.md','requirements-region-scan.txt','requirements-bootstrap.txt',
                'requirements-review.txt','slurm/region_score_lrc.sbatch',
                'results/pilot/calibration/thresholds/thresholds.json',
                'results/explainer_v3/EXPLAINER_V3.md']]
    # Small factual CPU validation tables; no private or raw genome material.
    files += [ROOT/s for s in ['results/region_scan_validation/length_sweep/comparison.tsv',
                'results/region_scan_validation/prepared/summary.json',
                'results/region_scan_validation/prepared/locus_retention.tsv',
                'results/region_scan_validation/validation.json'] if (ROOT/s).exists()]
    data={str(f.relative_to(ROOT)).replace('\\','/'):f.read_bytes() for f in files}
    # Shell scripts must be executable Linux text, even if checked out on Windows.
    for name in data:
        if name.endswith(('.sbatch','.sh')):
            data[name]=data[name].replace(b'\r\n',b'\n')
    data['MANIFEST.sha256']=''.join(hashlib.sha256(v).hexdigest()+'  '+k+'\n'
                                  for k,v in sorted(data.items())).encode()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(a.out,'w:gz') as tf:
        for name,content in sorted(data.items()):
            entry=tarfile.TarInfo('GLM_repeat_region_scan_v0.1/'+name)
            entry.size=len(content)
            entry.mode=0o755 if name.endswith('.sbatch') else 0o644
            entry.mtime=0
            tf.addfile(entry,io.BytesIO(content))
    print(json.dumps(dict(archive=str(a.out),files=len(data),bytes=a.out.stat().st_size,
                         sha256=hashlib.sha256(a.out.read_bytes()).hexdigest()),indent=2))


if __name__=='__main__':
    main()
