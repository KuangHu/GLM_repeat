#!/usr/bin/env python3
"""Select a bounded array shard, then invoke the unchanged frozen GPU scorer."""
import argparse
import subprocess
import sys
from pathlib import Path

from make_windows import read_fasta
from region_io import new_dir


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--snapshot',type=Path,required=True)
    p.add_argument('--chunks',type=int,required=True)
    p.add_argument('--task',type=int,required=True)
    p.add_argument('--mode',choices=['exact','random_partition'],default='random_partition')
    p.add_argument('--trust-reviewed-code',action='store_true')
    a = p.parse_args()
    if not 1<=a.task<=a.chunks:
        p.error('Need 1 <= task <= chunks')
    if not a.trust_reviewed_code:
        p.error('Review model snapshot code before passing --trust-reviewed-code')
    ids = (a.prepared/'records.txt').read_text(encoding='utf-8').splitlines()
    if len(ids)!=len(set(ids)):
        raise ValueError('Duplicate record IDs')
    wanted = set(ids[a.task-1::a.chunks])
    # Essential: empty --records in timing_score would otherwise score ALL windows.
    if not wanted:
        print('Empty shard: no GPU model loaded and no scores emitted.')
        return
    part = new_dir(a.out/f'part{a.task:05d}')
    found = set()
    with (part/'input.fa').open('w',encoding='utf-8') as fh:
        for name,seq in read_fasta(a.prepared/'windows.fa'):
            if name in wanted:
                fh.write(f'>{name}\n{seq.lower()}\n')
                found.add(name)
    if found!=wanted:
        raise ValueError('Record list and prepared FASTA disagree')
    subprocess.run([sys.executable,str(Path(__file__).with_name('timing_score.py')),
        '--snapshot',str(a.snapshot),'--fasta',str(part/'input.fa'),'--out',str(part/'scores'),
        '--mode',a.mode,'--groups','20','--repeats','2','--seed','20261014',
        '--attention-mask-policy','always','--both-strands','--batch-size','1',
        '--device','cuda:0','--trust-reviewed-code'],check=True)


if __name__ == '__main__':
    main()
