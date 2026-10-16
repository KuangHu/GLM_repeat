#!/usr/bin/env python3
"""CPU-only length sweep on frozen CDS predictions; labels are evaluation only."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

from region_io import new_dir, table, write_table, provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--gff',type=Path,required=True)
    p.add_argument('--loci',type=Path,required=True)
    p.add_argument('--units',type=Path,help='Optional annotated repeat copies, same coordinate columns')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--lengths',type=int,nargs='+',default=[0,40,60,100])
    a=p.parse_args()
    if any(x<0 for x in a.lengths) or len(a.lengths)!=len(set(a.lengths)):
        p.error('Lengths must be nonnegative and distinct')
    out=new_dir(a.out)
    rows=[]
    for length in a.lengths:
        dest=out/f'min{length}'
        subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_region_scan.py')),
            '--fasta',str(a.fasta),'--gff',str(a.gff),'--loci',str(a.loci),'--out',str(dest),
            '--min-noncoding',str(length),'--audit-fraction','0'],check=True,capture_output=True,text=True)
        summary=json.loads((dest/'summary.json').read_text())
        retention=table(dest/'locus_retention.tsv')
        row=dict(min_noncoding=length,n_loci=len(retention),
            loci_target_overlap=sum(int(r['policy_target_overlap_nt'])>0 for r in retention),
            loci_fully_window_covered=sum(r['policy_fully_covered']=='True' for r in retention),
            n_full_windows=summary['n_full_windows'],n_policy_windows=summary['n_policy_windows'],
            window_nt_reduction_fraction=summary['policy_window_nt_reduction_fraction'])
        if a.units:
            from region_io import SpanIndex, fasta, interval
            seqs=fasta(a.fasta)
            by={c:[] for c in seqs}
            for w in table(dest/'windows.tsv'):
                by[w['contig']].append((int(w['start0']),int(w['end0'])))
            indexes={c:SpanIndex(v) for c,v in by.items()}
            units=[]
            for u in table(a.units):
                c,s,e=interval(u,seqs)
                units.append({**u,'window_covered_nt':indexes[c].overlap(s,e),
                              'fully_window_covered':indexes[c].overlap(s,e)==e-s})
            write_table(dest/'unit_retention.tsv',units,list(units[0]) if units else ['contig','start0','end0'])
            row.update(n_units=len(units),units_fully_window_covered=sum(u['fully_window_covered'] for u in units))
        rows.append(row)
    write_table(out/'comparison.tsv',rows,list(rows[0]))
    provenance(out,a,[a.fasta,a.gff,a.loci,a.units],
               note='No threshold is selected by this script. Coverage is not detector recall. '
                    'Keep independent held-out genomes for generalization assessment.')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    main()
