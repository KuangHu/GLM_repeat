#!/usr/bin/env python3
"""Predict CDSs on whole contigs; emit GFF3 + proteins, without filtering DNA."""
from __future__ import annotations

import argparse
import importlib.metadata
import re
from pathlib import Path
from urllib.parse import quote

from region_io import fasta, new_dir, provenance, write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--caller', choices=['pyrodigal', 'pyrodigal-gv'], default='pyrodigal')
    p.add_argument('--mode', choices=['single', 'meta'], default='meta',
                   help='single: input must contain only contigs from ONE genome')
    p.add_argument('--translation-table', type=int, default=None,
                   help='single mode only; default 11. Meta mode selects its model/code.')
    a = p.parse_args()
    if a.caller == 'pyrodigal-gv' and a.mode != 'meta':
        p.error('Use pyrodigal-gv with meta mode')
    if a.mode == 'meta' and a.translation_table is not None:
        p.error('An explicit translation table is only supported in single mode')
    if a.caller == 'pyrodigal-gv':
        from pyrodigal_gv import ViralGeneFinder as GeneFinder
    else:
        from pyrodigal import GeneFinder
    seqs = fasta(a.fasta)
    # Replace ambiguous symbols with N in-place; never change sequence length.
    normalized = {c: re.sub('[^ACGT]', 'N', s) for c, s in seqs.items()}
    finder = GeneFinder(meta=a.mode == 'meta', mask=True, min_mask=1, closed=False)
    if a.mode == 'single':
        finder.train(*normalized.values(), translation_table=a.translation_table or 11)
    out = new_dir(a.out)
    contigs = []
    counter = 0
    with (out/'genes.gff3').open('w', encoding='utf-8') as gf, \
         (out/'proteins.faa').open('w', encoding='utf-8') as pf:
        gf.write('##gff-version 3\n')
        for contig, seq in normalized.items():
            genes = finder.find_genes(seq)
            code = genes.training_info.translation_table
            gf.write(f'##sequence-region {quote(contig, safe="")} 1 {len(seq)}\n')
            contigs.append(dict(contig=contig, length_nt=len(seq), n_CDS=len(genes),
                                translation_table=code))
            for gene in genes:
                counter += 1
                gid = f'cds{counter:08d}'
                s, e = gene.begin-1, gene.end
                if not 0 <= s < e <= len(seq):
                    raise ValueError('Caller returned invalid coordinates')
                attr = (f'ID={gid};partial={int(gene.partial_begin)}{int(gene.partial_end)};'
                        f'conf={gene.confidence():.6f};cscore={gene.cscore:.6f};'
                        f'transl_table={code};product=hypothetical%20protein')
                gf.write(f'{quote(contig, safe="")}\t{a.caller}\tCDS\t{s+1}\t{e}\t.\t'
                         f'{"+" if gene.strand == 1 else "-"}\t0\t{attr}\n')
                pf.write(f'>{gid}\n{gene.translate(include_stop=False)}\n')
    write_table(out/'contigs.tsv', contigs, ['contig','length_nt','n_CDS','translation_table'])
    provenance(out, a, [a.fasta], n_CDS=counter,
               versions={name: importlib.metadata.version(name) for name in
                         (['pyrodigal', 'pyrodigal-gv'] if a.caller.endswith('-gv') else ['pyrodigal'])},
               limitations=['gene confidence is not independently calibrated noncoding probability',
                            'linear contig coordinates; no circular-origin gene reconstruction'])
    print(f'{len(seqs)} contigs, {counter} CDSs -> {out}')


if __name__ == '__main__':
    main()
