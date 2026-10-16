#!/usr/bin/env python3
"""Cut genome FASTA into fixed-stride scoring windows with original coordinates.

Rules (WORKFLOW_SETUP.md section 5):
  * split every contig at each non-ACGT base; never delete ambiguous bases and
    join their flanks; never join across contigs (linear contigs only in v1)
  * valid segments shorter than --min-segment are reported, not windowed
  * windows start at 0, stride, 2*stride, ... of each segment; the last emitted
    window is the first one whose end reaches the segment end, so a short tail
    already covered by an earlier window is not emitted again
  * coordinates are 0-based half-open on the original contig
No CDS or low-complexity masking is applied.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path
from typing import Iterator

VALID = re.compile(r'[ACGT]+')

WINDOW_FIELDS = ['window_id', 'genome_id', 'contig', 'start0', 'end0', 'length_nt',
                 'segment_start0', 'segment_end0', 'contig_length', 'edge_flag',
                 'gc_fraction', 'window_sha256']
SEGMENT_FIELDS = ['genome_id', 'contig', 'start0', 'end0', 'length_nt', 'status']
GENOME_FIELDS = ['genome_id', 'path', 'source_sha256', 'n_contigs', 'total_nt',
                 'acgt_nt', 'non_acgt_nt', 'short_segment_nt', 'windowed_nt',
                 'n_segments', 'n_short_segments', 'n_windows']


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    name, parts, seen = None, [], set()
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts).upper()
                fields = line[1:].split()
                if not fields or fields[0] in seen:
                    raise ValueError(f'{path}: empty or duplicated FASTA identifier.')
                name, parts = fields[0], []
                seen.add(name)
            elif name is None:
                raise ValueError(f'{path}: sequence before first header.')
            else:
                parts.append(line)
    if name is not None:
        yield name, ''.join(parts).upper()


def valid_segments(sequence: str) -> list[tuple[int, int]]:
    """Maximal ACGT-only runs as [start0, end0) on the input sequence."""
    return [(m.start(), m.end()) for m in VALID.finditer(sequence.upper())]


def window_starts(length: int, window: int, stride: int) -> list[tuple[int, int]]:
    """Windows over [0, length); the last one is the first to reach the end."""
    if length <= 0 or window <= 0 or stride <= 0 or stride > window:
        raise ValueError('Need positive length/window/stride with stride <= window.')
    out, start = [], 0
    while True:
        end = min(start + window, length)
        out.append((start, end))
        if end == length:
            return out
        start += stride


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load_inputs(args) -> list[tuple[str, Path]]:
    if args.manifest:
        with args.manifest.open() as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        inputs = [(r['genome_id'], Path(r['path'])) for r in rows]
    else:
        inputs = []
        for p in args.fasta:
            stem = p.name
            for suffix in ('.gz', '.fna', '.fa', '.fasta'):
                stem = stem[:-len(suffix)] if stem.endswith(suffix) else stem
            inputs.append((stem, p))
    ids = [g for g, _ in inputs]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicated genome_id.')
    return inputs


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--manifest', type=Path, help='TSV with genome_id and path columns')
    src.add_argument('--fasta', type=Path, nargs='+', help='FASTA(.gz); genome_id = file stem')
    p.add_argument('--out', type=Path, required=True, help='new output directory')
    p.add_argument('--window', type=int, default=3000)
    p.add_argument('--stride', type=int, default=1500)
    p.add_argument('--min-segment', type=int, default=40)
    a = p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    a.out.mkdir(parents=True, exist_ok=True)

    with (a.out / 'windows.fa').open('w') as fa, \
         (a.out / 'windows.tsv').open('w', newline='') as wt, \
         (a.out / 'segments.tsv').open('w', newline='') as st, \
         (a.out / 'genomes.tsv').open('w', newline='') as gt:
        ww = csv.DictWriter(wt, WINDOW_FIELDS, delimiter='\t')
        sw = csv.DictWriter(st, SEGMENT_FIELDS, delimiter='\t')
        gw = csv.DictWriter(gt, GENOME_FIELDS, delimiter='\t')
        for w in (ww, sw, gw):
            w.writeheader()
        for genome_id, path in load_inputs(a):
            g = dict.fromkeys(GENOME_FIELDS, 0)
            g.update(genome_id=genome_id, path=str(path), source_sha256=sha256_file(path))
            for contig, seq in read_fasta(path):
                g['n_contigs'] += 1
                g['total_nt'] += len(seq)
                for s0, e0 in valid_segments(seq):
                    seg_len = e0 - s0
                    g['acgt_nt'] += seg_len
                    g['n_segments'] += 1
                    status = 'windowed' if seg_len >= a.min_segment else 'short_unscored'
                    sw.writerow(dict(genome_id=genome_id, contig=contig, start0=s0,
                                     end0=e0, length_nt=seg_len, status=status))
                    if status != 'windowed':
                        g['n_short_segments'] += 1
                        g['short_segment_nt'] += seg_len
                        continue
                    g['windowed_nt'] += seg_len
                    for ws, we in window_starts(seg_len, a.window, a.stride):
                        start0, end0 = s0 + ws, s0 + we
                        w = seq[start0:end0].lower()
                        edges = []
                        if start0 == 0 or end0 == len(seq):
                            edges.append('contig_edge')
                        if (start0 == s0 and s0 > 0) or (end0 == e0 and e0 < len(seq)):
                            edges.append('non_acgt_edge')
                        if end0 - start0 < a.window:
                            edges.append('short_window')
                        wid = f'{genome_id}:{contig}:{start0}-{end0}'
                        fa.write(f'>{wid}\n{w}\n')
                        ww.writerow(dict(
                            window_id=wid, genome_id=genome_id, contig=contig,
                            start0=start0, end0=end0, length_nt=end0 - start0,
                            segment_start0=s0, segment_end0=e0, contig_length=len(seq),
                            edge_flag=','.join(edges) or '.',
                            gc_fraction=f'{(w.count("g") + w.count("c")) / len(w):.4f}',
                            window_sha256=hashlib.sha256(w.encode()).hexdigest()))
                        g['n_windows'] += 1
            g['non_acgt_nt'] = g['total_nt'] - g['acgt_nt']
            gw.writerow(g)
            print(f"{genome_id}: {g['n_contigs']} contigs, {g['total_nt']} nt, "
                  f"{g['n_windows']} windows, non-ACGT {g['non_acgt_nt']} nt, "
                  f"short-unscored {g['short_segment_nt']} nt", flush=True)


if __name__ == '__main__':
    main()
