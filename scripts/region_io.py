"""Shared, strict coordinate I/O for the optional region-selection workflow."""
from __future__ import annotations

import csv
import hashlib
import json
from bisect import bisect_left, bisect_right
from pathlib import Path
from urllib.parse import unquote

from make_windows import read_fasta, sha256_file


def table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_table(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def new_dir(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'Output directory is nonempty: {path}')
    path.mkdir(parents=True, exist_ok=True)
    return path


def fasta(path):
    seqs = dict(read_fasta(Path(path)))
    if not seqs or any(not s for s in seqs.values()):
        raise ValueError('FASTA must have nonempty sequences with globally unique contig IDs')
    return seqs


def interval(row, seqs):
    c, s, e = row['contig'], int(row['start0']), int(row['end0'])
    if c not in seqs or not 0 <= s < e <= len(seqs[c]):
        raise ValueError(f'Invalid/unknown interval: {c}:{s}-{e}')
    return c, s, e


def union(spans):
    result = []
    for s, e in sorted(spans):
        if e <= s:
            continue
        if result and s <= result[-1][1]:
            result[-1] = (result[-1][0], max(e, result[-1][1]))
        else:
            result.append((s, e))
    return result


def complement(spans, start, end):
    cursor = start
    result = []
    for s, e in union((max(start, s), min(end, e)) for s, e in spans):
        if cursor < s:
            result.append((cursor, s))
        cursor = max(cursor, e)
    if cursor < end:
        result.append((cursor, end))
    return result


class SpanIndex:
    """Union index; interval queries never count overlapping bases twice."""
    def __init__(self, spans):
        self.spans = union(spans)
        self.starts = [s for s, _ in self.spans]
        self.ends = [e for _, e in self.spans]
        self.prefix = [0]
        for s, e in self.spans:
            self.prefix.append(self.prefix[-1] + e - s)

    def overlap(self, start, end):
        a = bisect_right(self.ends, start)
        b = bisect_left(self.starts, end)
        if a >= b:
            return 0
        return (self.prefix[b] - self.prefix[a] - max(0, start-self.starts[a])
                - max(0, self.ends[b-1]-end))


def read_gff(path, seqs):
    """Read GFF3 CDSs, preserving joined parts and strand. Never infer confidence.

    Repeated CDS IDs (or Parent if ID absent) denote parts of one gene.
    A discontinuous feature's intervening bases are NOT classified as coding.
    """
    genes = {}
    with Path(path).open(encoding='utf-8-sig') as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                raise ValueError(f'{path}:{line_no}: expected nine GFF3 columns')
            if f[2] != 'CDS':
                continue
            attrs = dict((unquote(k), unquote(v)) for k, v in
                         (item.split('=', 1) for item in f[8].split(';') if '=' in item))
            gid = attrs.get('ID') or attrs.get('Parent') or f'line{line_no}'
            c, s, e = interval(dict(contig=unquote(f[0]), start0=int(f[3])-1, end0=f[4]), seqs)
            if f[6] not in ('+', '-'):
                raise ValueError(f'CDS {gid}: strand must be + or -')
            key = (c, gid)
            g = genes.setdefault(key, dict(contig=c, gene_id=gid, strand=f[6], parts=[],
                attrs=attrs, product=attrs.get('product', ''), source=f[1]))
            if g['strand'] != f[6]:
                raise ValueError(f'Conflicting strands for CDS {key}')
            g['parts'].append((s, e))
    for g in genes.values():
        parts = sorted(g['parts'])
        if any(a[1] > b[0] for a, b in zip(parts, parts[1:])):
            raise ValueError(f'Overlapping parts of CDS {g["gene_id"]}')
        g.update(parts=parts, start0=parts[0][0], end0=parts[-1][1],
                 length_nt=sum(e-s for s, e in parts))
    return sorted(genes.values(), key=lambda g: (g['contig'], g['start0'], g['gene_id']))


def provenance(out, args, sources, **extra):
    data = dict(parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                input_sha256={str(p): sha256_file(Path(p)) for p in sources if p},
                code_sha256={p.name: sha256_file(p) for p in Path(__file__).parent.glob('*.py')},
                coordinates='0-based half-open; GFF3 converted from 1-based inclusive', **extra)
    (Path(out)/'provenance.json').write_text(json.dumps(data, indent=2, default=str)+'\n', encoding='utf-8')


def stable_id(prefix, text):
    return prefix + hashlib.sha256(text.encode()).hexdigest()[:16]
