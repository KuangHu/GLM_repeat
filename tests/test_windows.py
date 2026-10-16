"""CPU tests for genome windowing (make_windows.py)."""
import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

path = Path(__file__).parents[1] / 'scripts' / 'make_windows.py'
spec = importlib.util.spec_from_file_location('make_windows', path)
mw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mw)


@pytest.mark.parametrize('length,expected', [
    (4500, [(0, 3000), (1500, 4500)]),
    (5000, [(0, 3000), (1500, 4500), (3000, 5000)]),
    (3000, [(0, 3000)]),
    (3001, [(0, 3000), (1500, 3001)]),
    (40, [(0, 40)]),
    (6000, [(0, 3000), (1500, 4500), (3000, 6000)]),
])
def test_window_tail_rule(length, expected):
    assert mw.window_starts(length, 3000, 1500) == expected


def test_full_coverage():
    for length in range(1, 9000, 37):
        wins = mw.window_starts(length, 3000, 1500)
        covered = set()
        for s, e in wins:
            covered.update(range(s, e))
        assert covered == set(range(length))
        assert wins[-1][1] == length


def test_segments_split_at_non_acgt_and_keep_coordinates():
    assert mw.valid_segments('ACGTNNacgtRAC') == [(0, 4), (6, 10), (11, 13)]


def test_cli(tmp_path):
    seq1 = 'A' * 50 + 'N' * 3 + 'C' * 4000 + 'N' + 'G' * 10
    fa = tmp_path / 'g1.fa'
    fa.write_text(f'>c1 desc\n{seq1[:70]}\n{seq1[70:]}\n>c2\nACGT\n')
    out = tmp_path / 'out'
    subprocess.run([sys.executable, str(path), '--fasta', str(fa), '--out', str(out)],
                   check=True, capture_output=True)
    rows = list(csv.DictReader((out / 'windows.tsv').open(), delimiter='\t'))
    coords = [(r['contig'], int(r['start0']), int(r['end0'])) for r in rows]
    assert coords == [('c1', 0, 50), ('c1', 53, 3053), ('c1', 1553, 4053)]
    fasta = (out / 'windows.fa').read_text().split('\n')
    for r, (hdr, s) in zip(rows, zip(fasta[0::2], fasta[1::2])):
        assert hdr == '>' + r['window_id']
        assert s == seq1[int(r['start0']):int(r['end0'])].lower()
    segs = list(csv.DictReader((out / 'segments.tsv').open(), delimiter='\t'))
    assert [(s['contig'], s['status']) for s in segs] == [
        ('c1', 'windowed'), ('c1', 'windowed'), ('c1', 'short_unscored'),
        ('c2', 'short_unscored')]
    g = next(csv.DictReader((out / 'genomes.tsv').open(), delimiter='\t'))
    assert int(g['non_acgt_nt']) == 4 and int(g['short_segment_nt']) == 14
    assert int(g['acgt_nt']) == int(g['windowed_nt']) + int(g['short_segment_nt'])
