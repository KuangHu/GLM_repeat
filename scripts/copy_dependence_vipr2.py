#!/usr/bin/env python3
"""VIPR copy-dependence v2: the 60-nt YNNGG region as the primary test region.

v1 was uninterpretable because, with 5-nt units and a 10-nt flank, the perturbed
units abutted the preserved zone, so the 'copy' arm rewrote the target's immediate
context. v2 preserves the target unit(s) AND an immediate flanking context of
`FLANK` nt, and perturbs only units lying entirely outside that zone - i.e. at
least 3 periods away.

Two unit definitions over the same region (window w0003, [1487,1547)):
  single : 12 contiguous 5-nt YNNGG units, target = one unit
  pair   : 6 contiguous 10-nt blocks (two units each), target = one block

The 94-nt vrRNA test is NOT run: see VIPR_V2_NOTES.md. Everything else (mutation,
budget- and distance-matched non-copy arm, fixed seeds, scoring of unchanged
target positions only) is reused unchanged from copy_dependence.py.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import copy_dependence as cd  # noqa: E402

REGION = (1487, 1547)   # window-local coords of the mapped 60-nt YNNGG segment
FLANK = 15              # >= 3 periods: immediate flanking context is preserved
TARGETS = {'single': [0, 5, 10], 'pair': [0, 2, 5]}


def vipr_units(bench, seqs):
    r0, r1 = REGION
    single = [(r0 + 5 * k, r0 + 5 * k + 5) for k in range((r1 - r0) // 5)]
    pair = [(r0 + 10 * k, r0 + 10 * k + 10) for k in range((r1 - r0) // 10)]
    return {
        'SUSP1_VIPR_60nt_single': (
            'w0003', FLANK, single,
            f'12 contiguous 5-nt YNNGG units in [{r0},{r1}); flank {FLANK} nt so only '
            'units >= 3 periods away are perturbed'),
        'SUSP1_VIPR_60nt_pair': (
            'w0003', FLANK, pair,
            f'6 contiguous 10-nt blocks (2 YNNGG units each) in [{r0},{r1}); '
            f'flank {FLANK} nt'),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bench', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    # locus-specific target choice: ends and middle, so distant units exist on at
    # least one side of every target.
    cd.locus_units = vipr_units
    cd.pick_targets = lambda units: TARGETS['single' if len(units) == 12 else 'pair']
    cd.build(a)


if __name__ == '__main__':
    main()
