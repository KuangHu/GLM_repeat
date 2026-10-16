#!/usr/bin/env python3
"""Write synthetic engineering fixtures, explicitly NOT biological benchmarks."""
import argparse
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(20261014)
    dna = lambda n: ''.join(rng.choice('acgt') for _ in range(n))
    motif = dna(24)
    spaced = motif + ''.join(dna(36) + motif for _ in range(5))
    core = dna(18)
    long_spaced = core + ''.join(dna(150) + core for _ in range(4))
    short_periodic = ''.join(rng.choice('ct') + dna(2) + 'gg' for _ in range(14))
    examples = {
        'toy_short_partial_periodicity': dna(120) + short_periodic + dna(120),
        'toy_spaced_repeat': dna(120) + spaced + dna(120),
        'toy_long_spaced_repeat': dna(120) + long_spaced + dna(120),
        'toy_simple_repeat_positive': dna(120) + 'ac' * 40 + dna(120),
        'toy_homopolymer_positive': dna(120) + 'a' * 80 + dna(120),
        'toy_iid_background_not_a_biological_negative': dna(700),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w') as fh:
        for name, seq in examples.items():
            fh.write('>' + name + '\n' + seq + '\n')


if __name__ == '__main__':
    main()
