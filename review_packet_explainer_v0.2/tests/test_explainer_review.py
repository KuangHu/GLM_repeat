"""Checks that reproduce the review findings in review/FINDINGS.md.

These are regression tests for the REVIEWED behaviour of v0.2. Several of them
assert that a known defect is still present; when a defect is fixed in a later
version the corresponding test must be updated, which makes the fix visible.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    'ex2', Path(__file__).parents[1] / 'code' / 'explain_repeats_v2.py')
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


def test_B2_phase_information_is_confounded_by_unit_count():
    """All phases of a perfectly periodic sequence are circularly equivalent, yet
    the phase fitting one extra unit scores highest."""
    seq = 'ACGTT' * 20
    scores = {}
    for ph in range(5):
        units = [seq[ph + i * 5: ph + (i + 1) * 5] for i in range((len(seq) - ph) // 5)]
        units = [u for u in units if len(u) == 5]
        scores[ph] = (len(units), ex.information(ex.pwm_from(units)))
    assert scores[0][0] == 20 and all(scores[p][0] == 19 for p in (1, 2, 3, 4))
    assert scores[0][1] > scores[1][1]                      # defect present
    assert scores[0][1] - scores[1][1] < 0.02               # and small


def test_L1_self_support_is_inflated_on_random_sequence():
    rng = np.random.default_rng(0)
    for n, floor in ((2, 0.55), (3, 0.50)):
        vals = []
        for _ in range(200):
            copies = [''.join(rng.choice(list('ACGT'), 5)) for _ in range(n)]
            cons = ex.consensus(ex.pwm_from(copies))
            vals.append(np.mean([1 - ex.mismatches(c, cons) / 5 for c in copies]))
        assert np.mean(vals) > floor      # noise manufactures support


def test_L1_leave_one_out_support_is_at_chance():
    rng = np.random.default_rng(1)
    vals = []
    for _ in range(200):
        copies = [''.join(rng.choice(list('ACGT'), 5)) for _ in range(5)]
        for i, c in enumerate(copies):
            others = copies[:i] + copies[i + 1:]
            vals.append(1 - ex.mismatches(c, ex.consensus(ex.pwm_from(others))) / 5)
    assert 0.20 < np.mean(vals) < 0.32


def test_L2_information_increases_with_copy_count_alone():
    unit = 'ACGTT'
    vals = [ex.information(ex.pwm_from([unit] * n)) for n in (3, 10, 50)]
    assert vals[0] < vals[1] < vals[2]


def test_L3_selection_score_favours_long_periods_on_noise():
    rng = np.random.default_rng(2)
    def score(d, n):
        c = [''.join(rng.choice(list('ACGT'), d)) for _ in range(n)]
        return n * d * ex.information(ex.pwm_from(c))
    assert score(30, 2) > 3 * score(2, 30)


@pytest.mark.parametrize('seq,period', [('ACGTT' * 12, 5), ('AT' * 40, 2)])
def test_periodic_branch_finds_a_multiple_of_the_planted_period(seq, period):
    """Exact-repeat case: the returned offset is a MULTIPLE of the true period,
    not necessarily the period itself (B5 harmonics)."""
    d, ident = ex.best_period(seq)
    assert d % period == 0 and ident > 0.9


def test_B5_period_search_favours_large_offsets_on_random_sequence():
    rng = np.random.default_rng(5)
    picks = [ex.best_period(''.join(rng.choice(list('ACGT'), 1000)))[0]
             for _ in range(40)]
    assert np.mean([d > 100 for d in picks]) > 0.5      # defect present


def test_periodic_branch_needs_no_seed_kmer():
    """A degenerate period-5 pattern has no repeated 12-mer, yet is still modelled."""
    rng = np.random.default_rng(3)
    seq = ''.join('C' + ''.join(rng.choice(list('ACGT'), 2)) + 'GG' for _ in range(60))
    from collections import Counter
    # the degenerate pattern provides no strongly repeated seed k-mer
    top = Counter(seq[i:i + 12] for i in range(len(seq) - 11)).most_common(1)[0][1]
    assert top <= 2, f'unexpectedly repeated 12-mer x{top}'
    d, ident = ex.best_period(seq)
    # B5 (defect, documented): the search is dragged to large offsets because
    # identity at offset d uses only len-d comparisons. A clean period-5 pattern
    # does NOT return 5. Update this assertion when B5 is fixed.
    assert ident > 0.5
    assert d != 5, 'B5 appears fixed - update this test and FINDINGS.md'
    assert ex.period_identity(seq, 5)[0] > 0.5, 'true period still has high identity' 


def test_core_copies_are_preserved_when_extension_fails():
    models = ex.motif_models('AAACATTATATGAATTTT' + 'G' * 40 +
                             'CCCCATTATATGAATGGG' + 'T' * 40 +
                             'TTTCATTATATGAATCCC', min_copies=2)
    assert models
    m = models[0]
    assert m['core_copies'] >= 3
    assert m['source'] in ('motif_core', 'motif_extended')
    if m['source'] == 'motif_core':
        assert m['core_copies'] >= (m['extended_copies'] or 0)


def test_no_known_motif_is_hard_coded():
    src = (Path(__file__).parents[1] / 'code' / 'explain_repeats_v2.py').read_text()
    for forbidden in ('ATATGAATACGTAT', 'YNNGG', 'GTTCCAATAAGA', '29169', '29229'):
        assert forbidden not in src
