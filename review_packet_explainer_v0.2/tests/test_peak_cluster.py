"""Peak-cluster scorer: peak calling, clustering and cross-window merging."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

path = Path(__file__).parents[1] / 'scripts' / 'peak_cluster.py'
spec = importlib.util.spec_from_file_location('peak_cluster', path)
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)


def track_with_peaks(length, spans, high=3.0, low=0.0):
    v = np.full(length, low)
    for s, e in spans:
        v[s:e] = high
    return v


def test_smooth_preserves_length():
    for n in [1, 5, 15, 16, 17, 100]:
        assert len(pc.smooth(np.arange(n, dtype=float))) == n


def test_flat_track_has_no_peaks():
    peaks, *_ = pc.call_peaks(np.zeros(1000))
    assert peaks == []


def test_peaks_found_at_expected_places():
    v = track_with_peaks(1000, [(100, 140), (300, 340), (500, 540)])
    peaks, tau, med, sd = pc.call_peaks(v)
    assert len(peaks) == 3
    for pk, (s, e) in zip(peaks, [(100, 140), (300, 340), (500, 540)]):
        assert abs(pk['start'] - s) <= pc.SMOOTH_NT and abs(pk['end'] - e) <= pc.SMOOTH_NT


def test_short_peaks_rejected():
    v = track_with_peaks(1000, [(100, 103), (300, 303)])
    peaks, *_ = pc.call_peaks(v)
    assert all(p['end'] - p['start'] >= pc.MIN_PEAK_NT for p in peaks)


def test_single_peak_is_not_a_cluster():
    peaks = [{'start': 100, 'end': 140, 'mean': 2.0, 'summit': 120}]
    assert pc.form_clusters(peaks) == []


def test_two_near_peaks_cluster_far_peaks_do_not():
    near = [{'start': 0, 'end': 20, 'mean': 2.0, 'summit': 10},
            {'start': 200, 'end': 220, 'mean': 2.0, 'summit': 210}]
    assert len(pc.form_clusters(near)) == 1
    far = [{'start': 0, 'end': 20, 'mean': 2.0, 'summit': 10},
           {'start': 5000, 'end': 5020, 'mean': 2.0, 'summit': 5010}]
    assert pc.form_clusters(far) == []


def test_cluster_score_rewards_more_peaks():
    mk = lambda n: [{'start': i * 100, 'end': i * 100 + 20, 'mean': 2.0,
                     'summit': i * 100 + 10} for i in range(n)]
    s2 = pc.form_clusters(mk(2))[0]['score']
    s5 = pc.form_clusters(mk(5))[0]['score']
    assert s5 > s2 == pytest.approx(2.0 * np.log2(3))


def test_cluster_score_uses_median_not_max():
    g = [{'start': 0, 'end': 20, 'mean': 1.0, 'summit': 5},
         {'start': 100, 'end': 120, 'mean': 1.0, 'summit': 105},
         {'start': 200, 'end': 220, 'mean': 9.0, 'summit': 205}]
    assert pc.form_clusters(g)[0]['median_peak_mean'] == 1.0


def test_overlapping_windows_merge_into_one_candidate():
    base = {'contig': 'c1', 'n_peaks': 3, 'median_peak_mean': 1.0,
            'orientation': '+', 'peak_lengths': '', 'gaps': ''}
    cands = [{**base, 'window_id': 'w1', 'start0': 1000, 'end0': 1500, 'score': 2.0},
             {**base, 'window_id': 'w2', 'start0': 1400, 'end0': 1900, 'score': 3.0}]
    merged = pc.merge_genomic(cands)
    assert len(merged) == 1
    assert (merged[0]['start0'], merged[0]['end0']) == (1000, 1900)
    assert merged[0]['score'] == 3.0 and merged[0]['n_windows'] == 2
    assert merged[0]['best_window'] == 'w2'


def test_non_overlapping_stay_separate_and_contigs_isolated():
    base = {'n_peaks': 2, 'median_peak_mean': 1.0, 'orientation': '+',
            'peak_lengths': '', 'gaps': '', 'score': 1.0}
    cands = [{**base, 'contig': 'c1', 'window_id': 'w1', 'start0': 0, 'end0': 100},
             {**base, 'contig': 'c1', 'window_id': 'w2', 'start0': 500, 'end0': 600},
             {**base, 'contig': 'c2', 'window_id': 'w3', 'start0': 50, 'end0': 150}]
    assert len(pc.merge_genomic(cands)) == 3


def test_frozen_parameters_unchanged():
    assert (pc.SMOOTH_NT, pc.THRESH_K, pc.MIN_PEAK_NT, pc.MAX_GAP_NT, pc.MIN_PEAKS) \
        == (16, 3.0, 8, 500, 2)


def test_window_local_policy_is_unchanged_default():
    v = track_with_peaks(1000, [(100, 140), (300, 340)])
    a, *_ = pc.call_peaks(v)
    b, *_ = pc.call_peaks(v, 'window_local')
    assert a == b


def test_background_calibrated_uses_fixed_threshold():
    v = track_with_peaks(1000, [(100, 140), (300, 340)], high=1.0, low=0.0)
    few, tau_hi, *_ = pc.call_peaks(v, 'background_calibrated', 0.9)
    many, tau_lo, *_ = pc.call_peaks(v, 'background_calibrated', 0.1)
    assert tau_hi == 0.9 and tau_lo == 0.1
    assert len(many) >= len(few)


def test_background_calibrated_requires_threshold():
    with pytest.raises(ValueError):
        pc.call_peaks(np.zeros(100), 'background_calibrated')
    with pytest.raises(ValueError):
        pc.call_peaks(np.zeros(100), 'nonsense')


def test_unreachable_threshold_detection():
    assert pc.threshold_is_unreachable(1.3863) is True
    assert pc.threshold_is_unreachable(1.0) is False
    # the w0011 case that produced a zero by construction
    assert pc.threshold_is_unreachable(1.596) is True


def test_fixed_threshold_below_ceiling_can_fire():
    v = np.full(500, 1.3)
    peaks, tau, *_ = pc.call_peaks(v, 'background_calibrated', 1.0)
    assert not pc.threshold_is_unreachable(tau) and len(peaks) == 1


# --- pilot scan: unique scanned nt must not double-count stride overlap ---
import importlib.util as _ilu
_ps = Path(__file__).parents[1] / 'scripts' / 'pilot_scan.py'
_spec = _ilu.spec_from_file_location('pilot_scan', _ps)
ps = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(ps)


def _win(contig, s, e):
    return {'contig': contig, 'start0': str(s), 'end0': str(e)}


def test_unique_scanned_nt_merges_overlapping_windows():
    w = {'a': _win('c1', 0, 3000), 'b': _win('c1', 1500, 4500), 'c': _win('c1', 3000, 6000)}
    total, per = ps.unique_scanned_nt(w, list(w))
    assert total == 6000 and per['c1'] == 6000     # not 9000


def test_unique_scanned_nt_separates_contigs_and_gaps():
    w = {'a': _win('c1', 0, 3000), 'b': _win('c2', 0, 3000),
         'c': _win('c1', 10000, 11000)}
    total, per = ps.unique_scanned_nt(w, list(w))
    assert per == {'c1': 4000, 'c2': 3000} and total == 7000


def test_unique_scanned_nt_excludes_unscored_windows():
    w = {'a': _win('c1', 0, 3000), 'b': _win('c1', 1500, 4500)}
    total, _ = ps.unique_scanned_nt(w, ['a'])       # 'b' has no track
    assert total == 3000
