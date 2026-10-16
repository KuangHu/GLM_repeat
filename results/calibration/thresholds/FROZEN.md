# Frozen operating thresholds

Derived from calibration background ONLY (217 windows, 651 kb, disjoint from seed19
and from the validation panel). No positive, no seed19 window and no validation
window was read during calibration. Peak/cluster parameters are unchanged from the
original frozen design.

Freeze order, auditable from job records:
- RP thresholds frozen 2026-10-15 from job 26969803 (calibration RP), BEFORE the
  validation RP job 26969863 produced any track.
- Exact thresholds frozen from job 26969804 (calibration exact), before validation
  exact tracks were evaluated.

Target densities are the spec's operating points (WORKFLOW_SETUP.md section 8):
1, 10 and 100 merged background candidates per Mb.

## Resolution limit

With 0.651 Mb of background, one candidate already equals ~1.5 per Mb, so the
1/Mb operating point is below the resolution of this calibration set. In RP mode it
is additionally unreachable within the search grid: density is still 7.7/Mb at the
grid maximum tau = 1.35, and the grid is deliberately capped below the measurement
ceiling ln(4) = 1.386. The 1/Mb row is therefore reported as NOT ESTABLISHED rather
than assigned a threshold.

Thresholds in this directory are not to be re-derived against positives, and were
not re-tuned after any evaluation in this round.
0f3553bedf5be21ce68fd1eecdb8166114bc3164ea039ccd26f795e607f3e246  thresholds_rp_g20_r2.json
62d210701b3744a5953d604146fff12e7d13a3141cc39a2f70ab80c36d57d093  results/calibration/thresholds/thresholds_exact.json

## Exact-mode resolution (added after exact calibration)

The exact calibration background is 42 windows = 0.126 Mb, because exact scoring
costs ~75x RP. At that volume a single merged candidate already equals 7.9 per Mb,
so the 1/Mb and 10/Mb operating points are BELOW RESOLUTION and are reported as
NOT ESTABLISHED, not assigned a threshold. Only 100/Mb (tau = 0.92, achieved 87.3
per Mb, 11 candidates) is established for exact mode.

Establishing a 10/Mb exact operating point with >=10 candidates would need ~1 Mb of
exact background (~333 windows, ~150 GPU-hours on one A40). RP reaches the same
operating points at 1/75 of that cost, which is the practical argument for using RP
for calibration once RP-vs-exact equivalence has been demonstrated on the scorer in
use.
