# Copy-dependence diagnostic

Target unit and its local flank are byte-identical in every arm; only the unchanged target positions are re-scored (exact single-position masking, both orientations).
Model 08754cba59, dtype torch.bfloat16, policy always. Seeds [11, 23, 37, 41, 53]. Mutation: every base in a mutated span replaced by a uniformly chosen different base.

DIAGNOSTIC ONLY: unit definitions are annotation-derived and are never inputs to the discovery scorer.

## Per target (mean over seeds, nats above uniform)

| locus | target | other copies | orient | original | copy | noncopy | copy-original | noncopy-original | copy-dependence |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| SUSP1_VIPR_60nt_pair | 0 | 3 | + | +0.836 | +0.823 | +0.807 | -0.013 | -0.029 | **+0.016** |
| SUSP1_VIPR_60nt_pair | 0 | 3 | - | +0.814 | +0.807 | +0.801 | -0.008 | -0.013 | **+0.006** |
| SUSP1_VIPR_60nt_pair | 2 | 1 | + | +0.558 | +0.598 | +0.562 | +0.040 | +0.005 | **+0.035** |
| SUSP1_VIPR_60nt_pair | 2 | 1 | - | +0.634 | +0.642 | +0.640 | +0.007 | +0.006 | **+0.002** |
| SUSP1_VIPR_60nt_pair | 5 | 3 | + | +0.666 | +0.641 | +0.672 | -0.026 | +0.006 | **-0.032** |
| SUSP1_VIPR_60nt_pair | 5 | 3 | - | +0.855 | +0.872 | +0.863 | +0.017 | +0.008 | **+0.009** |
| SUSP1_VIPR_60nt_single | 0 | 8 | + | +1.012 | +0.987 | +0.949 | -0.025 | -0.063 | **+0.038** |
| SUSP1_VIPR_60nt_single | 0 | 8 | - | +0.886 | +0.750 | +0.868 | -0.136 | -0.018 | **-0.118** |
| SUSP1_VIPR_60nt_single | 5 | 5 | + | +0.376 | +0.436 | +0.391 | +0.060 | +0.015 | **+0.044** |
| SUSP1_VIPR_60nt_single | 5 | 5 | - | +0.229 | +0.259 | +0.246 | +0.030 | +0.017 | **+0.013** |
| SUSP1_VIPR_60nt_single | 10 | 7 | + | +1.082 | +0.861 | +1.066 | -0.222 | -0.016 | **-0.205** |
| SUSP1_VIPR_60nt_single | 10 | 7 | - | +1.219 | +1.194 | +1.212 | -0.026 | -0.007 | **-0.018** |

## Per locus (mean over targets and orientations)

| locus | targets x orient | copy-original | noncopy-original | copy-dependence | seed SD (copy) | seed SD (noncopy) |
|---|---:|---:|---:|---:|---:|---:|
| SUSP1_VIPR_60nt_single | 6 | -0.053 | -0.012 | **-0.041** | 0.041 | 0.023 |
| SUSP1_VIPR_60nt_pair | 6 | +0.003 | -0.003 | **+0.006** | 0.011 | 0.013 |

A strongly negative copy-dependence means the target unit is predictable largely *because* its other copies are present. A value near zero means the target is predictable on its own, or equally sensitive to any perturbation.

Caveats: the non-copy arm matches perturbed-nt budget and approximate distance to the target, not local composition or annotation status. Distance matching widens when placement is tight (see `distance_mismatch_max` in scores.tsv). Target counts are small and these loci are not independent (SA1/MarsHill share a lineage).

## Reading of VIPR v2 (2026-10-15)

With the v1 confound removed (nearest perturbed unit now 15-20 nt from the target
instead of 10 nt), the 60-nt YNNGG region shows **no array-like copy dependence**:
-0.041 for single 5-nt units and +0.006 for 10-nt pairs, against -1.17 to -1.36 for
the CRISPR and ART arrays measured under the same protocol. The `pair` definition is
the tighter of the two (all six rows within +/-0.035, seed SD 0.011-0.013).

Per-target signs remain inconsistent for the `single` definition (+0.044 to -0.205),
so even the small negative mean is not a stable effect.

This is not explained by a smaller perturbation budget: the ART arms perturbed only
56 nt and still produced -1.17, while VIPR single arms perturbed 25-40 nt.

Interpretation, stated conservatively: within this one locus, the predictability of
a YNNGG unit does not depend on the presence of the distant units of the same
periodic region. A degenerate 5-nt pattern is largely specified by local periodic
phase and composition, which the preserved flank retains. This says nothing about
whether the vrRNA is functional, and v1's apparent -0.413 on one target is best
explained by the immediate-context confound that v2 removes.

Caveat: the middle `pair` target had only one distant block (10 nt) available to
perturb, the smallest budget in the set.
