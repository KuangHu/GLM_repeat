# Copy-dependence diagnostic

Target unit and its local flank are byte-identical in every arm; only the unchanged target positions are re-scored (exact single-position masking, both orientations).
Model 08754cba59, dtype torch.bfloat16, policy always. Seeds [11, 23, 37, 41, 53]. Mutation: every base in a mutated span replaced by a uniformly chosen different base.

DIAGNOSTIC ONLY: unit definitions are annotation-derived and are never inputs to the discovery scorer.

## Per target (mean over seeds, nats above uniform)

| locus | target | other copies | orient | original | copy | noncopy | copy-original | noncopy-original | copy-dependence |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| LMD9_CRISPR_1377229 | 1 | 8 | + | +1.385 | +0.026 | +1.386 | -1.360 | +0.000 | **-1.360** |
| LMD9_CRISPR_1377229 | 1 | 8 | - | +1.386 | +0.024 | +1.386 | -1.362 | +0.000 | **-1.362** |
| LMD9_CRISPR_1377229 | 4 | 8 | + | +1.386 | +0.041 | +1.386 | -1.345 | -0.000 | **-1.345** |
| LMD9_CRISPR_1377229 | 4 | 8 | - | +1.386 | +0.029 | +1.386 | -1.357 | +0.000 | **-1.357** |
| LMD9_CRISPR_1377229 | 7 | 8 | + | +1.385 | +0.039 | +1.385 | -1.346 | -0.000 | **-1.346** |
| LMD9_CRISPR_1377229 | 7 | 8 | - | +1.386 | +0.021 | +1.386 | -1.365 | +0.000 | **-1.365** |
| MarsHill_ART | 1 | 4 | + | +1.282 | +0.024 | +1.298 | -1.258 | +0.016 | **-1.274** |
| MarsHill_ART | 1 | 4 | - | +1.342 | +0.045 | +1.338 | -1.297 | -0.004 | **-1.293** |
| MarsHill_ART | 2 | 4 | + | +1.332 | -0.065 | +1.326 | -1.398 | -0.006 | **-1.391** |
| MarsHill_ART | 2 | 4 | - | +1.296 | -0.068 | +1.280 | -1.363 | -0.016 | **-1.348** |
| MarsHill_ART | 3 | 4 | + | +1.300 | -0.068 | +1.313 | -1.368 | +0.012 | **-1.380** |
| MarsHill_ART | 3 | 4 | - | +1.289 | -0.064 | +1.298 | -1.353 | +0.008 | **-1.361** |
| SA1_ART | 1 | 4 | + | +0.955 | +0.208 | +1.004 | -0.746 | +0.050 | **-0.796** |
| SA1_ART | 1 | 4 | - | +1.073 | +0.240 | +1.115 | -0.833 | +0.042 | **-0.875** |
| SA1_ART | 2 | 4 | + | +1.308 | -0.039 | +1.282 | -1.347 | -0.026 | **-1.321** |
| SA1_ART | 2 | 4 | - | +1.319 | -0.028 | +1.305 | -1.347 | -0.014 | **-1.333** |
| SA1_ART | 3 | 4 | + | +1.287 | -0.088 | +1.290 | -1.375 | +0.003 | **-1.379** |
| SA1_ART | 3 | 4 | - | +1.310 | -0.053 | +1.312 | -1.364 | +0.001 | **-1.365** |
| SUSP1_VIPR_YNNGG | 1 | 8 | + | +0.659 | +0.742 | +0.654 | +0.083 | -0.005 | **+0.088** |
| SUSP1_VIPR_YNNGG | 1 | 8 | - | +0.741 | +0.916 | +0.727 | +0.175 | -0.013 | **+0.189** |
| SUSP1_VIPR_YNNGG | 5 | 7 | + | +0.376 | +0.359 | +0.440 | -0.017 | +0.064 | **-0.081** |
| SUSP1_VIPR_YNNGG | 5 | 7 | - | +0.229 | +0.251 | +0.272 | +0.023 | +0.044 | **-0.021** |
| SUSP1_VIPR_YNNGG | 10 | 8 | + | +1.082 | +0.648 | +1.061 | -0.435 | -0.021 | **-0.413** |
| SUSP1_VIPR_YNNGG | 10 | 8 | - | +1.219 | +1.061 | +1.210 | -0.158 | -0.009 | **-0.149** |
| SUSP1_tRNA_cluster | 1 | 9 | + | +0.808 | +0.812 | +0.814 | +0.004 | +0.006 | **-0.002** |
| SUSP1_tRNA_cluster | 1 | 9 | - | +0.778 | +0.801 | +0.784 | +0.023 | +0.006 | **+0.017** |
| SUSP1_tRNA_cluster | 5 | 10 | + | +0.825 | +0.776 | +0.829 | -0.048 | +0.004 | **-0.053** |
| SUSP1_tRNA_cluster | 5 | 10 | - | +0.836 | +0.817 | +0.829 | -0.018 | -0.007 | **-0.012** |
| SUSP1_tRNA_cluster | 10 | 10 | + | +0.930 | +0.922 | +0.931 | -0.009 | +0.001 | **-0.010** |
| SUSP1_tRNA_cluster | 10 | 10 | - | +0.956 | +0.946 | +0.960 | -0.010 | +0.005 | **-0.015** |

## Per locus (mean over targets and orientations)

| locus | targets x orient | copy-original | noncopy-original | copy-dependence | seed SD (copy) | seed SD (noncopy) |
|---|---:|---:|---:|---:|---:|---:|
| LMD9_CRISPR_1377229 | 6 | -1.356 | +0.000 | **-1.356** | 0.012 | 0.000 |
| MarsHill_ART | 6 | -1.339 | +0.002 | **-1.341** | 0.019 | 0.018 |
| SA1_ART | 6 | -1.169 | +0.009 | **-1.178** | 0.016 | 0.041 |
| SUSP1_VIPR_YNNGG | 6 | -0.055 | +0.010 | **-0.064** | 0.097 | 0.017 |
| SUSP1_tRNA_cluster | 6 | -0.010 | +0.003 | **-0.012** | 0.004 | 0.004 |

A strongly negative copy-dependence means the target unit is predictable largely *because* its other copies are present. A value near zero means the target is predictable on its own, or equally sensitive to any perturbation.

Caveats: the non-copy arm matches perturbed-nt budget and approximate distance to the target, not local composition or annotation status. Distance matching widens when placement is tight (see `distance_mismatch_max` in scores.tsv). Target counts are small and these loci are not independent (SA1/MarsHill share a lineage).

## Reading of this run (2026-10-15)

Three loci show large, seed-stable copy dependence: destroying the other copies
costs the target 1.17-1.36 nats, while destroying the same number of nucleotides of
matched non-copy sequence costs it nothing (+0.000 to +0.009). Seed SD is 0.004-0.041,
so the contrast is far outside seed noise.

The tRNA cluster is a clean null: -0.012 overall, consistent across all three
targets and both orientations (range -0.053 to +0.017). Its units are individually
predictable; they do not need each other. This separates it from the repeat arrays
on a mechanism rather than on score magnitude, and unlike the peak-cluster zero for
w0011 it is not an artefact of a threshold.

**The VIPR arm is not interpretable at this unit size.** Its three targets disagree
in sign (+0.189, -0.021, -0.413) with the largest seed SD in the panel (0.097). The
units are 5 nt with a 10-nt flank, so "other copies" include the pentamers directly
abutting the preserved zone: the copy arm necessarily rewrites the target's
immediate context, which the non-copy arm does not. A larger unit definition (the
whole 60-nt periodic segment, or the 94-nt vrRNA against other vrRNAs) would be
needed to test VIPR copy dependence properly.

Distance matching is good for the ART loci and VIPR (<=150 nt) but degrades where
placement is tight: CRISPR 485 nt, tRNA cluster 600 nt. For the tRNA cluster this
weakens the control, though its near-zero effect in *both* arms makes a
distance artefact an unlikely explanation.
