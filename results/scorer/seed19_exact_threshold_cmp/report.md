# Seed19, EXACT tracks: threshold policy comparison

Thresholds were calibrated on background only and frozen before evaluation (results/calibration/thresholds/FROZEN.md). The original `window_local` scorer and its report are preserved unchanged.

## Policy summary

| policy | threshold | merged candidates | positives non-zero | backgrounds non-zero | tracks with unreachable threshold |
|---|---:|---:|---:|---:|---:|
| window_local | per-window | 13 | 5/5 | 6/14 | 2 |
| bg100_per_Mb | 0.92 | 13 | 5/5 | 6/14 | 0 |

Unreachable thresholds (tau >= ln(4) = 1.3863). Such a window yields zero peaks **by construction**; this is never evidence that no repetitive organization exists:

- window_local: w0011 (+) tau = 1.5963
- window_local: w0011 (-) tau = 1.5962

## Window scores (all windows)

| record | family | window_local | bg100_per_Mb |
|---|---|---:|---:|
| w0001 | CRISPR | 5.465 | 5.454 |
| w0002 | CRISPR | 4.663 | 4.513 |
| w0011 **<-** | background | 0.000 | 4.382 |
| w0005 | ART | 3.362 | 3.499 |
| w0003 | VIPR | 2.930 | 3.006 |
| w0018 | background | 2.465 | 2.263 |
| w0019 | background | 2.363 | 2.153 |
| w0004 | ART | 2.332 | 2.332 |
| w0008 | background | 2.324 | 0.000 |
| w0006 | background | 2.119 | 2.272 |
| w0012 | background | 1.631 | 1.673 |
| w0007 | background | 1.661 | 0.000 |
| w0015 | background | 0.000 | 1.573 |
| w0009 | background | 0.000 | 0.000 |
| w0010 | background | 0.000 | 0.000 |
| w0013 | background | 0.000 | 0.000 |
| w0014 | background | 0.000 | 0.000 |
| w0016 | background | 0.000 | 0.000 |
| w0017 | background | 0.000 | 0.000 |

## Rank of each positive (of 19 windows)

| record | family | window_local | bg100_per_Mb |
|---|---|---:|---:|
| w0001 | CRISPR | 1 | 1 |
| w0002 | CRISPR | 2 | 2 |
| w0003 | VIPR | 4 | 5 |
| w0004 | ART | 7 | 6 |
| w0005 | ART | 3 | 4 |

## w0011 reported separately

Panel label: **background** (unchanged). Source SUSP1.

| policy | score | rank | non-zero |
|---|---:|---:|---|
| window_local | 0.000 | 14 | no |
| bg100_per_Mb | 4.382 | 3 | yes |

Backgrounds are unlabeled real DNA, not confirmed negatives. No threshold was tuned against positives and none was re-tuned after this evaluation.
