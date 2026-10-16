# Seed19, RP tracks: threshold policy comparison

Thresholds were calibrated on background only and frozen before evaluation (results/calibration/thresholds/FROZEN.md). The original `window_local` scorer and its report are preserved unchanged.

## Policy summary

| policy | threshold | merged candidates | positives non-zero | backgrounds non-zero | tracks with unreachable threshold |
|---|---:|---:|---:|---:|---:|
| window_local | per-window | 13 | 5/5 | 6/14 | 2 |
| bg10_per_Mb | 1.33 | 3 | 3/5 | 0/14 | 0 |
| bg100_per_Mb | 0.89 | 13 | 5/5 | 6/14 | 0 |

Unreachable thresholds (tau >= ln(4) = 1.3863). Such a window yields zero peaks **by construction**; this is never evidence that no repetitive organization exists:

- window_local: w0011 (+) tau = 1.4536
- window_local: w0011 (-) tau = 1.4871

## Window scores (all windows)

| record | family | window_local | bg10_per_Mb | bg100_per_Mb |
|---|---|---:|---:|---:|
| w0001 | CRISPR | 5.409 | 5.779 | 5.410 |
| w0002 | CRISPR | 4.705 | 4.602 | 4.638 |
| w0011 **<-** | background | 0.000 | 0.000 | 4.106 |
| w0005 | ART | 3.247 | 2.705 | 3.433 |
| w0003 | VIPR | 2.903 | 0.000 | 3.170 |
| w0018 | background | 2.302 | 0.000 | 2.390 |
| w0004 | ART | 2.320 | 0.000 | 2.307 |
| w0019 | background | 2.256 | 0.000 | 2.151 |
| w0006 | background | 2.068 | 0.000 | 1.788 |
| w0007 | background | 1.556 | 0.000 | 0.000 |
| w0012 | background | 1.503 | 0.000 | 1.532 |
| w0015 | background | 0.000 | 0.000 | 1.530 |
| w0008 | background | 1.242 | 0.000 | 0.000 |
| w0009 | background | 0.000 | 0.000 | 0.000 |
| w0010 | background | 0.000 | 0.000 | 0.000 |
| w0013 | background | 0.000 | 0.000 | 0.000 |
| w0014 | background | 0.000 | 0.000 | 0.000 |
| w0016 | background | 0.000 | 0.000 | 0.000 |
| w0017 | background | 0.000 | 0.000 | 0.000 |

## Rank of each positive (of 19 windows)

| record | family | window_local | bg10_per_Mb | bg100_per_Mb |
|---|---|---:|---:|---:|
| w0001 | CRISPR | 1 | 1 | 1 |
| w0002 | CRISPR | 2 | 2 | 2 |
| w0003 | VIPR | 4 | 4 | 5 |
| w0004 | ART | 5 | 5 | 7 |
| w0005 | ART | 3 | 3 | 4 |

## w0011 reported separately

Panel label: **background** (unchanged). Source SUSP1.

| policy | score | rank | non-zero |
|---|---:|---:|---|
| window_local | 0.000 | 14 | no |
| bg10_per_Mb | 0.000 | 11 | no |
| bg100_per_Mb | 4.106 | 3 | yes |

Backgrounds are unlabeled real DNA, not confirmed negatives. No threshold was tuned against positives and none was re-tuned after this evaluation.
