# Peak-cluster scorer on the seed panel

Scorer design is agent-authored and was frozen before any track was scored: results/scorer/PEAK_CLUSTER_DESIGN.md.
Parameters: {"SMOOTH_NT": 16, "THRESH_K": 3.0, "MIN_PEAK_NT": 8, "MAX_GAP_NT": 500, "MIN_PEAKS": 2, "cluster_score": "median(peak_mean) * log2(1 + n_peaks)", "window_score": "max over clusters and orientations"}

Tracks: 38 (exact mode, model 08754cba59, policy always).
seed19 is development data: no FDR, no calibrated threshold, no recall claim.

## 1. Candidate coordinates (merged across windows, genomic, 0-based half-open)

| contig | start0 | end0 | span | score | peaks | median peak mean | orient | windows |
|---|---:|---:|---:|---:|---:|---:|---|---|
| NC_008532.1 | 649126 | 650211 | 1085 | 5.4652 | 17 | 1.3106 | - | w0001 |
| NC_008532.1 | 1376866 | 1378094 | 1226 | 4.6633 | 12 | 1.2602 | + | w0002 |
| MW248466.1 | 213994 | 214785 | 791 | 3.362 | 5 | 1.3006 | - | w0005 |
| NC_028808.2 | 27682 | 28002 | 319 | 2.9303 | 5 | 1.1336 | - | w0003 |
| NC_000964.3 | 624359 | 625143 | 784 | 2.4649 | 4 | 1.0616 | + | w0018 |
| NC_000964.3 | 2817867 | 2818034 | 167 | 2.3634 | 5 | 0.9143 | + | w0019 |
| MW218148.1 | 9951 | 10370 | 419 | 2.3321 | 3 | 1.1661 | - | w0004 |
| NC_008532.1 | 745310 | 746042 | 732 | 2.3242 | 6 | 0.8279 | + | w0008 |
| NC_008532.1 | 1012272 | 1012434 | 160 | 2.1193 | 3 | 1.0596 | + | w0006 |
| NC_008532.1 | 1241033 | 1241128 | 94 | 1.6609 | 3 | 0.8305 | + | w0007 |
| NC_000964.3 | 626323 | 626376 | 53 | 1.6595 | 2 | 1.047 | - | w0018 |
| MW218148.1 | 207740 | 207883 | 143 | 1.6311 | 2 | 1.0291 | - | w0012 |
| NC_028808.2 | 29211 | 29264 | 52 | 1.538 | 2 | 0.9704 | - | w0003 |

25 per-window clusters merged into 13 candidates (overlapping windows never counted as independent peaks).

## 2. Peak support

| candidate | peaks | peak lengths (nt) | gaps between peaks (nt) |
|---|---:|---|---|
| NC_008532.1:649126-650211 | 17 | 31,29,30,29,30,30,30,31,30,31,30,31,31,30,31,30,26 | 35,36,37,36,36,36,37,35,36,35,37,35,36,35,36,36 |
| NC_008532.1:1376866-1378094 | 12 | 13,14,35,35,36,33,34,34,35,34,34,23 | 9,327,31,31,30,36,31,33,31,31,276 |
| MW248466.1:213994-214785 | 5 | 9,10,28,36,21 | 204,168,173,139 |
| NC_028808.2:27682-28002 | 5 | 21,37,37,45,41 | 15,57,10,56 |
| NC_000964.3:624359-625143 | 4 | 15,20,19,33 | 395,172,130 |
| NC_000964.3:2817867-2818034 | 5 | 17,18,16,22,16 | 13,3,19,35 |
| MW218148.1:9951-10370 | 3 | 45,45,39 | 156,134 |
| NC_008532.1:745310-746042 | 6 | 21,12,12,19,8,17 | 59,31,399,109,45 |
| NC_008532.1:1012272-1012434 | 3 | 21,24,26 | 83,6 |
| NC_008532.1:1241033-1241128 | 3 | 8,8,21 | 52,5 |
| NC_000964.3:626323-626376 | 2 | 25,15 | 13 |
| MW218148.1:207740-207883 | 2 | 17,21 | 105 |
| NC_028808.2:29211-29264 | 2 | 9,25 | 19 |

145 raw peaks over all windows and orientations; median per window-orientation 8.

## 3. Label localization (labels used only to read the result)

| record | family | label | best overlapping candidate | overlap nt | fraction of label | fraction of candidate |
|---|---|---|---|---:|---:|---:|
| w0001 | CRISPR | 649124-650217 | 649126-650211 | 1085 | 0.99 | 1.00 |
| w0002 | CRISPR | 1377228-1377794 | 1376866-1378094 | 566 | 1.00 | 0.46 |
| w0003 | VIPR | 29135-29229 | 29211-29264 | 18 | 0.19 | 0.34 |
| w0004 | ART | 9447-10646 | 9951-10370 | 419 | 0.35 | 1.00 |
| w0005 | ART | 213934-215088 | 213994-214785 | 791 | 0.69 | 1.00 |

5/5 labels overlapped by a candidate. Overlap is loose recovery, not boundary reconstruction.

Candidates in background windows (no label; not false positives by construction - backgrounds are unlabeled real DNA):

- NC_000964.3:624359-625143 score 2.4649 in w0018 (Bsub168)
- NC_000964.3:2817867-2818034 score 2.3634 in w0019 (Bsub168)
- NC_008532.1:745310-746042 score 2.3242 in w0008 (LMD9)
- NC_008532.1:1012272-1012434 score 2.1193 in w0006 (LMD9)
- NC_008532.1:1241033-1241128 score 1.6609 in w0007 (LMD9)
- NC_000964.3:626323-626376 score 1.6595 in w0018 (Bsub168)
- MW218148.1:207740-207883 score 1.6311 in w0012 (SA1)

## 4. Window ranking (peak-cluster vs unchanged block baselines)

| rank | record | family | peak-cluster | blk256_mean_max | blk1024_mean_max | blk256_top20_max | blk1024_top20_max |
|---:|---|---|---:|---:|---:|---:|---:|
| 1 | w0001 | CRISPR | 5.465 | 0.82161 | 0.78449 | 1.38611 | 1.38610 |
| 2 | w0002 | CRISPR | 4.663 | 0.84004 | 0.52811 | 1.38610 | 1.38597 |
| 3 | w0005 | ART | 3.362 | 0.61972 | 0.47114 | 1.36778 | 1.30960 |
| 4 | w0003 | VIPR | 2.930 | 0.88373 | 0.42510 | 1.35821 | 1.31090 |
| 5 | w0018 | background | 2.465 | 0.44037 | 0.28353 | 1.26181 | 1.11785 |
| 6 | w0019 | background | 2.363 | 0.62521 | 0.33490 | 1.36852 | 1.23665 |
| 7 | w0004 | ART | 2.332 | 0.56133 | 0.31671 | 1.36736 | 1.23519 |
| 8 | w0008 | background | 2.324 | 0.40481 | 0.24952 | 1.25352 | 1.07584 |
| 9 | w0006 | background | 2.119 | 0.47310 | 0.24159 | 1.30842 | 0.99732 |
| 10 | w0007 | background | 1.661 | 0.35764 | 0.23944 | 1.21591 | 0.95673 |
| 11 | w0012 | background | 1.631 | 0.52366 | 0.39776 | 1.20913 | 1.09931 |
| 12 | w0009 | background | 0.000 | 0.39365 | 0.29407 | 1.23955 | 0.96010 |
| 13 | w0010 | background | 0.000 | 0.27571 | 0.20538 | 1.02263 | 0.81869 |
| 14 | w0011 | background | 0.000 | 0.92398 | 0.67936 | 1.38232 | 1.37836 |
| 15 | w0013 | background | 0.000 | 0.41704 | 0.33512 | 1.09500 | 1.01079 |
| 16 | w0014 | background | 0.000 | 0.46616 | 0.35092 | 1.08914 | 0.98878 |
| 17 | w0015 | background | 0.000 | 0.58168 | 0.36636 | 1.25719 | 1.07822 |
| 18 | w0016 | background | 0.000 | 0.31014 | 0.27929 | 0.99287 | 0.94306 |
| 19 | w0017 | background | 0.000 | 0.29767 | 0.21555 | 1.06223 | 0.89731 |

Non-zero peak-cluster score: 5/5 positives, 6/14 backgrounds. P(positive > background) = 0.97 (ties at 0 counted as 0.5; descriptive only at n=5).

## 5. Random partition vs exact, evaluated on this scorer

| RP run | candidates | exact candidates recovered | RP candidates not in exact | window-score Spearman | positives with non-zero score |
|---|---:|---|---:|---:|---:|
| rp_seed20261014 | 13 | 13/13 | 0 | 0.984 | 5/5 |
| rp_seed1 | 13 | 13/13 | 0 | 0.996 | 5/5 |
| rp_seed2 | 13 | 13/13 | 0 | 0.996 | 5/5 |

Exact-supported candidates missed by RP: none

## 6. Limitations found in this run

**The w0011 zero is an artefact, not a discrimination result.** The window-local
threshold is `median + 3 x robust SD`. For w0011 (the 12-tRNA cluster) the smoothed
track is broadly elevated and highly dispersed (median 0.358, robust SD 0.413), so
tau = 1.596 — which is above the theoretical ceiling of this measurement, ln(4) =
1.386. No position can ever exceed it, so the window yields zero peaks by
construction. w0011 ranked #1-#3 on every block baseline and scores 0 here, but
that must not be read as the scorer correctly rejecting structured-RNA clusters.
2/38 tracks are affected, both of them w0011.

Fixing this (clipping tau below ln(4), or using an upper-quantile rather than a
MAD-based threshold) is a parameter change and was deliberately not made this
round. It is the first thing to settle before this scorer is used for anything
beyond development.

**Other limits.** The window-local baseline makes the score relative within a
window: a window that is uniformly predictable scores low regardless of its
content. Candidate boundaries are peak extents, not repeat boundaries (w0002's
candidate spans 1226 nt for a 566-nt label, fraction-of-candidate 0.46). seed19 is
development data; no threshold was calibrated and no FDR is claimed.
