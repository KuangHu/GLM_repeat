# Exact audit of sampled held-out windows

24 windows x 2 orientations = 48 tracks, scored exactly in their original untrimmed 3,000-nt contexts.

**The RP operating point was not transferred to exact.** Peaks below are called on both modes at a matched tau = 0.8 solely so that peak structure is comparable. That is not a calibrated exact operating point: exact has no calibration at this density, and establishing one would need ~1 Mb of exact background (~150 GPU-hours), which was explicitly declined.

**Exact is not treated as ground truth.** It is a reference measurement of the same model under a different masking scheme, not an independent observation of the genome.

## Agreement by stratum (median over tracks)

| stratum | tracks | bias RP-exact | MAE | p99 abs | Pearson | Pearson 16-nt | peak-span Jaccard | cluster-span Jaccard |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| high | 16 | -0.0016 | 0.0377 | 0.2560 | 0.9941 | 0.9989 | 0.993 | 1.000 |
| adjacent | 16 | -0.0091 | 0.0705 | 0.4282 | 0.9836 | 0.9841 | 0.590 | 0.987 |
| non_hit | 16 | -0.0054 | 0.0643 | 0.3694 | 0.9818 | 0.9829 | 1.000 | 1.000 |

## Cluster-level correspondence at matched tau

| stratum | tracks | exact has cluster | RP has cluster | both | neither |
|---|---:|---:|---:|---:|---:|
| high | 16 | 16 | 16 | 16 | 0 |
| adjacent | 16 | 6 | 8 | 5 | 7 |
| non_hit | 16 | 0 | 0 | 0 | 16 |

Tracks where only one mode forms a cluster: 4 of 48.
- NC_006461.1:NC_006461.1:411063-414063 (+, adjacent): exact 1 cluster(s) best 1.6892, RP 0 best 0.0
- MZ779063.1:MZ779063.1:172500-175500 (-, adjacent): exact 0 cluster(s) best 0.0, RP 1 best 1.3347
- NC_006461.1:NC_006461.1:1404468-1407468 (-, adjacent): exact 0 cluster(s) best 0.0, RP 1 best 1.3147
- NC_006461.1:NC_006461.1:399063-402063 (-, adjacent): exact 0 cluster(s) best 0.0, RP 1 best 1.3046

## Reading

The `adjacent` stratum is the one that matters: it was drawn at the cluster-formation boundary (lowest surviving RP cluster scores, 1.305-1.337), so small track differences can flip a detection either way. Disagreement there is expected and is the quantity of interest, not a defect in either mode.

Sampling-rule note: the `adjacent` stratum is NOT "closest in value to tau"; see results/pilot/heldout/exact_audit/sampling.json.

---

## CORRECTIONS 2026-10-15

### 1. Exact-only clusters were not absent from the audit

The closeout said "no window called by exact that RP missed in this sample". That
is wrong for the audit as a whole. Exact-only clusters by stratum:

| stratum | exact-only | RP-only |
|---|---:|---:|
| high | 0 | 0 |
| adjacent | **1** | 3 |
| non_hit | 0 | 0 |

The correct statement: **in the random non-hit stratum** (and in the high stratum)
no exact-only cluster appeared. In the threshold-adjacent stratum one track
produced a cluster under exact and none under RP
(NC_006461.1:411063-414063 (+), exact best 1.6892). So exact did call something RP
missed; it happened only at the detection boundary, not among random non-hits.

### 2. The non-hit Jaccard of 1.000 is mostly empty-set agreement

`span_overlap` returns 1.0 when both interval lists are empty. In the non_hit
stratum **12 of 16 tracks had zero peaks in both modes**, so their Jaccard is
vacuous, not evidence of overlapping real peaks.

| stratum | tracks | both-zero-peak tracks | median peaks (exact / RP) |
|---|---:|---:|---|
| high | 16 | 0 | 29 / 29 |
| adjacent | 16 | 2 | 2 / 2 |
| non_hit | 16 | **12** | 0 / 0 |

Only the `high` stratum's Jaccard of 0.993 reflects agreement between substantial
real peak sets. The adjacent figure (0.590) is over small peak sets, two of which
are empty on both sides.
