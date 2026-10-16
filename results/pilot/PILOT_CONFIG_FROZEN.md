# Frozen scorer and RP configuration for the 10-Mb pilot

Frozen 2026-10-15, before any pilot window was scored. Nothing in this file is
re-tuned during the pilot.

## Model and inference

| item | value |
|---|---|
| model | tattabio/gLM2_650M @ 08754cba59a1f97d517f873fad6c672d2b1abdc7 |
| dtype | bfloat16 (A40) |
| attention_mask_policy | **always** (the validated reference path; `auto` remains available but unused) |
| scoring mode | random_partition, groups 20, repeats 2, seed 20261014 |
| microbatch | 1 |
| orientations | both, scored independently, mapped to forward coordinates |

## Windowing (WORKFLOW_SETUP.md section 5 / NEXT_STEPS section 5)

window 3000, stride 1500, min valid segment 40, split at every non-ACGT base,
0-based half-open, no CDS or low-complexity masking, no redundant contained tail.
The same `scripts/make_windows.py` pipeline is used for calibration and held-out
material.

## Peak-cluster scorer (unchanged from results/scorer/PEAK_CLUSTER_DESIGN.md)

SMOOTH_NT 16, THRESH_K 3.0 (window_local only), MIN_PEAK_NT 8, MAX_GAP_NT 500,
MIN_PEAKS 2; cluster score = median(peak_mean) * log2(1 + n_peaks);
window score = max over clusters and orientations.
Threshold policy for the pilot: `background_calibrated`.

## Candidate merging (identical for calibration and held-out)

Peaks and clusters are emitted in genomic coordinates. Clusters overlapping on the
same contig are merged into one candidate regardless of which window or which
strand produced them. Overlapping windows and the two strands therefore never
contribute independent candidates.

## Density definition

Candidate density is per **unique scanned genomic Mb**: the union of genomic
positions covered by emitted windows on a contig, counted once. It is NOT the sum
of window lengths (which double-counts the 50% stride overlap) and NOT per
orientation.

## Eligibility

A candidate is eligible on the peak-cluster score alone. **Copy dependence is a
diagnostic and is NOT a candidate eligibility criterion.** No annotation, motif,
repeat finder, CDS call or low-complexity filter gates any candidate. Annotation is
applied only after scoring, for description.

## Archived references, unchanged

- results/benchmark/timing6_exact (immutable exact reference)
- results/scorer/exact and results/scorer/report_exact (original window_local run)
- results/calibration/thresholds/FROZEN.md (round-6 thresholds)
- provenance/reference_score.bundle_v0.1.py, provenance/copy_dependence.uploaded_truncated_186lines.py

## Data separation

Calibration and held-out material share no genome and no genus. Separation is by
genome / related-sequence cluster, never by random window assignment.

| set | genomes | approx Mb |
|---|---|---|
| calibration | *Methanocaldococcus jannaschii* NC_000909.1; *Deinococcus radiodurans* NC_001263.1; *Staphylococcus* phage K NC_005880.2; phage Twort NC_007021.1 | 4.59 |
| held-out | *Pyrococcus furiosus* NC_003413.1; *Thermus thermophilus* HB8 NC_006461.1; Madawaska MW349129.1; LY01 OR836606.1; S6 LC680885.1; PALS_2 MN091626.1; UFV_DC4 MZ779063.1 | 5.07 |

Type matching: each held-out genome class (hyperthermophilic archaeon,
Deinococcus-Thermus bacterium, Staphylococcus phage) has a counterpart class in
calibration, drawn from a different genus. Calibration genomes are ordinary genomic
DNA and may themselves contain real repeats, including CRISPR arrays; they are
calibration material, not verified negatives. Reported quantities are candidate
yield per Mb, never FDR.
