# Seed panel (19 windows): exact masked scoring

Development/smoke panel only: not calibration or test data. Two CRISPR loci share
the LMD9 genome; both ART loci share one Staphylococcus phage lineage. Backgrounds
are unlabeled real DNA (not confirmed negatives). Scores are masked predictability
(nats above uniform), not repeat probabilities. No thresholds were selected.

## Provenance

- model 08754cba59a1f97d517f873fad6c672d2b1abdc7, dtype torch.bfloat16, attention_mask_policy always, mode exact, batch 1, code not recorded: bundle v0.1 reference_score.py sha256 012c0e1a...
- model 08754cba59a1f97d517f873fad6c672d2b1abdc7, dtype torch.bfloat16, attention_mask_policy always, mode exact, batch 1, code {"reference_score.py": "2357e1dcffbd4b053ee0a2a78464991f87c8c3a43d5c8445af4a5e13ccb62eb0", "timing_score.py": "27f3fe918922267b8318bf33484a87d674cdf920e7c3581b612c858d00b5c030"}
- tracks: 38 (19 windows x orientations); missing: none

## Timing

- per 3,000-nt orientation: median 809 s (range 804-836); peak GPU 2204 MB
- total GPU inference 8.58 h for 38 orientations

## Label-level signal (uses labels)

| record | family | locus | orient | label nt | label mean | outside mean | difference |
|---|---|---|---|---:|---:|---:|---:|
| w0001 | CRISPR | LMD9_CRISPR_649125 | + | 1093 | 0.789 | 0.167 | +0.623 |
| w0001 | CRISPR | LMD9_CRISPR_649125 | - | 1093 | 0.786 | 0.156 | +0.630 |
| w0002 | CRISPR | LMD9_CRISPR_1377229 | + | 566 | 0.822 | 0.166 | +0.656 |
| w0002 | CRISPR | LMD9_CRISPR_1377229 | - | 566 | 0.812 | 0.162 | +0.650 |
| w0003 | VIPR | SUSP1_vrRNA_35TR_2 | + | 94 | 0.455 | 0.270 | +0.185 |
| w0003 | VIPR | SUSP1_vrRNA_35TR_2 | - | 94 | 0.419 | 0.277 | +0.141 |
| w0004 | ART | SA1_ART | + | 1199 | 0.327 | 0.269 | +0.057 |
| w0004 | ART | SA1_ART | - | 1199 | 0.324 | 0.274 | +0.050 |
| w0005 | ART | MarsHill_ART | + | 1154 | 0.433 | 0.336 | +0.097 |
| w0005 | ART | MarsHill_ART | - | 1154 | 0.414 | 0.339 | +0.075 |

Annotated sub-regions (orientation-averaged track):

| record | feature | genomic [start0,end0) | sub-region mean | rest of label mean |
|---|---|---|---:|---:|
| w0003 | 12_YNNGG_units_of_mapped_94nt_vrRNA | 29169-29229 | 0.650 | 0.061 |
| w0004 | ART_core | 9597-..-10371 (5 parts) | 1.006 | 0.283 |
| w0005 | ART_core | 213994-..-214801 (5 parts) | 1.226 | 0.372 |

Localisation: position of the best 256-nt block (orientation-averaged track, chosen label-blind) relative to the label:

- w0001 (CRISPR): best block [1600,1856) mean 0.819; overlaps label by 256 nt
- w0002 (CRISPR): best block [1216,1472) mean 0.836; overlaps label by 255 nt
- w0003 (VIPR): best block [0,256) mean 0.883; overlaps label by 0 nt
- w0004 (ART): best block [1600,1856) mean 0.549; overlaps label by 256 nt
- w0005 (ART): best block [1536,1792) mean 0.617; overlaps label by 256 nt

## Label-blind window ranking (labels used only to read the result)

Features fixed in advance: mean, q90, blk256_mean_max, blk1024_mean_max, blk256_top20_max, blk1024_top20_max. 5 positives, 14 backgrounds.

| record | family | source | mean | q90 | blk256_mean_max | blk1024_mean_max | blk256_top20_max | blk1024_top20_max |
|---|---|---|---:|---:|---:|---:|---:|---:|
| w0001 | CRISPR | LMD9 | 0.393 (#2) | 1.386 (#1) | 0.822 (#4) | 0.784 (#1) | 1.386 (#1) | 1.386 (#1) |
| w0011 | background | SUSP1 | 0.476 (#1) | 1.370 (#3) | 0.924 (#1) | 0.679 (#2) | 1.382 (#3) | 1.378 (#3) |
| w0002 | CRISPR | LMD9 | 0.289 (#9) | 1.384 (#2) | 0.840 (#3) | 0.528 (#3) | 1.386 (#2) | 1.386 (#2) |
| w0005 | ART | MarsHill | 0.373 (#3) | 1.032 (#5) | 0.620 (#6) | 0.471 (#4) | 1.368 (#5) | 1.310 (#5) |
| w0003 | VIPR | SUSP1 | 0.282 (#10) | 1.121 (#4) | 0.884 (#2) | 0.425 (#5) | 1.358 (#7) | 1.311 (#4) |
| w0012 | background | SA1 | 0.337 (#4) | 0.996 (#6) | 0.524 (#9) | 0.398 (#6) | 1.209 (#14) | 1.099 (#9) |
| w0015 | background | MarsHill | 0.335 (#5) | 0.957 (#9) | 0.582 (#7) | 0.366 (#7) | 1.257 (#10) | 1.078 (#10) |
| w0014 | background | MarsHill | 0.304 (#6) | 0.881 (#11) | 0.466 (#11) | 0.351 (#8) | 1.089 (#16) | 0.989 (#14) |
| w0013 | background | SA1 | 0.294 (#8) | 0.879 (#12) | 0.417 (#13) | 0.335 (#9) | 1.095 (#15) | 1.011 (#12) |
| w0019 | background | Bsub168 | 0.234 (#12) | 0.918 (#10) | 0.625 (#5) | 0.335 (#10) | 1.369 (#4) | 1.237 (#6) |
| w0004 | ART | SA1 | 0.294 (#7) | 0.958 (#8) | 0.561 (#8) | 0.317 (#11) | 1.367 (#6) | 1.235 (#7) |
| w0009 | background | LMD9 | 0.250 (#11) | 0.852 (#14) | 0.394 (#15) | 0.294 (#12) | 1.240 (#12) | 0.960 (#15) |
| w0018 | background | Bsub168 | 0.232 (#13) | 0.974 (#7) | 0.440 (#12) | 0.284 (#13) | 1.262 (#9) | 1.118 (#8) |
| w0016 | background | MG1655 | 0.220 (#14) | 0.840 (#16) | 0.310 (#17) | 0.279 (#14) | 0.993 (#19) | 0.943 (#17) |
| w0008 | background | LMD9 | 0.190 (#17) | 0.864 (#13) | 0.405 (#14) | 0.250 (#15) | 1.254 (#11) | 1.076 (#11) |
| w0006 | background | LMD9 | 0.211 (#15) | 0.842 (#15) | 0.473 (#10) | 0.242 (#16) | 1.308 (#8) | 0.997 (#13) |
| w0007 | background | LMD9 | 0.194 (#16) | 0.814 (#17) | 0.358 (#16) | 0.239 (#17) | 1.216 (#13) | 0.957 (#16) |
| w0017 | background | MG1655 | 0.165 (#19) | 0.774 (#18) | 0.298 (#18) | 0.216 (#18) | 1.062 (#17) | 0.897 (#18) |
| w0010 | background | SUSP1 | 0.190 (#18) | 0.761 (#19) | 0.276 (#19) | 0.205 (#19) | 1.023 (#18) | 0.819 (#19) |

| feature | positive ranks (of 19) | P(pos > bg) | best-ranked background |
|---|---|---:|---|
| mean | CRISPR:2, ART:3, ART:7, CRISPR:9, VIPR:10 | 0.77 | w0011 (SUSP1, #1) |
| q90 | CRISPR:1, CRISPR:2, VIPR:4, ART:5, ART:8 | 0.93 | w0011 (SUSP1, #3) |
| blk256_mean_max | VIPR:2, CRISPR:3, CRISPR:4, ART:6, ART:8 | 0.89 | w0011 (SUSP1, #1) |
| blk1024_mean_max | CRISPR:1, CRISPR:3, ART:4, VIPR:5, ART:11 | 0.87 | w0011 (SUSP1, #2) |
| blk256_top20_max | CRISPR:1, CRISPR:2, ART:5, ART:6, VIPR:7 | 0.91 | w0011 (SUSP1, #3) |
| blk1024_top20_max | CRISPR:1, CRISPR:2, VIPR:4, ART:5, ART:7 | 0.94 | w0011 (SUSP1, #3) |

## Positive vs its paired backgrounds (same genome, GC within 2 points)

| positive | family | mean | q90 | blk256_mean_max | blk1024_mean_max | blk256_top20_max | blk1024_top20_max |
|---|---|---:|---:|---:|---:|---:|---:|
| w0001 vs w0006,w0007 | CRISPR | 0.393 vs 0.211/0.194 | 1.386 vs 0.842/0.814 | 0.822 vs 0.473/0.358 | 0.784 vs 0.242/0.239 | 1.386 vs 1.308/1.216 | 1.386 vs 0.997/0.957 |
| w0002 vs w0008,w0009 | CRISPR | 0.289 vs 0.190/0.250 | 1.384 vs 0.864/0.852 | 0.840 vs 0.405/0.394 | 0.528 vs 0.250/0.294 | 1.386 vs 1.254/1.240 | 1.386 vs 1.076/0.960 |
| w0003 vs w0010,w0011 | VIPR | 0.282 vs 0.190/0.476 | 1.121 vs 0.761/1.370 | 0.884 vs 0.276/0.924 | 0.425 vs 0.205/0.679 | 1.358 vs 1.023/1.382 | 1.311 vs 0.819/1.378 |
| w0004 vs w0012,w0013 | ART | 0.294 vs 0.337/0.294 | 0.958 vs 0.996/0.879 | 0.561 vs 0.524/0.417 | 0.317 vs 0.398/0.335 | 1.367 vs 1.209/1.095 | 1.235 vs 1.099/1.011 |
| w0005 vs w0014,w0015 | ART | 0.373 vs 0.304/0.335 | 1.032 vs 0.881/0.957 | 0.620 vs 0.466/0.582 | 0.471 vs 0.351/0.366 | 1.368 vs 1.089/1.257 | 1.310 vs 0.989/1.078 |

## Separate audit (NOT ground truth, not used in any metric above)

- SUSP1_offlabel_c1: NC_028808.2:[27682,27906) sequence_supported_repeat_candidate; vrRNA identity unconfirmed. exact 39-nt direct repeat x2 (27714,27857); 16 YNNGG matches on minus strand in 240 nt vs 3.7 genome-wide expectation; exact-track 64-nt top-6 blocks all here; orientation-avg mean 0.922 vs 0.279 window

n is tiny (5 positives, 3 families, 2 independent lineages for CRISPR/ART): all numbers are descriptive.
- Top-ranked background w0011 is an annotated 12-tRNA cluster: see
  results/audit/background_high_score_audit.md. It stays labelled background.

Attention policy: these tracks use policy=always (validated reference path). The
auto-policy tracks for the same 13 windows (exact_auto/) failed the pre-declared
per-base p99 criterion (0.0316 > 0.03) while passing feature, label and ranking
criteria; see results/regression/.
