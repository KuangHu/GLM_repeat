# Exact-score diagnostics (development only, not calibrated)

model tattabio/gLM2_650M@08754cba59  mode exact  batch 1  dtype torch.bfloat16  device NVIDIA A40
model load 9.1s

## Timing per orientation

| record | family | orient | variants | inference s | variants/s | peak GPU MB |
|---|---|---|---:|---:|---:|---:|
| w0001 | CRISPR | + | 3000 | 806.9 | 3.7 | 2204 |
| w0001 | CRISPR | - | 3000 | 807.5 | 3.7 | 2204 |
| w0006 | background | + | 3000 | 807.2 | 3.7 | 2204 |
| w0006 | background | - | 3000 | 808.1 | 3.7 | 2204 |
| w0003 | VIPR | + | 3000 | 809.4 | 3.7 | 2204 |
| w0003 | VIPR | - | 3000 | 809.4 | 3.7 | 2204 |
| w0010 | background | + | 3000 | 809.3 | 3.7 | 2204 |
| w0010 | background | - | 3000 | 809.1 | 3.7 | 2204 |
| w0004 | ART | + | 3000 | 809.2 | 3.7 | 2204 |
| w0004 | ART | - | 3000 | 807.7 | 3.7 | 2204 |
| w0012 | background | + | 3000 | 809.0 | 3.7 | 2204 |
| w0012 | background | - | 3000 | 808.9 | 3.7 | 2204 |

Total inference 9702s for 36000 scored nt-orientations (3.7 nt/s).

## Label region vs rest of window (above_uniform, nats)

| record | family | orient | label mean | outside mean | label max16 |
|---|---|---|---:|---:|---:|
| w0001 | CRISPR | + | 0.789 | 0.167 | 1.386 |
| w0001 | CRISPR | - | 0.786 | 0.156 | 1.386 |
| w0003 | VIPR | + | 0.455 | 0.270 | 0.838 |
| w0003 | VIPR | - | 0.419 | 0.277 | 0.894 |
| w0004 | ART | + | 0.327 | 0.269 | 1.358 |
| w0004 | ART | - | 0.324 | 0.274 | 1.364 |

## Positive vs paired background (window score = max over orientations)

| positive | background | mean | max16 | blk64_mean_max | blk256_mean_max | blk1024_mean_max | blk64_top20_max | blk256_top20_max | blk1024_top20_max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w0001 (CRISPR) | w0006 | 0.393 vs 0.211 | 1.386 vs 1.320 | 0.902 vs 0.989 | 0.822 vs 0.473 | 0.784 vs 0.242 | 1.386 vs 1.367 | 1.386 vs 1.308 | 1.386 vs 0.997 |
| w0003 (VIPR) | w0010 | 0.282 vs 0.190 | 1.337 vs 1.245 | 1.131 vs 0.487 | 0.884 vs 0.276 | 0.425 vs 0.205 | 1.364 vs 1.334 | 1.358 vs 1.023 | 1.311 vs 0.819 |
| w0004 (ART) | w0012 | 0.294 vs 0.337 | 1.364 vs 1.163 | 1.095 vs 0.747 | 0.561 vs 0.524 | 0.317 vs 0.398 | 1.373 vs 1.307 | 1.367 vs 1.209 | 1.235 vs 1.099 |

n=3 pairs: descriptive only. Backgrounds are unlabeled real DNA, not verified negatives.

## Notes (manual inspection, 2026-10-14)

- Exact track is reproducible bit-for-bit at batch 1 across separate jobs/GPUs.
- Timing: 3000 nt exact = 3000 variants, ~808 s per orientation on A40, bf16, batch 1,
  peak 2.2 GB. Batch 4 gives no speedup (GPU saturated at 3k tokens).
- Omitting the all-ones attention mask (timing_score.py --drop-full-attention-mask)
  runs 2.6x faster (308 s bs1, 284 s bs8); differences vs reference (max 0.05,
  mean 0.005 nats) are the same size as the batch-1 vs batch-4 bf16 noise.
- CRISPR w0001: all top-6 64-nt blocks lie inside the annotated array.
- ART w0004: all top-6 64-nt blocks inside the label; the five 14-nt cores average
  1.006 vs 0.283 for the rest of the label region (signal concentrates on cores, so
  the whole-label mean is only slightly above background).
- VIPR w0003: the 60-nt YNNGG subregion averages 0.65 vs 0.06 for the rest of the
  vrRNA. The window's strongest blocks are at [0,224), outside the label; that
  region contains a ~38-nt element present twice and YNNGG-like units on the minus
  strand: plausibly unlabelled vrRNA/repeat units (label covers one RNA only).
  Not confirmed.
- Short-scale statistics (max16, blk64) do not separate positives from paired
  background; 256-1024-nt block statistics did for CRISPR and VIPR, not for ART
  (blk1024_mean 0.317 vs 0.398). n=1 per family: no generalisation.
