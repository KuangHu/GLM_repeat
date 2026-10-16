# attention_mask_policy=auto regression: acceptance criteria
Declared 2026-10-14, BEFORE the auto-policy regression run was scored.

Reference: results/benchmark/timing6_exact (policy=always, batch 1, bf16, A40,
model 08754cba59a1f97d517f873fad6c672d2b1abdc7). Immutable; never overwritten.
Candidate: same windows, both orientations, policy=auto, batch 1, same model/dtype.

Observed bf16 noise floor (same code, batch 1 vs batch 4, w0001 +):
per-base max |diff| 0.064, mean |diff| 0.005 nats.

PASS requires, over all 6 windows x 2 orientations (36,000 positions):
1. per-base |diff| p99 <= 0.03 and max <= 0.15 nats (about 2x noise floor)
2. every block feature (scales 64-2048, mean_max and top20_max), window mean,
   q90 and max16: |diff| <= 0.01 nats
3. label-region mean |diff| <= 0.005 nats for each positive and orientation
4. label-blind ranking of the 6 windows unchanged for every summary feature
   (ties within 0.01 excepted)
5. all tracks finite; identical n_variants and mask coverage

If PASS: DEFAULT_ATTENTION_MASK_POLICY -> 'auto'. If FAIL: keep 'always', report.

## Outcome (2026-10-14)
Verdict: FAIL on criterion 1 only (per-base p99 0.0303 > 0.03; max 0.091 ok).
Criteria 2-5 passed (features <= 0.006, label means <= 0.0034, no ranking flips).
Context: on w0001 + the auto-vs-reference p99 (0.0252) is below the batch-1-vs-4
noise p99 (0.0268); the p99 limit had ~10% headroom over a single-track noise
estimate. Per the pre-declared rule the default stays 'always'. Changing the
criterion is left to the project owner. Details: timing6_auto_vs_reference/.
