# Track comparison

reference: results/benchmark/timing6_exact  policy=always (pre-policy code) mode=exact batch=1 dtype=torch.bfloat16
candidate: results/regression/timing6_auto  policy=auto mode=exact batch=1 dtype=torch.bfloat16
model revision: ref 08754cba59a1f97d517f873fad6c672d2b1abdc7 / cand 08754cba59a1f97d517f873fad6c672d2b1abdc7
candidate code: {"timing_score.py": "27f3fe918922267b8318bf33484a87d674cdf920e7c3581b612c858d00b5c030", "reference_score.py": "2357e1dcffbd4b053ee0a2a78464991f87c8c3a43d5c8445af4a5e13ccb62eb0"}
tracks compared: 12; missing: none

## Verdict: **FAIL**

- 1_perbase_p99<=0.03: FAILED
- 1_perbase_max<=0.15: ok
- 2_features<=0.01: ok
- 3_label_mean<=0.005: ok
- 4_ranking_unchanged: ok
- 5_finite_and_coverage: ok

## Per-base |diff| (nats), all positions

| n | mean | p50 | p90 | p99 | p99.9 | max |
|---:|---:|---:|---:|---:|---:|---:|
| 36000 | 0.00738 | 0.00515 | 0.01766 | 0.03033 | 0.04873 | 0.09118 |

## Per track

| record | orient | mean | p99 | max | argmax | corr |
|---|---|---:|---:|---:|---:|---:|
| w0001 | + | 0.00563 | 0.02521 | 0.05349 | 2869 | 0.99991 |
| w0001 | - | 0.00562 | 0.02561 | 0.04582 | 2410 | 0.99992 |
| w0003 | + | 0.00721 | 0.03277 | 0.05569 | 215 | 0.99986 |
| w0003 | - | 0.00784 | 0.03415 | 0.06578 | 168 | 0.99983 |
| w0004 | + | 0.00720 | 0.02818 | 0.05969 | 1665 | 0.99987 |
| w0004 | - | 0.00747 | 0.03163 | 0.09118 | 1665 | 0.99986 |
| w0006 | + | 0.00783 | 0.02843 | 0.07194 | 2075 | 0.99983 |
| w0006 | - | 0.00759 | 0.02931 | 0.05113 | 1065 | 0.99981 |
| w0010 | + | 0.00765 | 0.02823 | 0.08954 | 1959 | 0.99977 |
| w0010 | - | 0.00762 | 0.02863 | 0.07893 | 459 | 0.99978 |
| w0012 | + | 0.00833 | 0.03408 | 0.06369 | 935 | 0.99986 |
| w0012 | - | 0.00854 | 0.03676 | 0.06736 | 912 | 0.99983 |

## Max |diff| of window/block features (over tracks)

| feature | max abs diff |
|---|---:|
| mean | 0.00035 |
| sd | 0.00054 |
| q50 | 0.00330 |
| q90 | 0.00340 |
| max16 | 0.00582 |
| blk64_mean_max | 0.00242 |
| blk64_top20_max | 0.00300 |
| blk128_mean_max | 0.00175 |
| blk128_top20_max | 0.00155 |
| blk256_mean_max | 0.00139 |
| blk256_top20_max | 0.00090 |
| blk512_mean_max | 0.00083 |
| blk512_top20_max | 0.00107 |
| blk1024_mean_max | 0.00067 |
| blk1024_top20_max | 0.00081 |
| blk2048_mean_max | 0.00046 |
| blk2048_top20_max | 0.00041 |

## Label-region mean

| record | family | orient | ref | cand | abs diff |
|---|---|---|---:|---:|---:|
| w0001 | CRISPR | + | 0.7892 | 0.7889 | 0.00025 |
| w0001 | CRISPR | - | 0.7857 | 0.7858 | 0.00007 |
| w0003 | VIPR | + | 0.4552 | 0.4548 | 0.00034 |
| w0003 | VIPR | - | 0.4189 | 0.4224 | 0.00343 |
| w0004 | ART | + | 0.3268 | 0.3271 | 0.00038 |
| w0004 | ART | - | 0.3240 | 0.3240 | 0.00001 |

## Label-blind ranking (window score = max over orientations)

Order flips beyond 0.01 tie tolerance: none
