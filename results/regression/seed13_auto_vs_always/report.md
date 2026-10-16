# Track comparison

reference: results/benchmark/seed19/exact_always  policy=always mode=exact batch=1 dtype=torch.bfloat16
candidate: results/benchmark/seed19/exact_auto  policy=auto mode=exact batch=1 dtype=torch.bfloat16
model revision: ref 08754cba59a1f97d517f873fad6c672d2b1abdc7 / cand 08754cba59a1f97d517f873fad6c672d2b1abdc7
candidate code: {"timing_score.py": "27f3fe918922267b8318bf33484a87d674cdf920e7c3581b612c858d00b5c030", "reference_score.py": "2357e1dcffbd4b053ee0a2a78464991f87c8c3a43d5c8445af4a5e13ccb62eb0"}
tracks compared: 26; missing: none

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
| 78000 | 0.00764 | 0.00534 | 0.01810 | 0.03157 | 0.05136 | 0.12116 |

## Per track

| record | orient | mean | p99 | max | argmax | corr |
|---|---|---:|---:|---:|---:|---:|
| w0002 | + | 0.00633 | 0.02663 | 0.05101 | 2781 | 0.99990 |
| w0002 | - | 0.00636 | 0.02668 | 0.05148 | 1847 | 0.99989 |
| w0005 | + | 0.00727 | 0.03000 | 0.06286 | 1705 | 0.99987 |
| w0005 | - | 0.00753 | 0.03193 | 0.07406 | 1432 | 0.99986 |
| w0007 | + | 0.00789 | 0.02857 | 0.04920 | 1048 | 0.99982 |
| w0007 | - | 0.00795 | 0.02914 | 0.05157 | 2577 | 0.99981 |
| w0008 | + | 0.00765 | 0.02983 | 0.06790 | 2206 | 0.99983 |
| w0008 | - | 0.00748 | 0.02914 | 0.05880 | 2544 | 0.99983 |
| w0009 | + | 0.00779 | 0.03285 | 0.05188 | 2726 | 0.99980 |
| w0009 | - | 0.00762 | 0.02758 | 0.05621 | 508 | 0.99979 |
| w0011 | + | 0.00751 | 0.04659 | 0.12116 | 600 | 0.99986 |
| w0011 | - | 0.00762 | 0.04436 | 0.09690 | 867 | 0.99986 |
| w0013 | + | 0.00783 | 0.03168 | 0.04669 | 23 | 0.99982 |
| w0013 | - | 0.00786 | 0.03167 | 0.05577 | 2908 | 0.99985 |
| w0014 | + | 0.00772 | 0.03138 | 0.04437 | 1900 | 0.99985 |
| w0014 | - | 0.00789 | 0.03436 | 0.05461 | 2944 | 0.99985 |
| w0015 | + | 0.00798 | 0.03397 | 0.07753 | 2819 | 0.99984 |
| w0015 | - | 0.00777 | 0.03056 | 0.08645 | 1189 | 0.99986 |
| w0016 | + | 0.00789 | 0.03066 | 0.04272 | 2222 | 0.99975 |
| w0016 | - | 0.00763 | 0.03003 | 0.05755 | 1273 | 0.99980 |
| w0017 | + | 0.00775 | 0.02849 | 0.04517 | 767 | 0.99979 |
| w0017 | - | 0.00789 | 0.02874 | 0.05652 | 795 | 0.99974 |
| w0018 | + | 0.00801 | 0.03246 | 0.05388 | 1420 | 0.99983 |
| w0018 | - | 0.00772 | 0.03262 | 0.08149 | 920 | 0.99984 |
| w0019 | + | 0.00795 | 0.02941 | 0.05754 | 2369 | 0.99983 |
| w0019 | - | 0.00774 | 0.03109 | 0.05411 | 2427 | 0.99981 |

## Max |diff| of window/block features (over tracks)

| feature | max abs diff |
|---|---:|
| mean | 0.00036 |
| sd | 0.00038 |
| q50 | 0.00392 |
| q90 | 0.00428 |
| max16 | 0.00541 |
| blk64_mean_max | 0.00331 |
| blk64_top20_max | 0.00173 |
| blk128_mean_max | 0.00223 |
| blk128_top20_max | 0.00294 |
| blk256_mean_max | 0.00187 |
| blk256_top20_max | 0.00235 |
| blk512_mean_max | 0.00121 |
| blk512_top20_max | 0.00190 |
| blk1024_mean_max | 0.00064 |
| blk1024_top20_max | 0.00120 |
| blk2048_mean_max | 0.00046 |
| blk2048_top20_max | 0.00081 |

## Label-region mean

| record | family | orient | ref | cand | abs diff |
|---|---|---|---:|---:|---:|
| w0002 | CRISPR | + | 0.8216 | 0.8215 | 0.00009 |
| w0002 | CRISPR | - | 0.8116 | 0.8116 | 0.00000 |
| w0005 | ART | + | 0.4333 | 0.4331 | 0.00017 |
| w0005 | ART | - | 0.4141 | 0.4137 | 0.00035 |

## Label-blind ranking (window score = max over orientations)

Order flips beyond 0.01 tie tolerance: none
