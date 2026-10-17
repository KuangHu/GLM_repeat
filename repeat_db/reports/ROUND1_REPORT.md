# Exploration round 1: host + phage, frozen configuration

Config frozen throughout: window 3000 / stride 1500, random_partition 20x2,
scoring seed 20261014, attention-mask-policy always, both orientations,
model 08754cba59a1f97d517f873fad6c672d2b1abdc7, threshold file
results/pilot/calibration/thresholds/thresholds.json key bg100_per_Mb (tau 0.80).
`bg100_per_Mb` is a historical configuration name, NOT a density or FDR guarantee
for this round. Scores were computed once over the full grid (`--audit-fraction 1`);
the CDS policy was applied only as a counterfactual at candidate-calling time.

## Completeness

| set | windows | both orientations | incomplete | scoring settings observed |
|---|---:|---:|---:|---|
| host | 3,061 | 3,061 | 0 | one frozen setting |
| phage | 2,996 | 2,996 | 0 | one frozen setting |

No tracks for unknown windows; no shard retried or double-counted.

## Yield and cost, reported separately per data type

| | host round-1 | phage round-1 |
|---|---:|---:|
| genomes | 4 (4 lineages) | 72 |
| taxonomy (verified against manifest) | 4 lineages | 72 distinct named families, 11 branches, one genome per family |
| unique scanned | 4.616 Mb | 4.551 Mb |
| tracks | 6,122 | 5,992 |
| GPU hours (measured) | 18.38 | 17.90 |
| **GPU-h per Mb** | **3.98** | **3.93** |
| seconds per track | 10.8 | 10.8 |
| candidates | 170 | 399 |
| **candidates per scanned Mb** | **36.8** | **87.7** |
| target-region density (denominator) | 117.4/Mb (1.082 Mb) | 236.0/Mb (1.373 Mb) |
| target-associated candidates | 127/170 | 324/399 |

Phage yield is 2.4x the host yield on near-identical scanned length and identical
GPU cost. These must not be averaged. Per-genome phage density spans 0 to 343/Mb
(median 82/Mb); 5 of 72 phage genomes produced no candidate, which is a selection
and threshold outcome, not evidence of absence.

The earlier 143/Mb figure from the single-genome smoke test (SUSP1) did not
transfer to either panel.

## CDS-policy audit (counterfactual only; policy was NOT used to scan)

| | host | phage |
|---|---:|---:|
| window nt reduction offered by policy | 11.5% | 14.0% |
| candidate envelopes lost | 6/170 (3.5%) | 17/399 (4.3%) |
| copy records not retained | 287/5,837 -> 287/2,079 (13.8%) | 392/5,837 (6.7%) |
| of those, no overlapping policy copy at any threshold | 271/287 | 389/392 |
| attributable to windows never scored | 287/287 | 392/392 |
| attributable to model change with sequence available | 0 | 0 |

Copy relation is deduplicated by genomic coordinate and reported as
retained / lost / gained with coordinates in each run's
`policy_attribution/copy_relation.tsv`.

Host: retained 1,792, full-only 287, policy-only 15, net -272.
Phage: retained 5,445, full-only 392, policy-only 32, net -360.

> Under the CDS policy, 13.8% (host) and 6.7% (phage) of full-grid explainer copy
> records were not retained. The large majority have no overlapping policy copy at
> any threshold and are attributable to supporting windows never being scored; a
> small minority are boundary-shifted re-representations. These proportions are not
> verified biological miss rates.

Focus case `NC_000908.2:485272-487045` (score 4.26, near the top of the host panel)
shows both effects at once: 2 of its 4 supporting windows are policy-excluded, and
of 27 "lost" copies, 12 match a policy copy at >=50% reciprocal overlap as a 4-nt
phase shift while 15 have no counterpart. Envelope overlap did not imply core
retention. Full detail in `region_round1/policy_attribution/FINDINGS.md`.

**Decision: full grid is the default for BOTH host and phage routes.** Window
saving is limited on both (11.5% and 14.0%) and both show copy-evidence loss. CDS
is retained for neighbourhood interpretation and as a counterfactual audit, not as
a compute saving. Whether the filter can be repaired into a reliable one is
separate work.

### Loss attribution on two independent axes

Upstream trigger (which stage changed) and representation outcome (whether the copy
still appears in any form) are different questions and are reported separately.

| upstream trigger | host | phage |
|---|---:|---:|
| supporting window excluded | 270 | 351 |
| window excluded, copy still inside policy explainer input | 17 | 41 |
| candidate range change alone | 0 | 0 |
| model change with sequence available to both | 0 | 0 |

| representation outcome | host | phage |
|---|---:|---:|
| no corresponding copy | 271 | 389 |
| phase / boundary shift (>=50% reciprocal) | 12 | 0 |
| partial overlap only | 4 | 3 |

"Window exclusion is the upstream trigger for all of them" does NOT mean nothing
else changed: 12 host copies reappear with shifted boundaries from a different
window. Only `no_corresponding_copy` means the evidence is absent from the policy
output, and even that is a computational evidence difference, not a verified
biological miss.

## Explainer outcome (v0.3, all hypotheses retained)

| | host | phage |
|---|---:|---:|
| candidate-associated | 74/170 | 215/399 |
| flank-only | 18/170 | 36/399 |
| unresolved | 78/170 | 148/399 |

On the unbiased random stratum (40 candidates each, the only stratum that supports
a yield estimate): host 17 candidate-associated / 7 flank-only / 16 no model;
phage 21 / 3 / 16. Architectures are dominated by `repeat_pair` and
`dispersed_repeat`, i.e. weak two-copy evidence, in both panels. A `repeat_pair`
does not by itself establish functional significance.

### No-model candidates are not simply false positives

All 226 no-model candidates were re-probed with an indel-tolerant method
(exact 10-mer seeds + affine-gap local extension). A minimap2 `-x asm20` probe was
tried first and DISCARDED: its positive control found nothing for candidates with
11-34 verified copies, because minimizer seeding cannot seed on 8-33 nt units. The
retained probe reports 90% (host) and 85% (phage) sensitivity on modelled
candidates, so its negatives are interpretable.

| probe verdict | host (78) | phage (148) |
|---|---:|---:|
| seeds present but no extension | 50 | 89 |
| ungapped similarity the explainer did not model | 21 | 31 |
| no repeated seed at all | 6 | 26 |
| gapped similarity the ungapped explainer cannot represent | 1 | 2 |

So `no_model` splits into at least three different situations: a small indel
coverage gap (3 of 226), a larger group with detectable ungapped similarity that
the explainer nevertheless did not model (52 of 226), and weak or absent seed
evidence (171 of 226). None of these is established as a detector false positive
without human review.

## Manual inspection sheets

`manual_review/manual_review.tsv` per run: host 65 rows, phage 77 rows. Strata are
labelled and reported separately: 40 uniform-random each (unbiased), 10 high,
10 boundary-adjacent, 12 structure representatives, and ALL policy-lost candidates
(6 host, 17 phage) with no subsampling. Verdict columns are empty and awaiting
human review.

## Limitations

- Exploration round. No recall, no FDR, no calibrated biological probability.
- The frozen configuration missed the tested annotated VIPR locus
  (NC_028808.2:[29169,29229) produced no candidate at tau 0.80). VIPR detection
  capability is not established.
- The whole-virus metadata query failed; the catalogue was assembled from 18
  verified branches. Branch coverage is enumerated, not claimed exhaustive.
- Cross-circular-origin detection is not implemented. This matters more for phages
  than for host chromosomes and is an uncovered region of every circular genome.
- `host` in virus metadata is the isolation source, not a verified replication
  host; domain was taken from viral taxonomy instead.
- Backgrounds and calibration genomes are unlabelled real DNA, never verified
  negatives.
- pretraining_overlap_status is `unknown` for every sequence.

## Next round proposal (NOT started)

Budget from measured cost, ~3.95 GPU-h per Mb, and from each panel's own yield:

| composition | Mb | GPU-h | expected candidates |
|---|---:|---:|---:|
| host, new lineages | 30 | ~119 | ~1,100 |
| phage, new families | 20 | ~79 | ~1,750 |
| **total** | **50** | **~198** | **~2,850** |

Host and phage budgets are kept separate because their yields differ 2.4x. The
~2,850 candidate estimate implies the review strategy must stay stratified
sampling plus full inspection of policy-lost cases; exhaustive manual review is not
feasible at this scale. Recommend proceeding only after the round-1 manual sheets
are reviewed, since the random stratum suggests ~40% of candidates carry no
sequence-level model at all.
