# CDS-policy audit, host round-1 (4.616 Mb, 4 lineages)

Scores computed ONCE over the full grid (`--audit-fraction 1`); candidates then
called separately from full-grid windows and from CDS-policy windows. All
hypotheses are compared, not only the primary display model (copy evidence covers
459 full / 408 policy hypotheses).

## Copy-set relation, deduplicated by genomic coordinate

| relation | copies |
|---|---:|
| retained (exact coordinates) | 1,792 |
| present only in full grid | 287 |
| present only under policy | 15 |
| net change | **-272** |

287 is the number of full-grid copy records not retained, not a net loss, and not a
verified biological miss rate.

## How much of the 287 is genuine absence vs re-representation

Interval matching against every policy copy, at explicit reciprocal-overlap
thresholds:

| reciprocal overlap | lost copies with a policy match |
|---:|---:|
| >= 0.9 | 0 |
| >= 0.7 | 8 |
| >= 0.5 | 12 |
| >= 0.1 | 16 |
| any overlap at all | 16 |

**271 of 287 lost copies have no overlapping policy copy at all.** The 16 that do
are concentrated in a single candidate (below); every other affected candidate has
zero matches at any threshold.

A 90% reciprocal threshold alone is misleading here: copies are often 12-13 nt, so
a 4-nt boundary shift falls to 0.67 and is scored as lost. Both the strict and the
permissive views are reported above for that reason.

## Attribution of the 287 by cause

| cause | copies |
|---|---:|
| supporting window excluded by policy | 270 |
| window excluded, copy still inside the policy explainer input range | 17 |
| outside explainer input (boundary move only) | 0 |
| model change with sequence available to both | 0 |

No copy loss is attributable to the explainer switching representative models while
the same sequence was available to both runs.

## Per candidate

| full candidate | score | envelope | windows excluded | copies lost | lost with >=50% policy match |
|---|---:|---|---|---:|---:|
| NC_001318.1:844629-845024 | 2.38 | lost | 2/3 | 107 | 0 |
| NC_001318.1:615756-615807 | 1.77 | lost | 2/3 | 59 | 0 |
| NC_001318.1:677841-677874 | 1.42 | lost | 2/3 | 37 | 0 |
| NC_002607.1:204659-204728 | 1.47 | lost | 1/2 | 28 | 0 |
| **NC_000908.2:485272-487045** | **4.26** | **overlapped but shrank** | **2/4** | **27** | **12** |
| NC_000963.1:119257-119900 | 2.40 | lost | 3/3 | 12 | 0 |
| NC_002607.1:493043-493202 | 1.40 | lost | 2/3 | 12 | 0 |
| NC_001318.1:771842-771933 | 2.08 | overlapped but shrank | 1/2 | 5 | 0 |

## Focus case: NC_000908.2:485272-487045 (Mycoplasmoides, score 4.26)

Near the top of the 170-candidate panel, and the clearest illustration that
envelope overlap does not imply core retention.

Supporting windows and their policy status:

```
483000-486000  policy=True   target_overlap=769
484500-487500  policy=False  target_overlap=0
486000-489000  policy=False  target_overlap=0
487500-490500  policy=True   target_overlap=436
```

Copy sets: full 43, policy 41 nearby, retained 16, lost 27, **gained 25**.

The lost and gained coordinates are largely the same physical repeats with shifted
unit boundaries:

```
lost:   (485389,485401) (485390,485402) (485391,485403)
gained: (485385,485397) (485386,485398) (485387,485399)
```

A 4-nt phase shift. 12 of the 27 "lost" copies here match a policy copy at >= 50%
reciprocal overlap. The remaining 15 have no counterpart. So this candidate shows
both effects at once: real loss of sequence coverage AND re-representation of
retained repeats under different unit boundaries.

## Statement for the report

> Under the CDS policy, 13.8% of full-grid explainer copy records were not
> retained. 271 of 287 have no overlapping policy copy at any threshold and are
> attributable to supporting windows never being scored; 16, concentrated in one
> candidate, are boundary-shifted re-representations. This proportion is not a
> verified biological miss rate.

## Decision support

Window saving under the policy was 11.5% on this panel. The cost is 3.5% of
candidates and ~13% of copy records, including a top-scoring locus whose envelope
survived. Full-grid scanning is the appropriate default for the host route; the CDS
policy is retained as a post-hoc annotation layer and as a counterfactual audit.
