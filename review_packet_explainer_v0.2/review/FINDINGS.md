# Explainer v0.2 review findings

Written **before** any algorithm parameter was changed. v0.2 is frozen as reviewed;
nothing in this round was tuned in response to these findings. Detector outputs and
thresholds are untouched. No GPU inference was used for this review.

Findings are separated into implementation bugs (wrong code, fixable), method
limitations (the approach behaving as designed, but not meaning what a reader would
assume), and ambiguous biology (not resolvable from sequence here).

---

## Headline correction: "15/15 controls explained" is not "15/15 validated"

v0.2 emitted a model for all 15 controls. Emitting a model is not validation. Graded
by manual inspection:

| control | verdict | basis |
|---|---|---|
| 6 Pyrococcus CRISPR | **strong** | 30-nt unit, 20-52 copies, pitch 67-68 matching annotation; consensus is the real repeat |
| 2 Thermus CRISPR | **strong** | 36-nt unit, pitch 75-76 matching annotation |
| 5 ART phages | **good, with a threshold caveat** | all 5 annotated cores matched in each; 2 extra copies per locus (see A1) |
| SUSP1 VIPR 60-nt | **partial** | period recovered, phase not supported (B2, A2) |
| SUSP1 off-label | **weak** | rests on 2 copies, where self-support is ~0.63 on pure noise (L1) |

So: 13 of 15 are well supported, 1 partial, 1 weak. No control is independent
validation of the detector, and none is evidence about VIPR generalization.

---

## Implementation bugs (confirmed, reproducible)

### B1. Periodic alternatives are emitted without an identity gate (regression from v0.1)

v0.1 appended a `periodic_self_match` alternative only when identity at the period
was >= 0.5. v0.2 appends it unconditionally whenever a period is returned. Random
DNA has ~0.25 identity at any offset, so most of these are noise.

- 110 of 124 unresolved rows carry a periodic alternative below 0.5 identity
- 293 of 343 explained rows carry one too

Effect: a reader sees "period 12, identity 0.36" and may read it as a finding. It is
the argmax of a noisy search. **Fix: restore the >= 0.5 gate (or a per-length null).**

### B2. Phase scoring is confounded by unit count, so phase choice is not meaningful

Phases are scored by PWM information, but the number of units differs by phase
(`(len(region) - phase) // d`), and information rises with copy count (L2). On a
perfectly periodic sequence, where all phases are circularly equivalent by
construction:

```
ACGTT x 20, period 5:
  phase 0: 20 units, info 1.7211   <- always wins
  phase 1: 19 units, info 1.7099
  phase 2: 19 units, info 1.7099
  phase 3: 19 units, info 1.7099
  phase 4: 19 units, info 1.7099
```

The winner is the phase that fits one more unit, not the phase with better-conserved
columns. The VIPR control's reported `phase_info_margin_bits` is **0.0035**, smaller
than this 0.0112 artifact, so its reported phase is unsupported.

**Fix: score all phases over a common unit count and a common region.**

### B3. The same physical repeat copy is emitted under more than one candidate

Candidate spans are disjoint after merging, but the explainer adds 500-nt flanks, and
**73 candidate pairs have overlapping flanked regions**. 51% of all emitted copies
(1259 of 2466) lie in flanks rather than in the candidate span.

Result: **20 copy pairs overlap physically but are emitted under different
`candidate_id`s**, typically offset by ~3 nt, e.g.

```
NC_006461.1 [530488,530500) from candidate 530965-531005
NC_006461.1 [530491,530503) from candidate 530404-530430
```

Exact-coordinate deduplication finds nothing (0 identical triples), because the two
candidates infer slightly different unit boundaries. Any downstream count of
"distinct repeat copies" is therefore inflated.

**Fix: deduplicate copies globally by genomic overlap, not by exact coordinates.**

Coordinate round-tripping itself is **correct**: all 263 control copies and all 2466
held-out copies reproduce their stored sequence when re-sliced from the source
genome (0 mismatches).

### B5. The period search is dragged to large offsets and ignores harmonics

`best_period` returns the offset `d` maximising identity, but identity at offset `d`
is computed over only `len(seq) - d` comparisons. Larger offsets therefore have
noisier estimates and win by variance. There is no length normalisation and no null
correction.

On **200 pure-random 1000-nt sequences** (no repeat at all):

| quantity | value |
|---|---:|
| median period picked | **200** (search range 2-300) |
| fraction picking period > 100 | **0.82** |
| mean identity at the pick | 0.294 (chance is ~0.25) |

It also has no harmonic handling: multiples of the true period score as well as the
period itself. On a clean period-5 pattern (`CNNGG` x 60):

```
offset   5: identity 0.692 over 295 comparisons   <- true period
offset  10: identity 0.697 over 290 comparisons   <- harmonic, scores HIGHER
offset  95: identity 0.683 over 205 comparisons
```

A unit test written for this review asserted that a clean period-5 pattern returns
d = 5; it returned **d = 95**. In the shipped run, 61 of 467 reported periods exceed
100 and are most likely this artifact.

The VIPR control is not affected in practice: its span is 60 nt, so the search range
collapses to 2-20 and the true period wins. Long candidates are affected.

**Fix: normalise identity by comparison count (or score against a per-offset null),
and prefer the smallest offset among harmonics.**

### B4. `n_alternative_models` is misnamed

Its range is {0, 1} (320 and 23 of 343). It records whether a second *branch*
(motif vs periodic) produced anything, not a count of alternative models. Rename or
redefine.

---

## Method limitations (working as designed; the field names oversell)

### L1. `median_unit_support` is self-supporting and inflated at low copy counts

Support is measured against a consensus the copy itself helped define. Measured on
**pure random sequence with no repeat at all** (300 trials, 5-mers):

| copies | reported self-support | leave-one-out support |
|---:|---:|---:|
| 2 | **0.625** | 0.248 |
| 3 | 0.562 | 0.236 |
| 5 | 0.493 | 0.231 |
| 10 | 0.419 | 0.262 |

Leave-one-out sits at chance (~0.25), as it should. So a reported support of 0.6 on a
2-copy model is what noise produces. **97 of 343 models (28%) rest on <= 3 copies**
(47 on exactly 2). Those support values carry essentially no evidence.

Mitigation (not applied this round): report leave-one-out support instead, and
require a minimum copy count before quoting support at all.

### L2. `information_bits_per_position` conflates conservation with copy count

Identical copies, varying only n: info = 1.01 (n=3), 1.26 (n=5), 1.53 (n=10),
1.72 (n=20), 1.87 (n=50). The quantity is not comparable between models with
different copy counts, and it is not a measure of conservation alone.

### L3. Model selection is biased toward long periods with few copies

Selection uses `n_copies * unit_len * information_bits`. On **random sequence**, at
equal explained length (60 nt):

| model | selection score |
|---|---:|
| period 2 x 30 copies | 2.3 |
| period 5 x 12 copies | 7.9 |
| period 10 x 6 copies | 17.1 |
| period 30 x 2 copies | **27.2** |

A 12x bias toward the long-unit/few-copy end, on noise. This is why the motif branch
usually beats the periodic branch when both fire, independently of which is right.

### L4. Ungapped alignment only

Substitution distance throughout, as documented. Indel-containing arrays are
under-reported: copies split, truncated or missed; spacer lengths absorb the indel.

### L5. Copy counts are not a property of the candidate

51% of copies lie in the 500-nt flanks. `n_copies` describes the flanked region, not
the detector candidate. This is intended (discovery is deliberately not restricted to
peak intervals) but the field name invites the wrong reading.

### L6. Boundary identity is very low power

`boundary_identity_left/right` is computed over exactly `period` positions — 5 for the
VIPR control. The reported 0.0/0.0 is entirely consistent with chance (0.75^5 = 0.24
per side). It should not be read as a sharp boundary.

---

## Ambiguous biology (not resolvable from sequence here)

### A1. ART: 7 modelled copies vs 5 annotated cores

The 2 extra copies per ART locus are real sequence matches, not errors:

```
S6: 166620 TATTATATGAAT support 0.917   (between annotated cores 1 and 2)
    167463 CATTATAATAAT support 0.833   (downstream of core 5)
```

All 5 annotated cores are matched in every ART control. The extras arise because the
de novo core is a 12-mer with a 2-mismatch budget, which is looser than the
manuscript's 14-nt <= 1-mismatch rule. This is a threshold difference, not a
contradiction, and it is not evidence that the annotation is incomplete.

### A2. VIPR: period right, phase unsupported, and only ~3 informative positions

12 of the 15 modelled units fall inside the labelled 60 nt (12 x 5 = 60, matching the
annotated unit count). But the unit boundaries are offset from the annotated YNNGG
boundaries, and the phase is unsupported (B2).

Per-position information of the VIPR PWM: **[0.17, 1.11, 1.38, 1.11, 0.15]**. Two of
five positions are near-uninformative. The model is effectively a 3-nt conserved core
with 2 free positions, reported as a 5-nt unit, and `median_unit_support` averages
over the uninformative positions.

**This is sequence explanation of a development locus. It is not independent VIPR
validation, which remains unavailable.**

### A3. AT-rich consensus in phage genomes

Several models in AT-rich phage genomes have consensus such as `TTCTTTTTTAATT` and
`CATTATTATTACC` with 2-4 copies. These may be compositional artifacts rather than
discrete repeat elements. Composition-matched null models would be needed to tell.

---

## What the unresolved cases look like

The 8 highest-scoring unexplained candidates all have best-period identity 0.36-0.41
(near the 0.25 chance level), and three of them report periods of 70, 101 and 112,
which are the B5 artifact and dominant-base fractions 0.34-0.45. They are
correctly left unresolved. Their reported periodic alternatives are the B1 artifact
and should be ignored.

---

## Recommended fixes, in priority order (none applied this round)

1. B5 normalise the period search by comparison count and collapse harmonics; this
   also removes most spurious long periods.
2. B1 restore the periodic-alternative identity gate.
3. L1 replace self-support with leave-one-out support; suppress support below a
   minimum copy count.
4. B2 score phases over a common unit count.
5. L3 replace the selection score with one that is calibrated per unit length, or
   compare models against a composition-matched null.
6. B3 deduplicate copies globally by genomic overlap.
7. B4 rename `n_alternative_models`; L2/L5/L6 rename or document the fields in place.

After fixes, rerun controls and all 467 candidates uniformly and regenerate the
review set.
