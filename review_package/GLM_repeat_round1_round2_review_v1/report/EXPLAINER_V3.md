# Explainer v0.3 evidence fixes

v0.3 fixes evidence preservation, candidate attribution, final-consensus distances
and reporting identified in the independent review of commit
`028c2926878c4f98cc1335cd803fba04431f0177`. Detector code, thresholds and results,
explainer v0.2, and the original review packet are preserved. This is a CPU-only
explanation revision, not a new detector benchmark.

## Results on the same inputs

All 15 controls and all 467 primary held-out candidates were processed uniformly.
Labels, known motifs and expected copy counts were not used by the explainer.

| Held-out status | Count | Meaning |
|---|---:|---|
| candidate_associated | 259 | At least one emitted model has a copy overlapping the candidate |
| flank_only | 84 | Models were found, but every copy is outside the candidate |
| unresolved | 124 | No supported seed/period interval model was emitted |

The old 343 `explained` rows now split into 259 candidate-associated and 84
flank-only cases. This is not a sensitivity decrease or an FDR estimate: the
status now describes evidence attribution. The independent v0.2 audit found 100
primary models with no overlap; v0.3 can find candidate-associated alternatives
for some of those cases because all proposed models are retained and checked.

There are 1,641 held-out hypotheses, with 230 candidates having more than one
candidate-associated hypothesis. `models.tsv` selects a display representative;
`hypotheses.tsv` and its linked evidence keep all alternatives. Neither a display
representative nor a candidate-associated status is biological validation.

## Core and extension contract

The motif branch retains the original 12-mer seed length, six-seed budget and
20% rounded mismatch rule. It proposes seed-supported cores, constructs the core
search consensus, and emits the entire resulting core coordinate set. Each copy
records its search motif, search distance and budget separately from its distance
to the final consensus.

Extensions are optional annotations on each core at fixed, linked coordinates.
The original 0.75 column-conservation and maximum 60-nt extension settings remain.
An extension template is the consensus of core-anchored contexts, rather than
the sequence of whichever core happened to occur first. Each context is checked
once against that template using the original length-scaled mismatch budget.
Extension attempts cannot create, replace or delete cores. Accepted and rejected
attempts both record the core copy ID and the actual decision evidence.

This intentionally changes v0.2 extension recruitment: a permissive extended
template no longer sweeps the whole region for extra unlinked copies. Additional
cores require independent core evidence. If no flank meets the column rule,
there is no extension attempt; the core remains in `copies.tsv`.

`copies.tsv` contains core units for motif models, even when all extensions are
accepted. Full extended-repeat sequences are in `extensions.tsv` and the
extension layer of `alignments.tsv`/`pwms.tsv`. Counts of these two layers must
not be added together as independent copies.

## Support and model choice

Every `mismatches_to_consensus` and `consensus_identity` is calculated against
the final saved consensus. The latter is an in-sample descriptive identity.
`loo_consensus_identity` and `loo_log2_enrichment` fit the profile using the other
copies; the enrichment uses the composition of the flanked input region as its
background. A quarter-count per base is used for both profiles and background.
These are descriptive evidence measures, not calibrated probabilities. Leaving
a copy out of its profile does not remove discovery/selection bias.

The old `n * unit_len * information` winner score is removed. Display ordering
prefers candidate association, then non-harmonic models, then candidate bases
covered, mean descriptive leave-one-out enrichment, and shorter units. All
models remain available. This ordering chooses an example to display; it does
not resolve biological ambiguity or make differently sized models statistically
equivalent.

Two copies are reported as `repeat_pair`, with spacing CV unavailable. One gap
cannot establish regularity. This does not dismiss strong evidence such as two
identical 39-nt segments.

## Period, phase and boundaries

All periods in 2–300, limited to at least three units in the candidate segment,
are recorded with identity, number of comparisons, local composition expectation
and excess identity. An offset is not promoted merely because it is the maximum
of a noisy search. The original 0.55 identity floor is used to propose periodic
intervals, and every supported interval is retained instead of choosing the
longest interval somewhere in the flanks.

Supported divisor periods whose intervals overlap are explicitly linked. A
non-harmonic alternative is displayed before a harmonic, but both remain in the
output. This removes the forced raw-identity argmax; it does not calibrate the
identity floor for composition or multiple testing. Long periods should not be
rejected merely for being long.

The union of supported local windows has its right endpoint computed without
the previous extra base. Interval boundaries still come from a heuristic sliding
identity rule and remain uncertain. All complete units starting at the interval
origin are displayed; trailing incomplete bases are not called copies.

`display_phase=0` is an interval-relative coordinate convention.
`phase_status=unanchored_cyclic` explicitly leaves the biological origin
unresolved. The old phase-information margin and weak outside-boundary identity
confidence fields are not emitted. Circularly equivalent column rotations cannot
establish a unique biological start.

Unknown bases split input into A/C/G/T segments; they are never replaced with A.
Coordinates and stored copy sequences are always genomic-forward, irrespective
of the biological RNA direction.

## Control evidence and corrections to the historical review

The primary model's units overlap all currently annotated units in each of the
eight CRISPR controls, but several controls also have additional core matches.
The 30-nt Pyrococcus and 36-nt Thermus extended repeats are in the extension
layer. A 12-nt core is not a complete repeat boundary.

| ART control | Primary cores | Annotated cores overlapped | Extra core matches | Accepted extensions |
|---|---:|---:|---:|---:|
| Madawaska | 5 | 5/5 | 0 | 5 of length 14 |
| LY01 | 8 | 5/5 | 3 | 0 |
| S6 | 6 | 5/5 | 1 | 5 of length 17 |
| PALS_2 | 6 | 5/5 | 1 | 5 of length 17 |
| UFV_DC4 | 5 | 5/5 | 0 | 5 of length 14 |

These are generated coordinate-overlap results, not assertions of exact boundary
recovery. Madawaska's five-core hypothesis was available among the seed proposals;
retaining proposals and changing the documented display rule exposes it without
hard-coding the ART motif or expected count. The historical v0.2 `FINDINGS.md`
incorrectly said all ART controls recovered all five cores and each had two
extras. Its actual saved v0.2 outputs were 3, 8, 7, 7 and 5 copies respectively.
Its pooled `13/15 well supported` statement should not be reused.

VIPR displays period 5, 16 units across the inferred interval, of which 12 are
fully contained in the labelled 60-nt span. Four units are outside that span.
The coincidence of the display convention with these label boundaries was not
an input or target. It is not evidence of a calibrated boundary or phase, and
does not validate all 16 units or establish independent VIPR generalization.

The off-label SUSP1 repeat pair retains its two cores and two identical 39-nt
extensions. Sequence repetition is supported; RNA function remains unknown.

## Annotation and copy identity

`add_cds_context.py` appends CDS overlap, union coverage, products, locations and
record availability without changing candidates, detector scores or original
selected RNA/repeat categories. `unclassified` means no selected RNA/repeat
category, not no GenBank annotation. Missing records yield unavailable fields,
not a zero-overlap assertion. Joined feature parts are evaluated separately.
The GenBank files used here were retrieved during the independent review; their
sequences match the frozen region FASTAs. They are not presented as the original
historical annotation snapshot.

`interval_id` is a shared contig/start/end key. Multiple hypothesis assignments
may reference it. Shifted overlapping intervals are linked in
`copy_relations.tsv` with unresolved identity; they are not automatically merged
into one biological copy. Counts of evidence rows or unique coordinate intervals
are not counts of unique biological repeats.

## Validation and review set

All 70 repository tests pass, including 17 new semantic regression tests. The
CPU-only 17-test suite runs independently with `requirements-review.txt`.
Original-DNA validation checks all 1,526 control and 9,313 held-out copy evidence
rows, including alternatives. All coordinates, consensus distances, association
flags and counts pass. The full core search sets were reconstructed for 84 control
and 1,464 held-out motif hypotheses. All 5,991 extension attempts have valid core
links and decision evidence; the 147 rejected attempts retain their cores.

The refreshed review set contains 48 distinct candidates: current representatives,
high-scoring unresolved/flank-only cases, every previous review representative
and the same fixed random10 (seed 20261015). Its 219 exact consensus navigation
groups use a different definition from v0.2's approximate 180 clusters and must
not be compared as numbers of biological systems.

The v0.2 explainer code and evidence remain preserved. Commit `0c2e15e` corrected
its packaging-only test entry point, removed the accidental detector test and
excluded caches from its manifest; its standalone explainer tests now pass 11/11.
All 46 manifest entries match their Git blobs, which contain mixed line endings.
The v0.3 package uses a complete `scripts/` + `tests/` layout, has its own CPU
requirements and an exact-archived-byte manifest excluding caches.

The published `releases/GLM_repeat_explainer_v0.3_review.zip` is self-contained:
it includes all original public DNA/GenBank records, pre-annotation candidates,
saved outputs, validation references, fixed review selections and 17 tests.
After extraction and CPU dependency installation, run
`python scripts/verify_explainer_packet.py --run-tests --rerun rerun` from the
packet directory to regenerate and compare all scientific evidence, including
the CDS context. No repository checkout, external sequences or GPU are required.
The archived output/provenance files remain frozen; relocated annotation-source
paths and regenerated source-code byte hashes are recorded separately from the
scientific-table comparison. The packet README gives the exact comparison rules.

## Reproduce locally

From the repository root, with the original public sequences available:

```text
python -m pip install -r requirements-review.txt
python -m pytest -q tests/test_explainer_v3.py
python scripts/add_cds_context.py --candidates results/pilot/heldout/annotated_primary/candidates_annotated.tsv --genbank independent_review/source/genomes.gb --out results/explainer_v3/inputs/heldout.tsv
python scripts/add_cds_context.py --candidates review_packet_explainer_v0.2/data/controls_definition.tsv --genbank independent_review/source/genomes.gb --out results/explainer_v3/inputs/controls.tsv
python scripts/explain_repeats_v3.py --candidates results/explainer_v3/inputs/heldout.tsv --genomes independent_review/source/genomes.fa --out results/explainer_v3/heldout_primary
python scripts/explain_repeats_v3.py --candidates results/explainer_v3/inputs/controls.tsv --genomes independent_review/source/genomes.fa --out results/explainer_v3/controls
python scripts/validate_explainer_v3.py --controls results/explainer_v3/controls --heldout results/explainer_v3/heldout_primary --genomes independent_review/source/genomes.fa --control-definitions review_packet_explainer_v0.2/data/controls_definition.tsv --annotated-units review_packet_explainer_v0.2/data/annotated_units.tsv --v2-heldout review_packet_explainer_v0.2/outputs/heldout_primary --out results/explainer_v3/validation
python scripts/make_review_set_v3.py --results results/explainer_v3/heldout_primary --genomes independent_review/source/genomes.fa --previous-review review_packet_explainer_v0.2/outputs/review_set/review_set.tsv --fixed-random independent_review/random10_fixed.tsv --out results/explainer_v3/review_set
```

Every output table has a header even when empty and is replaced on rerun. Missing
peak input yields blank peak-support fields; no GPU tracks were used in this
revision. Input checksums are saved per run.

Remaining limitations include ungapped alignment, exact-seed dependence in the
motif branch, the bounded period range, descriptive rather than calibrated model
evidence, and lack of independent natural VIPR validation. None is converted into
a claim of general completeness by this repair.
