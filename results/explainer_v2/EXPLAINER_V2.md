# Repeat explainer v0.2

v0.1 preserved at `provenance/explain_repeats.v0.1.py` (sha256 4d19a33e2887...) with
its outputs intact under `results/explainer/`. v0.2 applied uniformly to the
controls and to all 467 primary held-out candidates. No GPU inference, no new
background calibration, detector outputs and thresholds untouched.

Ungapped limitation unchanged and explicit: substitution distance only. Arrays whose
copies differ by indels are under-reported — copies may be split, truncated or
missed, and spacer lengths absorb the indel. Period estimates are likewise
substitution-based. Carried per model as `indel_handling`.

## Change 1: periodic profile as a first-class model source

The periodic branch needs **no exact k-mer seed**. Per candidate it searches periods
2-300, takes the period maximising identity at that offset, finds the maximal
interval where local identity stays above threshold, and chooses phase by PWM
information. Nothing about period, motif or boundary is hard-coded.

It emits the **same schema as the motif branch**: per-unit coordinates and
sequences, alignment, PWM position probabilities, per-unit support, plus
periodic-specific uncertainty — `periodic_phase`, `phase_info_margin_bits` (how much
better the chosen phase is than the runner-up) and `boundary_identity_left/right`
(identity at the period just outside the called interval).

Result on the VIPR control, which v0.1 could not model: the labelled 60-nt YNNGG
region is now explained as **period 5, 15 units**, discovered de novo. Its interval
extends slightly past the labelled 60 nt, which is exactly what the boundary fields
are there to expose. No YNNGG, no period 5 and no labelled VIPR boundary appears
anywhere in the code.

**This is sequence explanation only.** It is kept separate from independent VIPR
validation, which remains unavailable: no independent natural VIPR locus was
obtained, so VIPR generalization is still unresolved. Explaining the sequence of a
development locus is not evidence that the detector generalizes to VIPR.

## Change 2: core vs extended repeat evidence

Core and extension are now separate. The seed-supported core motif and its copies
are kept regardless of whether flank extension succeeds; the extended motif is an
optional refinement that never discards core copies. `core_copies` and
`extended_copies` are reported separately, and every core copy the extended motif
rejects is written to `rejected_cores.tsv` with its core mismatch count, its
mismatch count against the extended consensus, and the budget it exceeded.

Effect on the ART controls, where v0.1 over-extended and lost copies:

| locus | v0.1 copies | v0.2 source | v0.2 core | v0.2 extended |
|---|---:|---|---:|---:|
| LY01 | 3 | motif_core | **8** | 0 |
| S6 | 5 | motif_core | **7** | 0 |
| PALS_2 | 5 | motif_core | **7** | 0 |
| UFV_DC4 | 5 | motif_extended | 5 | 5 |
| Madawaska | 3 | motif_extended | 3 | 3 |

For LY01, S6 and PALS_2 the extension matched nothing and the core set is now
preserved instead of being discarded. Note these core counts (7-8) exceed the five
annotated cores: the retained core is a 12-mer with a 2-mismatch budget, which
matches more loosely than the manuscript's 14-nt core. That is a property of the
de novo core, not a corrected count of annotated cores. Madawaska still recovers
only 3 copies; its seed k-mer is the limiting factor, not extension.

Two Pyrococcus arrays report `motif_core` with 52 and 22 copies against 50 and 20
extended, with 2 rejected cores each — recorded with evidence rather than silently
dropped. No motif sequence or core length is hard-coded anywhere.

## Applied to all 467 candidates

| | v0.1 | v0.2 |
|---|---:|---:|
| explained | 337 | 343 |
| unresolved | 130 | 124 |

Model source in v0.2: motif_extended 164, motif_core 163, periodic 16.
Architectures: dispersed 201, regular-spaced 66, variable-spacer 56, tandem 14,
short-period tandem 6.

Unresolved cases are not forced into a model: 124 candidates still have none, and
they are emitted with their alternatives.

## Review set, refreshed

343 explained candidates deduplicate into **180 repeat-pattern clusters**. A
repeat-pattern cluster is a group of candidates whose consensus sequences resemble
one another. It is **not** a biological family, system or classification; sequence
similarity alone does not establish shared origin or function.

The refreshed review set holds 15 interpretable `unclassified` representatives
spanning 5 architectures plus 8 high-scoring candidates with no sequence
explanation. `unclassified` means the span carries no overlapping annotation in the
GenBank records used here, and nothing more.
