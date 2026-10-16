# Detector pilot closeout

Frozen scorer unchanged throughout. No detector retuning, no classifier training,
no Prodigal exclusion, no new genome-wide inference, no expanded background
calibration.

## 1. Exact audit of 24 held-out windows

48 tracks, original untrimmed 3,000-nt contexts. Peaks called on both modes at a
matched tau = 0.80 **only** to make structure comparable; that is not a calibrated
exact operating point, and exact is not treated as ground truth.

| stratum | tracks | MAE | Pearson | Pearson 16-nt | peak-span Jaccard | exact has cluster | RP has cluster | both | neither |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| high | 16 | 0.038 | 0.994 | 0.999 | 0.993 | 16 | 16 | 16 | 0 |
| adjacent | 16 | 0.071 | 0.984 | 0.984 | 0.590 | 6 | 8 | 5 | 7 |
| non_hit | 16 | 0.064 | 0.982 | 0.983 | 1.000 | 0 | 0 | 0 | 16 |

- **Strong detections reproduce exactly**: 16/16 agreement, cluster-span Jaccard 1.000.
- **Non-hits reproduce**: no window called by exact that RP missed in this sample.
- **Marginal detections do not always reproduce**: 4 of 48 tracks disagree, all in
  the `adjacent` stratum. Three are RP-calls / exact-does-not (RP scores 1.305-1.335,
  at the 2-peak floor of 1.268); one is exact-calls / RP-does-not (1.689).

Marginal RP calls near the cluster-formation boundary are therefore not stable
under a change of masking scheme. This bounds how much weight a single
low-scoring candidate can carry. It does not make exact correct and RP wrong:
exact is a different masking scheme of the same model, not an independent
observation.

Sampling-rule caveat: the `adjacent` stratum is the 8 windows with the **lowest
surviving cluster score**, not the 8 closest in value to tau. Sample preserved as
drawn; description corrected (results/pilot/heldout/exact_audit/sampling.json).

## 2. Terminology corrections (coordinate-level)

`results/pilot/MERGING_AUDIT.md` section 5. Withdrawn: "~230-nt gaps matching
170-215-nt spacing". It compared edge gaps to start-to-start pitch and quoted
maxima as typical values.

Reconciled: CRISPR peak pitch equals annotated repeat pitch (67-68 / 75-76 nt) with
**0 of 230 annotated units missing** across 8 arrays; ART peak pitch runs
consistently below annotated pitch (171-200 vs 203-204 nt), 25/25 cores matched.
Extra peaks (0-2 per array) lie inside array spans on no annotated unit; one peak
merges multiple copies. The ART pitch shortfall is real and unexplained here.

## 3. Explainer module

Separate module, detector outputs read-only, no gLM2 inference. See
`results/explainer/EXPLAINER_SUMMARY.md`.

- Controls, de novo with no motifs supplied: CRISPR repeats recovered at correct
  pitch; every ART consensus contains `ATATGAATACGTAT`. Over-extension cost copies
  in two phages (3 of 5 cores).
- All 467 primary candidates processed: 337 explained, 130 unresolved. Nothing forced.
- Architectures: dispersed 181, regular-spaced 100, variable-spacer 52, tandem 4.
- Emitted: consensus, per-copy coordinates with mismatches, alignments, PWM position
  probabilities, spacing and dispersion, peak-to-copy correspondence, alternatives,
  unresolved cases.
- Indel limitation recorded per model: ungapped substitution alignment only.

**VIPR, kept separate from detector recovery:** the labelled 60-nt YNNGG region is
the one control the motif explainer could not model (a degenerate 5-nt pattern
shares no exact 12-mer seed). Its alternative fired de novo: periodic self-match,
period 5, identity 0.673. Explainer outcome is a periodic profile only, with no
copy-level model. VIPR generalization remains unresolved.

## 4. Review set

337 explained candidates deduplicate to 185 families. Review set: 11 interpretable
`unclassified` representatives spanning 4 architectures plus 8 high-scoring
candidates with no sequence explanation.

`unclassified` means the span carries no overlapping annotation in the GenBank
records used here. It is not a claim of novelty, function or a new biological
system.

## Open items for any future round

- ART peak pitch shortfall (171-200 vs 203-204 nt) is unexplained.
- Explainer over-extends motifs and loses copies; no indel support.
- Exact has no calibrated operating point at any density (declined: ~150 GPU-h).
- VIPR has no independent held-out locus; generalization untested.
