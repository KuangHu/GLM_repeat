# Review packet — repeat explainer v0.2 (FROZEN)

Everything needed to review explainer v0.2 independently. **This is a review
artefact, not a release.** v0.2 is frozen exactly as reviewed; no parameter was
changed in response to the findings. Detector outputs and thresholds are untouched
and are not part of this packet except as frozen configuration.

Start with **`review/FINDINGS.md`**.

## Layout

```
code/      explain_repeats_v2.py   the reviewed explainer (v0.2)
           make_review_set.py      deduplication into repeat-pattern clusters
           explain_repeats.v0.1.py, make_review_set.v0.1.py   preserved predecessors
tests/     test_explainer_review.py  checks reproducing every finding
rules/     EXPLAINER_V2.md         model-selection and acceptance rules as implemented
data/      controls_definition.tsv  the 15 control regions (coordinates only)
           known_loci.tsv, annotated_units.tsv  annotation used for inspection only
           regions/*.fa            ORIGINAL region FASTAs: candidate span + 500-nt
                                   flanks, exactly the sequence the explainer saw
outputs/   controls/               models, copies, alignments, PWMs, unresolved,
                                   rejected cores for all 15 controls
           heldout_primary/        same for all 467 primary candidates (PWMs omitted
                                   for size; regenerate from code/)
           review_set/             the representative review set
environment/  pip_freeze.txt, python_version.txt, detector_config_frozen.md
review/    FINDINGS.md             the review write-up
MANIFEST.sha256                    checksums of every file in this packet
```

## How to re-run

```bash
python code/explain_repeats_v2.py \
  --candidates <annotated candidates tsv> --genomes <dir of .fna> \
  --out <outdir> --peaks <peaks tsv>
python -m pytest -q tests/
```

Genome FASTAs are not shipped (they are public NCBI records; accessions are in
`data/known_loci.tsv` and in every region header). The region FASTAs in
`data/regions/` are self-contained for inspecting any control or review-set entry.

## What this packet does and does not establish

- A model was emitted for 15 of 15 controls. **That is not 15/15 validated.**
  `review/FINDINGS.md` grades each: 13 well supported, 1 partial (VIPR phase), 1 weak
  (2 copies).
- Five implementation bugs are confirmed and reproduced by tests, including one found
  by a test written for this review (B5, period search).
- Six method limitations are quantified against random-sequence nulls.
- Nothing here is independent VIPR validation; that material does not exist yet.
- `unclassified` means no overlapping annotation in the records used. It is not a
  claim of novelty, function, or a new biological system.
- Alignment is ungapped throughout; indel-containing arrays are under-reported.

## Field meanings (see FINDINGS.md L1-L6 for the caveats)

| field | meaning | caveat |
|---|---|---|
| `median_unit_support` | 1 - mismatches/len vs the consensus | self-supporting; ~0.63 on noise at 2 copies |
| `information_bits_per_position` | mean per-column information of the PWM | rises with copy count alone |
| `phase_info_margin_bits` | information gap to the runner-up phase | confounded by unit count (B2) |
| `boundary_identity_left/right` | identity at the period just outside the interval | only `period` comparisons; very low power |
| `n_copies` | copies in span **plus flanks** | 51% of copies lie in flanks |
| `n_alternative_models` | whether a second branch fired | range {0,1}, misnamed (B4) |
