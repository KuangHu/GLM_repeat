# GLM_repeat — gLM2 masked-predictability repeat detection

Detecting repetitive sequence organization in DNA from a genomic language model
(`tattabio/gLM2_650M`), without motifs, protein anchors or family labels.

**Status: detector pilot closed and frozen; explainer v0.3 evidence fixes available.**
The v0.2 code and results remain frozen; its packet received a packaging-only
test/manifest fix in `0c2e15e`. v0.3 keeps core copies and
optional extensions separately, distinguishes candidate-associated models from
flank-only findings, and retains alternative hypotheses with descriptive support.
It does not establish general recall, a false-discovery rate, or VIPR generalization.
This repository holds code, frozen configuration and reports. Model weights and
per-base score tracks are not included. The standalone v0.3 review ZIP includes
the archived public genome sequences needed to regenerate explainer evidence.

## Start here

Optional CDS-guided scanning and post-scan gene neighborhoods are available in
[REGION_SCAN_zh.md](REGION_SCAN_zh.md), including cluster commands, audit controls,
CPU validation and limitations. This is a new selectable workflow; the frozen
full-genome pilot remains unchanged. CDS filtering is heuristic, and does not
establish noncoding identity or complete repeat recall.
The downloadable [cluster source bundle](releases/README.md) includes the code,
tests and runbook without model weights or raw genomes.
The separate [standalone explainer v0.3 review packet](releases/GLM_repeat_explainer_v0.3_review.zip)
includes original DNA, GenBank annotations, frozen inputs/outputs and CPU tests.
Its verifier can rebuild annotations and all explainer evidence without this
repository, genome downloads, detector inference or a GPU.

| | |
|---|---|
| `results/explainer_v3/EXPLAINER_V3.md` | current behavior, evidence contract, limits and reproducible commands |
| `results/explainer_v3/validation/validation.json` | original-DNA checks for every v0.3 core/periodic copy and extension |
| `results/explainer_v3/validation/controls_overlap.tsv` | actual control recovery; annotation overlap is not full validation |
| `review_packet_explainer_v0.2/review/FINDINGS.md` | historical review; its ART/control summary is corrected in the v0.3 notes |
| `results/PILOT_CLOSEOUT.md` | what the detector pilot established, and what it did not |
| `results/pilot/MERGING_AUDIT.md` | candidate geometry, merging audit and coordinate reconciliation |
| `WORKFLOW_SETUP.md` | the original workflow specification |

## Method, briefly

A nucleotide is masked one at a time and scored as log P(true base) over {a,c,g,t}
minus log(1/4) — "nats above uniform", capped at ln 4 = 1.386. Exact scoring needs
one forward pass per base; a random-partition approximation (20 groups x 2 repeats)
costs ~1/75 of that and was validated against exact on two independent panels
before being adopted as the screening path.

A frozen, label-blind peak-cluster scorer calls peaks on the smoothed track against
a background-calibrated threshold, groups peaks within 500 nt, and requires >= 2
peaks. Thresholds are calibrated on background material only, to a target candidate
yield per unique scanned Mb, and frozen before any held-out scan.

## Pilot result

A 10-Mb pilot split by genome (no shared genus between calibration and held-out)
recovered **13 of 13 held-out loci** — 8 CRISPR arrays from *Pyrococcus furiosus*
(Archaea) and *Thermus thermophilus*, 5 ART arrays from *Staphylococcus* phages — at
a yield of 91.9 candidates per unique scanned Mb. CRISPR peak pitch matched the
annotated repeat pitch (67-76 nt) with 0 of 230 annotated units missed.

## What this does not establish

- Candidate yields are **yields, not false-discovery rates**. Background material is
  unlabelled real DNA, not verified negatives.
- VIPR generalization is **unresolved**: no independent natural VIPR locus was
  obtainable, so VIPR is untested outside development data.
- Scores measure masked predictability, not repeat probability, and not function.
- `unclassified` means no matching selected RNA/repeat category. The original
  classifier did not inspect CDS annotations. v0.3 appends explicit CDS overlap
  and annotation availability without filtering or reranking detector candidates.
- Explainer alignment is ungapped; indel-containing arrays are under-reported.

Corrections made during development are kept in place rather than edited away — see
the CORRECTION sections in `results/pilot/calibration/thresholds/FROZEN.md`,
`results/pilot/MERGING_AUDIT.md` and `results/pilot/heldout/exact_audit/report/`.

## Layout

```
scripts/     scoring, windowing, calibration, scanning, audit, explainer
slurm/       Slurm job scripts (LBNL Lawrencium, A40 GPUs)
tests/       CPU tests, no GPU required
config/      workflow defaults
provenance/  preserved earlier versions of code
results/     frozen reports, thresholds, audits (no score tracks)
review_packet_explainer_v0.2/   self-contained explainer review packet
```

Run the repository tests with `python -m pytest -q tests/` (requires PyTorch for
the original scorer tests). For CPU-only v0.3 review, install
`requirements-review.txt` and run `python -m pytest -q tests/test_explainer_v3.py`.
Use `scripts/explain_repeats_v3.py` for new explanations; the v0.1/v0.2 entry points
remain available to reproduce historical results.
