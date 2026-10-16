# GLM_repeat — gLM2 masked-predictability repeat detection

Detecting repetitive sequence organization in DNA from a genomic language model
(`tattabio/gLM2_650M`), without motifs, protein anchors or family labels.

**Status: detector pilot closed and frozen; explainer v0.2 frozen for review.**
This repository holds code, frozen configuration and reports. Model weights, raw
genome downloads and per-base score tracks are not included.

## Start here

| | |
|---|---|
| `review_packet_explainer_v0.2/review/FINDINGS.md` | independent review of the explainer: 5 confirmed bugs, 6 quantified limitations |
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
- `unclassified` candidates carry no overlapping annotation in the records used.
  That is a statement about available annotation, not a claim of novelty.
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

Run the tests with `python -m pytest -q tests/`.
