# High-scoring background audit — NOT benchmark ground truth

w0011 (NC_028808.2:[86620,89620), SUSP1, labelled background, paired with VIPR w0003)
ranks #1-#3 of 19 on every pre-declared exact-score feature. It remains labelled
"background" in all metrics; this note only explains the score.

- GenBank (provenance/SUSP1.gb.gz) annotates 12 tRNAs inside the window
  (87462-89388; Leu, Lys, Ala, Gly, Thr, Val, Leu, Arg, Gln, Leu, Gln, Phe):
  a phage tRNA cluster.
- tRNA bases: 870/3000 (29%). Exact track (orientation-averaged): mean 0.907
  inside tRNAs vs 0.298 outside. 77% of bases with 16-nt smoothed score > 1.0
  lie inside annotated tRNAs.
- No exact direct repeat >= 20 nt in the window; the signal reflects conserved
  structured-RNA genes (e.g. GGTTCGAATCC T-arm motif), not copy-level repetition.

Implication: masked predictability also ranks tRNA arrays highly (as anticipated
in WORKFLOW_SETUP.md section 7: "predictable coding or structured single-copy DNA
can also score highly"). Whether tRNA clusters count as "repetitive organization"
is a scope decision for the project owner; until decided they stay background.
