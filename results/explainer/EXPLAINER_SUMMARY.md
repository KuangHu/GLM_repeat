# Repeat explainer: first results

Separate module. Detector outputs are read-only and unmodified; no gLM2 inference
is used. Discovery is from native sequence only — no motif, family name or
annotation is supplied — and searches the candidate span plus 500-nt flanks in
full, not the gLM2 peak intervals.

**Indel limitation (applies to everything below):** alignment is ungapped,
substitution-distance only. Arrays whose copies differ by insertions or deletions
are under-reported: copies may be split, truncated or missed, and reported spacer
lengths absorb the indel. Recorded per model as
`indel_handling: ungapped_substitutions_only`.

## Controls, discovered de novo

Known motifs were NOT supplied. 14 of 15 controls produced a model.

| control | unit | copies | pitch | consensus (discovered) |
|---|---:|---:|---:|---|
| P. furiosus CRISPR x6 | 30 | 20-50 | 67-68 | GTTCCAATAAGACTAAAATAGAATTGAAAG (and strand variants) |
| T. thermophilus CRISPR x2 | 36 | 5-15 | 76 | GTTGCAAGGGATTGAGCCCCGTAAGGGGATTGCGAC |
| Madawaska / LY01 ART | 30-33 | 3 | 191-195 | TTTATAGACATTATATGAATACGTATAATGATT |
| S6 / PALS_2 ART | 23 | 5 | 203 | CATTATATGAATACGTATAATGA |
| UFV_DC4 ART | 16 | 5 | 203 | TATGAATACGTATAAT |
| SUSP1 off-label (unconfirmed) | 39 | 2 | 143 | AGTCTCATTTTAATGTAGCTGAGACTAACTACCCATACC |

The recovered CRISPR pitches (67-68, 76 nt) equal the annotated pitches, and every
ART consensus contains the manuscript's 14-nt core `ATATGAATACGTAT` without it
having been provided.

Weakness: Madawaska and LY01 recover only 3 of 5 cores. Motif extension produced a
33-nt and 30-nt unit, long enough that the two most divergent copies exceed the
mismatch budget. Over-extension costs copies.

## VIPR, kept separate from detector recovery

The labelled 60-nt YNNGG region (NC_028808.2:[29169,29229)) is the **one control the
motif explainer could not model**: no motif with >=2 ungapped copies. Seed-based
discovery needs an exact 12-mer shared between copies, and a degenerate 5-nt YNNGG
pattern has none.

Its alternative explanation did fire, de novo: `periodic_self_match`, **period 5,
identity 0.673** on the candidate span. So the periodicity is detectable, but as a
profile, not as a motif model with per-copy coordinates.

**Explainer outcome for VIPR is therefore: periodic profile only, no copy-level
model.** This is reported separately from, and does not inherit, the detector's
recovery of that region. VIPR generalization remains unresolved.

The off-label candidate NC_028808.2:[27682,27906) is modelled (39-nt unit, 2 copies)
but stays an unconfirmed repeat candidate; it is never counted as a VIPR recovery.

## All 467 primary held-out candidates

| outcome | count |
|---|---:|
| explained | 337 |
| unresolved (no model) | 130 |

| architecture | count |
|---|---:|
| dispersed_repeat | 181 |
| regular_spaced_array | 100 |
| variable_spacer_array | 52 |
| tandem_repeat | 4 |

Nothing was forced: 28% of candidates have no sequence model, which is a reported
outcome. All 467 were processed, not only the unclassified ones.

Emitted per model: consensus, unit length, per-copy coordinates with mismatch
counts, copy alignments, PWM position probabilities (`pwm_*.npy`), start-to-start
spacing and edge gaps with their dispersion, peak-to-copy correspondence
(`copies.tsv: overlapping_peaks`), alternative explanations, and the unresolved
list with its alternatives.

## Review set

337 explained candidates deduplicate to **185 families** by strand-canonical
consensus similarity. The review set holds 11 interpretable `unclassified`
representatives spanning 4 architectures plus 8 high-scoring candidates with no
sequence explanation.

`unclassified` means only that the span carries no overlapping annotation in the
GenBank records used here. **It is not a claim of novelty, function, or a new
biological system.** One example worth noting factually and nothing more:
NC_003413.1:1091159-1091856 carries a 25-nt unit in 12 copies at 67-nt pitch — the
same pitch as the annotated CRISPR arrays in that genome, outside their annotated
spans.
