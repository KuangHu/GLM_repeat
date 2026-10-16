# Round 6 summary: threshold correction, VIPR v2, independent validation

All GPU work: model tattabio/gLM2_650M @08754cba, bf16, A40, attention policy
`always`. No Prodigal, no classifier, no genome-wide scan.

## 1. Threshold correction

The original `window_local` scorer and `results/scorer/report_exact/` are preserved
unchanged; a fresh `window_local` re-run reproduces its window scores byte-for-byte.
`background_calibrated` was added alongside it, with thresholds set from background
candidate DENSITY (the spec's 1/10/100 per Mb operating points) on 217 NEW
background windows (651 kb) disjoint from seed19 and from the validation panel.
No positive was read during calibration; thresholds were frozen (with job IDs) in
`results/calibration/thresholds/FROZEN.md` before evaluation.

| mode | 1/Mb | 10/Mb | 100/Mb |
|---|---|---|---|
| RP (0.651 Mb bg) | not established | tau 1.33 (9.2/Mb) | tau 0.89 (93.7/Mb) |
| exact (0.126 Mb bg) | not established | not established | tau 0.92 (87.3/Mb) |

Operating points are reported NOT ESTABLISHED rather than assigned a value where the
background volume cannot resolve them: one candidate already equals 1.5/Mb (RP) or
7.9/Mb (exact). The search grid is capped below ln(4) so calibration can never
return an unreachable threshold.

### w0011 reported separately

Its `window_local` thresholds were 1.596 (exact) and 1.454/1.487 (RP), all above the
ln(4) = 1.386 ceiling, so its zero was produced BY CONSTRUCTION.

| policy | exact score (rank) | RP score (rank) |
|---|---|---|
| window_local | 0.000 (14) | 0.000 (14) |
| bg100_per_Mb | 4.382 (3) | 4.106 (3) |

**Last round's reading was wrong.** The peak-cluster score does not reject the tRNA
cluster: with a reachable threshold it ranks #3 of 19, above three of the five
positives. What actually separates it is the copy-dependence diagnostic, not the
discovery score. At the stricter 10/Mb point w0011 is zero, but so are all 14
backgrounds and also VIPR and SA1 ART, so that zero is stringency, not
discrimination.

## 2. VIPR copy-dependence v2

94-nt vrRNA arm NOT run: no confidently homologous vrRNA exists (one perfect
self-match; best non-self 63% identity, inside the unconfirmed off-label candidate).
Recorded `missing_reference`. The off-label candidate stays unconfirmed.

60-nt YNNGG region, confound fixed (nearest perturbed unit 15-20 nt vs 10 nt in v1):

| definition | copy-dependence | seed SD |
|---|---:|---:|
| 12 x 5-nt units | -0.041 | 0.041 |
| 6 x 10-nt pairs | +0.006 | 0.011 |
| (CRISPR / ART, same protocol) | -1.18 to -1.36 | 0.012-0.019 |

VIPR shows **no array-like copy dependence**. Not a budget artefact: ART perturbed
only 56 nt and still gave -1.17. v1's apparent -0.413 is best explained by the
immediate-context confound that v2 removes.

## 3. Independent validation

New panel, scored exact and RP, frozen parameters, seed19 untouched as development
data. 13 positives + 13 GC-matched backgrounds.

- CRISPR: 6 loci *Pyrococcus furiosus* (NC_003413.1, **Archaea**) + 2 *Thermus
  thermophilus* (NC_006461.1), from GenBank `repeat_region /rpt_family=CRISPR`.
- ART: 5 manuscript-listed Staphylococcus phages. Extent defined by the 14-nt core
  rule; RT-start boundary NOT applied (records annotate only hypothetical proteins).
  Isolated singleton core matches excluded. Same lineage as SA1/MarsHill, grouped.
- VIPR: `missing_reference`, none added.

Results at the frozen exact threshold (tau 0.92):

- **13/13 labels recovered.** Label coverage 0.91-1.00; candidate precision
  (fraction of candidate inside the label) 0.70-1.00.
- **Ranking: all 13 positives occupy ranks 1-13; all 13 backgrounds 14-26.**
  P(positive > background) = 1.00. One background window (v0023, LY01) carries a
  candidate; it is unlabeled real DNA, not a confirmed false positive.
- RP vs exact on the same scorer: 14/14 candidates recovered, Spearman 0.989,
  13/13 positives non-zero. At the stricter RP 10/Mb point: 10/14 candidates,
  10/13 positives.

### Background candidate density does not transfer

| background set | nt | exact tau 0.92 | RP tau 0.89 | RP tau 1.33 |
|---|---:|---:|---:|---:|
| calibration (threshold set here) | 651k / 126k | 87.3/Mb | 93.7/Mb | 9.2/Mb |
| seed19 backgrounds | 42k | 190.5/Mb | 190.5/Mb | 0/Mb |
| validation backgrounds | 39k | 25.6/Mb | 25.6/Mb | 0/Mb |

A threshold calibrated on one background set gives ~2x the density on seed19
backgrounds and ~3.5x less on validation backgrounds. Counts are 1-11 candidates, so
intervals are very wide, but the operating point clearly does not transfer across
background composition. Any future density target must be calibrated on background
matched to the material actually being scanned. No FDR is claimed.
