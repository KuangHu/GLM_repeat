# Label-blind peak-cluster scorer — frozen design

**Authorship: this design is mine (the agent), not a user specification.** The
user's spec file was unavailable; they authorised me to pre-declare one. Declared
2026-10-15 BEFORE running the scorer on any track, exact or RP.

## Rationale

Block mean and block top-20% (the existing baselines, retained unchanged) reward a
window for being *predictable somewhere*. They cannot distinguish one long
predictable stretch (a conserved gene, a structured RNA) from several short
predictable stretches separated by unpredictable spacer (a repeat array). The
peak-cluster score is intended to reward the second pattern specifically, using
only the score track.

## Inputs

Per-base `above_uniform` tracks (nats), both scoring orientations, as produced by
`reference_score.py` / `timing_score.py` in either `exact` or `random_partition`
mode. No labels, annotations, family names, motif sequences, repeat finders, CDS
calls or low-complexity filters are inputs at any stage.

## Frozen parameters

| name | value | basis |
|---|---|---|
| `SMOOTH_NT` | 16 | the manuscript's smoothing width for this measurement |
| `THRESH_K` | 3.0 | peak = median + 3 robust SD, a conventional outlier cut |
| `MIN_PEAK_NT` | 8 | shorter than the shortest repeat unit class in scope (5 nt) would admit single-base noise; 8 keeps 2+ smoothed positions |
| `MAX_GAP_NT` | 500 | known array spacings in scope are < 500 nt between unit starts |
| `MIN_PEAKS` | 2 | a cluster must contain repetition; one peak is not a cluster |

`MIN_PEAK_NT` and `MAX_GAP_NT` are informed by general repeat architecture (unit
lengths ~5-100 nt, inter-unit spacing < 500 nt), **not** by any gLM2 score observed
on the seed panel. No parameter was chosen by inspecting seed19 scores, and none
will be re-tuned in this round.

## Algorithm

Per window and per orientation, independently:

1. Smooth the track with a centred running mean of width `SMOOTH_NT`.
2. Window-local robust baseline: `med = median(smoothed)`,
   `sd = 1.4826 * median(|smoothed - med|)`. Threshold `tau = med + THRESH_K * sd`.
   The baseline is per window, so no global threshold is fitted on any panel.
3. A **peak** is a maximal run of >= `MIN_PEAK_NT` consecutive positions with
   `smoothed > tau`. Record its start, end, mean and summit.
4. Peaks whose gap (end of one to start of next) is <= `MAX_GAP_NT` form a
   **cluster**. A cluster is kept only if it holds >= `MIN_PEAKS` peaks.
5. **Cluster score** (frozen):
   `score = median(peak_mean) * log2(1 + n_peaks)`
   Median peak strength resists one dominant peak; the log term credits repetition
   with diminishing returns.
6. **Window score** = max cluster score over clusters and over the two orientations
   (0 if no cluster survives), matching the existing max-over-orientations
   convention. Both orientation tracks are retained.

## Overlapping windows

Peaks and clusters are emitted in **genomic** coordinates
(`window.start0 + offset`). Before any counting, clusters from different windows on
the same contig that overlap are merged into one candidate, keeping the maximum
score and the union span. Overlapping windows therefore never contribute
independent peaks or independent candidates. Window ranking is reported separately
from merged candidates and is explicitly a per-window quantity.

## Reporting (kept separate, per the instruction)

1. candidate coordinates (merged, genomic, BED-style 0-based half-open)
2. peak support per candidate (peak count, lengths, spacings, median peak mean)
3. label localization: overlap of candidates with labels — computed **after**
   scoring, never used to form or score candidates
4. window ranking: peak-cluster score alongside the unchanged `blk*_mean_max` and
   `blk*_top20_max` baselines

## Declared limits

seed19 is development data. No FDR, no calibrated threshold, no recall-loss claim
is produced from it. A window score of 0 means no cluster passed the frozen
criteria in that window; it is not evidence of absence of repeats.

## Post-run addendum (2026-10-15, after first execution)

Observed defect, recorded without changing any parameter this round: the threshold
`median + THRESH_K * robust SD` can exceed the measurement's ceiling ln(4) = 1.386
in windows that are broadly predictable and highly dispersed. In that case the
window yields zero peaks by construction. On seed19 this affects w0011 only
(tau = 1.596, both orientations). Any future use must clip tau below ln(4) or
replace the MAD-based threshold; parameters above remain as frozen for this round.
