# Candidate geometry and merging audit

CPU only: recomputed from existing score tracks. **No new model inference was run
for this audit.** Thresholds were not retuned; both frozen operating points
(tau 0.80 primary, tau 1.32 stringent) are exactly as frozen before the held-out
scan.

## 1. Geometry at both frozen operating points

| quantity | calibration 0.80 | calibration 1.32 | held-out 0.80 | held-out 1.32 |
|---|---:|---:|---:|---:|
| merged candidates | 457 | 42 | 467 | 30 |
| unique scanned Mb | 4.593 | 4.593 | 5.081 | 5.081 |
| candidates per unique Mb | 99.51 | 9.15 | 91.91 | 5.90 |
| union candidate span, % of scanned | 2.19% | 0.26% | 3.25% | 0.43% |
| span median | 162 | 263 | 194 | 398 |
| span p95 | 1,068 | 1,057 | 1,230 | 2,555 |
| span max | 2,626 | 1,669 | 3,458 | 3,445 |
| pre-merge clusters per candidate, median / max | 3 / 11 | 3 / 8 | 3 / 9 | 3.5 / 8 |
| source windows per candidate, median / max | 2 / 4 | 2 / 3 | 2 / 4 | 2 / 4 |
| peak-support fraction, median | 0.43 | 0.25 | 0.415 | 0.244 |
| max internal peak gap | 499 | 462 | 499 | 448 |

Candidates draw on a median of 3 pre-merge clusters from 2 source windows, which is
what the 50% stride plus two strands should produce for one underlying locus.
Maximum span slightly exceeds one 3,000-nt window only on held-out (3,458 nt from at
most 4 source windows), so merging extends a locus across adjacent windows but does
not chain along the contig. Internal peak gaps never exceed MAX_GAP_NT = 500, as the
clustering rule requires.

Peak support is a minority of each candidate's span (median 0.42 primary, 0.24
stringent): candidates are mostly *gap* between peaks. That is expected for repeat
architecture and is the reason span overlap alone is not treated as localization
evidence (section 4).

## 2. The low-tau calibration example: 6,236 clusters -> 50 candidates

Resolved. The collapse is **transitive interval merging**, not MIN_PEAKS rejection.

At tau = 0.30 over 6,106 calibration tracks:

| | value |
|---|---:|
| tracks rejected by MIN_PEAKS (exactly 1 peak) | **0** |
| tracks with 0 peaks | 1 |
| per-window clusters formed | 6,236 |
| merged candidates | 50 |
| pre-merge clusters per candidate, median / max | 80.5 / 704 |
| span median / max | 58,316 / 520,469 nt |
| union candidate span | 99.4% of scanned |

Nearly every window produces a cluster; each overlaps its neighbours, and the chain
propagates transitively until a handful of intervals span whole contigs. The low
yield per Mb is a small *count* of enormous intervals, not sparse signal.

This corrects the mechanism stated in the original freeze note, which attributed the
low branch to MIN_PEAKS. The correction is recorded in
results/pilot/calibration/thresholds/FROZEN.md. **The original selection, its sweep
and the frozen thresholds are preserved unchanged**; only the explanation was wrong,
and the conclusion (never operate on the low branch) is unaffected. Candidate
generation was not modified, so no calibration claim is restated or versioned.

## 3. Cross-tier correspondence (nesting NOT assumed)

Counted from each side separately, since overlaps are many-to-one.

| set | both tiers (primary side) | primary_only | both tiers (stringent side) | stringent_only |
|---|---:|---:|---:|---:|
| calibration | 38 | 419 | 39 | **3** |
| held-out | 30 | 437 | 30 | 0 |

The tiers are **not** nested in general: calibration has 3 stringent-only candidates
with no overlapping primary candidate. Raising tau can split one broad peak into two
resolvable peaks, forming a cluster that did not exist at lower tau. Held-out happens
to be nested (0 stringent-only), which is an observation about this material, not a
property of the method.

## 4. Label overlap judged together with peak support

Reported in results/pilot/heldout/report_primary/report.md. Broad span overlap alone
is not accepted as precise localization, so each locus is reported with the peak
support, peak-block count and internal peak gaps of its best-overlapping candidate.

| family | loci | locus coverage | localization precision | peak support | peak blocks | internal gaps |
|---|---:|---:|---:|---:|---|---|
| CRISPR | 8/8 | 1.00 | 1.00 | 0.43 | 15-51 | 44-60 nt |
| ART | 5/5 | 1.00 | 0.70 | 0.18 | 7-9 | ~230 nt |

The evidence for precise localization is not the overlap figure but the peak
geometry: CRISPR candidates resolve 15-51 peak blocks separated by 44-60 nt, and ART
candidates 7-9 blocks separated by ~230 nt, matching the known repeat-spacer and
core-array periodicities respectively. One short locus
(Thermus_thermophilus_HB8_CRISPR_874397, 338 nt) is the weakest case: precision 0.49,
peak support 0.25, max internal gap 349 nt.

VIPR's labelled 60-nt YNNGG region is development material and is reported
separately; it is never pooled with these held-out families, and the off-label
candidate NC_028808.2:[27682,27906) is not counted as a VIPR recovery anywhere.

---

## 5. CORRECTION 2026-10-15: peak spacing terminology and coordinate reconciliation

Section 4 above compared peak gaps to annotated unit spacings. That comparison was
**invalid**, for two reasons, and the numbers quoted there should not be used.

1. **Edge gap is not start-to-start spacing.** An annotated unit spacing is
   start-to-start (pitch). An edge gap is end-of-peak to start-of-next-peak, which
   is smaller than the pitch by roughly the peak width. They are different
   quantities and were compared directly.
2. **The figures quoted were maxima, not typical values.** "~230 nt" for ART and
   "44-60 nt" for CRISPR were *maximum* internal edge gaps from the geometry audit,
   not medians.

Coordinate-level reconciliation (results/pilot/heldout/reconciliation/, same frozen
tracks and tau = 0.80, no new inference). Annotated units: 14-nt ART cores from the
manuscript rule; CRISPR repeat copies derived de novo as the most frequent periodic
k-mer within each annotated array.

| locus | units | annotated pitch | peak pitch | peak edge gap | matched | merged | extra | missing |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P. furiosus CRISPR 27091 | 48 | 67 | 68 | 39 | 48 | 0 | 1 | 0 |
| P. furiosus CRISPR 260714 | 20 | 67 | 68 | 39 | 20 | 0 | 0 | 0 |
| P. furiosus CRISPR 312405 | 21 | 68 | 68 | 39 | 21 | 0 | 2 | 0 |
| P. furiosus CRISPR 623119 | 31 | 68 | 68 | 40 | 28 | 1 | 0 | 0 |
| P. furiosus CRISPR 695937 | 44 | 67 | 67 | 39 | 44 | 0 | 2 | 0 |
| P. furiosus CRISPR 1064076 | 22 | 67 | 67 | 38 | 22 | 0 | 0 | 0 |
| T. thermophilus CRISPR 872101 | 14 | 76 | 75 | 41 | 14 | 0 | 1 | 0 |
| T. thermophilus CRISPR 874397 | 5 | 76 | 76 | 41 | 5 | 0 | 0 | 0 |
| Madawaska ART | 5 | 204 | 192 | 154 | 5 | 0 | 0 | 0 |
| LY01 ART | 5 | 203 | 195 | 162 | 5 | 0 | 1 | 0 |
| S6 ART | 5 | 203 | 171 | 152 | 5 | 0 | 1 | 0 |
| PALS_2 ART | 5 | 203 | 173 | 152 | 5 | 0 | 1 | 0 |
| UFV_DC4 ART | 5 | 203 | 200 | 159 | 5 | 0 | 0 | 0 |

Reconciled statements, on matched quantities only:

- **CRISPR: peak pitch equals annotated repeat pitch.** 67-68 nt against 67-68 nt in
  Pyrococcus, 75-76 against 76 in Thermus. Every annotated repeat copy is covered:
  **0 missing across all 230 units in all 8 arrays.** One peak in
  CRISPR_623119 spans more than one repeat copy (merged); 0-2 extra peaks per array
  fall inside the array span but on no annotated copy, most plausibly spacer
  sequence that is itself predictable.
- **ART: peak pitch is close to, but consistently below, the annotated core pitch.**
  171-200 nt against 203-204 nt. 5/5 cores matched in every phage, none missing,
  0-1 extra peaks. The shortfall is real and unexplained by this audit; peaks are
  wider than the 14-nt cores and their starts need not align with core starts.
- Edge gaps (CRISPR 38-41 nt, ART 152-162 nt) are inter-peak spacers and are
  **not** comparable to annotated pitch. They are reported here only as geometry.

The earlier phrasing "~230-nt gaps matching 170-215-nt spacing" is withdrawn.

---

## 6. CORRECTION 2026-10-15: ART pitch on matched coordinate sets

Section 5 compared ART peak pitch (171-200 nt) with annotated core pitch
(203-204 nt) and called the shortfall real and unexplained. That comparison used
**different coordinate sets**: the peak pitch was computed over *all* peaks inside
the locus, including `extra` peaks that overlap no annotated core, which shortens
the apparent spacing.

Recomputed over matched sets only — peaks classified `matched`, and the cores they
actually matched, in the same order:

| locus | core pitch | all-peak pitch | matched-peak pitch | pitch of the matched cores |
|---|---:|---:|---:|---:|
| Madawaska | 204 | 192 | 192 | 204 |
| LY01 | 203 | 195 | 202 | 203 |
| S6 | 203 | 171 | 191 | 203 |
| PALS_2 | 203 | 173 | 192 | 203 |
| UFV_DC4 | 203 | 200 | 200 | 203 |

Most of the apparent discrepancy was the coordinate-set mismatch: S6 and PALS_2
move from 171/173 to 191/192 once extra peaks are excluded. The residual shortfall
is 1-12 nt, not ~32 nt. A plausible account of the residual is that peaks are wider
than the 14-nt cores and a peak's start need not align with its core's start, so
peak starts drift slightly relative to core starts. That remains an account, not a
demonstration.

## 7. Terminology: repeat-pattern clusters

The 185 groups produced by deduplicating explainer models are **repeat-pattern
clusters** — groups of candidates sharing a similar consensus sequence. They were
previously called "repeat families", which implies biological relatedness that
sequence similarity alone does not establish. They are not biological systems and
no claim of shared origin, function or classification is made.
