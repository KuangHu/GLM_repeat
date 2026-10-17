# Frozen RP operating points for the held-out pilot scan

Derived from CALIBRATION material only: *M. jannaschii* NC_000909.1,
*D. radiodurans* NC_001263.1, *Staphylococcus* phage K NC_005880.2 and phage Twort
NC_007021.1 — 3,053 windows, 4.593 unique scanned Mb. No held-out genome and no
known locus was read. Frozen 2026-10-15, while held-out scoring (job 26976966) was
still running and before any held-out candidate was examined.

| operating point | tau | achieved yield | merged candidates |
|---|---:|---:|---:|
| primary discovery, ~100/Mb | **0.80** | 99.51 /Mb | 457 |
| stricter tier, ~10/Mb | **1.32** | 9.15 /Mb | 42 |

Both are below the measurement ceiling ln(4) = 1.3863, so neither is unreachable.

## The yield/threshold curve is not monotonic

Yield rises with tau to a maximum of 428.7 /Mb at tau = 0.53, then falls. Below
that maximum the scorer degenerates: almost every position exceeds the threshold,
the window collapses into one continuous peak, and MIN_PEAKS = 2 is not satisfied,
so yield is low for a reason that has nothing to do with repeat content. A naive
"smallest tau meeting the target" search returns tau = 0.30 from that degenerate
branch.

Operating points are therefore taken from the **stringent branch only**, at or
above the density maximum. This was found and fixed during calibration, before the
freeze; the first (incorrect) sweep was discarded and recomputed.

Yields are candidate yields per unique scanned genomic Mb. They are not FDRs.
Calibration genomes are ordinary DNA and may contain real repeats, including CRISPR
arrays; they are calibration material, not verified negatives. Per-genome yields are
reported as they fall and are not forced to match.
f6bc48cf522cfa1f185bb6fe751d66ad08304308f05c97680a533ef54a9e14e8  thresholds.json

---

## CORRECTION 2026-10-15 (explanation only; thresholds and selection UNCHANGED)

The mechanism stated above for the degenerate low-tau branch was **wrong**. A
geometry audit of the same tracks (results/pilot/calibration/audit_lowtau/, no new
inference) shows that at tau = 0.30:

- windows rejected by MIN_PEAKS (exactly one peak): **0**
- tracks with zero peaks: 1 of 6106

So MIN_PEAKS rejection plays no part. The real mechanism is **transitive interval
merging**: 6,236 per-window clusters chain into 50 candidates, a median of 80.5
(max 704) pre-merge clusters per candidate, median span 58,316 nt, maximum 520,469
nt, covering 99.4% of the scanned material. Low yield per Mb at low tau reflects a
small *count* of genome-spanning intervals, not scarce signal.

The conclusion is unchanged and if anything stronger: the low branch must not be
used as an operating point. **The frozen thresholds (0.80 and 1.32), the
stringent-branch selection rule, and the sweep they were derived from are retained
exactly as recorded above.** Nothing in candidate generation was altered, so no
calibration claim is restated or versioned.

Geometry at the frozen operating points is well behaved and shows no runaway
merging:

| point | tau | candidates | median span | p95 span | max span |
|---|---:|---:|---:|---:|---:|
| primary | 0.80 | 457 | 162 nt | 1,068 nt | 2,626 nt |
| stringent | 1.32 | 42 | 263 nt | 1,057 nt | 1,669 nt |

Maximum span at both points is smaller than one 3,000-nt window, so no candidate
chains across windows at the operating points.
