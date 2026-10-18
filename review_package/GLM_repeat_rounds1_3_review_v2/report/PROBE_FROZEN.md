# No-model repeat probe: frozen parameters and validated operating characteristics

Diagnostic only. Does not modify the detector, the explainer or any frozen output,
and is never fed back into candidate selection.

## Frozen parameters

| parameter | value |
|---|---|
| seed | exact k-mer, k = 10 |
| seed filter | >= 3 distinct bases; occurrences paired only if >= k apart |
| self-diagonal | excluded by construction (only distinct occurrences are paired) |
| extension | Biopython `PairwiseAligner`, mode = local |
| match / mismatch | +2 / -3 |
| gap open / extend | -5 / -2 (affine) |
| extension flank | 40 nt either side of each seed |
| max seed pairs per region | 40 |
| minimum alignment | 20 nt, defined as alignment score / match_score |
| identity threshold | none; extension score is the sole criterion |
| probed region | candidate span + 500 nt flank each side |

## Sensitivity (synthetic panel, planted repeats, 4 copies, 20 replicates per cell)

Background is shuffled real genomic DNA, so composition is realistic and the
planted repeat is the only signal.

By unit length:

| unit nt | 8 | 12 | 20 | 40 | 80 |
|---|---:|---:|---:|---:|---:|
| sensitivity | 0.20 | 0.22 | 0.33 | 0.77 | 0.87 |

By per-base divergence:

| divergence | 0% | 5% | 10% | 20% |
|---|---:|---:|---:|---:|
| sensitivity | 0.58 | 0.56 | 0.50 | 0.26 |

By indel type:

| indel | none | 1 nt | 3 nt | mixed |
|---|---:|---:|---:|---:|
| sensitivity | 0.52 | 0.47 | 0.43 | 0.48 |

**Unit length, not indel content, is the limiting factor.** At 12 nt units - the
length that dominates the real panels - sensitivity is 0.22.

The earlier 90% / 85% figures were measured on real modelled candidates with >= 5
copies, an easier stratum than this panel. Both are valid; they measure different
things, and the synthetic table is the one to use when interpreting a negative.

## False-positive rate (shuffled real DNA, no planted repeat, n = 120 per length)

| region nt | 600 | 1,200 | 2,000 |
|---|---:|---:|---:|
| positive calls | 4.2% | 13.3% | 25.8% |

A 10-mer recurs by chance in a kilobase-scale window, so positives must always be
compared against this baseline at the matching region length.

## Interpretation rules

1. A negative means **no qualifying local similarity was detected under these
   parameters**. It does NOT mean the region contains no repeat. A negative is only
   as informative as the matching cell of the sensitivity table; at 8-20 nt units
   that is 20-33%.
2. A positive must be compared against the false-positive baseline for the same
   region length before being called evidence.
3. The count of `gapped_similarity` calls is a floor, not an estimate of the indel
   coverage gap. With ~45% sensitivity for indel-containing repeats overall and
   ~20% at short unit lengths, the detected count understates the true number by an
   unknown factor. **3 detected calls do not show that the indel gap is small.**
4. The shuffled-DNA null destroys native composition structure, so the true chance
   rate in real sequence may be higher than measured here, which would further
   reduce any apparent excess.
