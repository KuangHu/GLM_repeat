# Held-out pilot scan

Configuration frozen before scanning: results/pilot/PILOT_CONFIG_FROZEN.md. Mode random_partition (groups 20, repeats 2, seed 20261014), policy always, model 08754cba59.
Threshold 0.8 (primary_100perMb), calibrated on calibration material only and frozen before this scan.

## Scan summary

- windows scored: 3381
- unique scanned genomic span: 5.0813 Mb
- merged candidates: 467 (from 1290 per-window clusters)
- **candidate yield: 91.91 per unique scanned Mb** (a yield, not an FDR)

### Yield per genome (densities are not forced to match)

| contig | unique Mb | candidates/Mb |
|---|---:|---:|
| MZ779063.1 | 0.263 | 205.18 |
| LC680885.1 | 0.267 | 187.23 |
| MN091626.1 | 0.269 | 182.33 |
| OR836606.1 | 0.259 | 181.54 |
| MW349129.1 | 0.265 | 180.83 |
| NC_006461.1 | 1.850 | 85.96 |
| NC_003413.1 | 1.908 | 31.44 |

## Held-out loci (grouped by genome / family)

Overlap alone is not evidence of precise localization, so peak support inside the best-overlapping candidate is reported alongside it.

| locus | family | family group | locus nt | candidates overlapping | locus coverage | best localization precision | peak support of that candidate | its peak blocks | max internal peak gap |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Pyrococcus_furiosus_DSM3638_CRISPR_27091 | CRISPR | CRISPR_Pyrococcus_genome | 3310 | 1 | 1.00 | 0.96 | 0.41 | 51 | 57 |
| Pyrococcus_furiosus_DSM3638_CRISPR_260714 | CRISPR | CRISPR_Pyrococcus_genome | 1332 | 1 | 1.00 | 0.96 | 0.42 | 21 | 51 |
| Pyrococcus_furiosus_DSM3638_CRISPR_312405 | CRISPR | CRISPR_Pyrococcus_genome | 1527 | 1 | 1.00 | 0.99 | 0.44 | 23 | 47 |
| Pyrococcus_furiosus_DSM3638_CRISPR_623119 | CRISPR | CRISPR_Pyrococcus_genome | 2058 | 1 | 1.00 | 1.00 | 0.46 | 29 | 44 |
| Pyrococcus_furiosus_DSM3638_CRISPR_695937 | CRISPR | CRISPR_Pyrococcus_genome | 3056 | 1 | 0.99 | 1.00 | 0.42 | 46 | 44 |
| Pyrococcus_furiosus_DSM3638_CRISPR_1064076 | CRISPR | CRISPR_Pyrococcus_genome | 1468 | 1 | 1.00 | 1.00 | 0.43 | 22 | 60 |
| Thermus_thermophilus_HB8_CRISPR_872101 | CRISPR | CRISPR_Thermus_genome | 1100 | 1 | 0.99 | 1.00 | 0.47 | 15 | 45 |
| Thermus_thermophilus_HB8_CRISPR_874397 | CRISPR | CRISPR_Thermus_genome | 338 | 1 | 1.00 | 0.49 | 0.25 | 6 | 349 |
| Madawaska_ART | ART | ART_Staphylococcus_lineage | 866 | 1 | 1.00 | 0.70 | 0.19 | 7 | 228 |
| LY01_ART | ART | ART_Staphylococcus_lineage | 871 | 1 | 1.00 | 0.70 | 0.16 | 9 | 229 |
| S6_ART | ART | ART_Staphylococcus_lineage | 864 | 1 | 1.00 | 0.70 | 0.18 | 9 | 232 |
| PALS_2_ART | ART | ART_Staphylococcus_lineage | 865 | 1 | 1.00 | 0.70 | 0.19 | 9 | 231 |
| UFV_DC4_ART | ART | ART_Staphylococcus_lineage | 866 | 1 | 1.00 | 0.70 | 0.15 | 8 | 231 |

### Outcome by family (reported separately, never pooled)

| family | family groups | loci | loci with >=1 overlapping candidate | median locus coverage | median localization precision | median peak support |
|---|---:|---:|---:|---:|---:|---:|
| CRISPR | 2 | 8 | 8 | 1.00 | 1.00 | 0.43 |
| ART | 1 | 5 | 5 | 1.00 | 0.70 | 0.18 |
| VIPR | 0 | 0 | not tested | - | - |

Loci within one family group are not independent: the Pyrococcus CRISPR arrays share a genome, and all ART loci share one Staphylococcus lineage. Counts of loci overstate the number of independent observations.

## Annotation of candidates (applied after scoring, never a filter)

| category | candidates | median score | median span nt |
|---|---:|---:|---:|
| unclassified | 393 | 1.582 | 205 |
| tRNA_single | 42 | 1.751 | 86 |
| tRNA_cluster | 16 | 2.305 | 250 |
| CRISPR_array | 8 | 5.673 | 1541 |
| ART_core_array_motif | 5 | 3.365 | 1236 |
| other_annotated_ncRNA | 3 | 2.967 | 399 |

`unclassified` candidates are retained and reported (393 of 467); they are not discarded and are not implied to be artefacts. Highest-scoring unclassified candidates:

- OR836606.1:151262-152315 score 4.428734621860781, 21 peaks, span 1053 nt
- NC_003413.1:1091159-1091856 score 4.405313256280017, 11 peaks, span 697 nt
- MW349129.1:141073-142114 score 4.224919012310553, 15 peaks, span 1041 nt
- OR836606.1:137600-138800 score 4.017384766435808, 14 peaks, span 1200 nt
- NC_006461.1:1609897-1611762 score 3.7402103464392393, 16 peaks, span 1865 nt
- NC_006461.1:1229552-1230975 score 3.6838867180825523, 15 peaks, span 1423 nt
- OR836606.1:40441-41091 score 3.622532223493779, 12 peaks, span 650 nt
- MZ779063.1:185368-186871 score 3.5542745928878876, 9 peaks, span 1503 nt
- LC680885.1:96729-97568 score 3.52504788798321, 12 peaks, span 839 nt
- OR836606.1:232905-233559 score 3.500803859072405, 8 peaks, span 654 nt

## VIPR, reported separately

No VIPR locus is present in the held-out material: no independent natural VIPR locus could be obtained (recorded `missing_reference`). **VIPR generalization is unresolved by this pilot and is not claimed.**

The original VIPR 60-nt YNNGG region is development data (seed19 w0003, NC_028808.2:[29169,29229)) and is reported apart from any off-label hit. The off-label candidate NC_028808.2:[27682,27906) remains a sequence-supported repeat candidate with vrRNA identity unconfirmed; it is not counted as a VIPR recovery anywhere in this report.

See results/diagnostics/vipr_v2/report/report.md for the copy-dependence outcome on that region.

## Limits

- Candidate yields are yields. No FDR, no calibrated biological probability.
- Calibration genomes are ordinary DNA and may contain real repeats; they are not verified negatives.
- Held-out loci are database or manuscript-derived annotations, not experimentally mapped RNA boundaries.
- Copy dependence was not used for candidate eligibility.
