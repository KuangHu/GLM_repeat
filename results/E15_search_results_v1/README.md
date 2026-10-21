# E15 route-B search: blastp + jackhmmer over 232 Mb of our own panels

Run 20 October 2026. Searched collection: the 122,704-protein predicted proteome of
all eight scanned panels (`protein_bundle_v1`), NOT a public database.

## Queries, validated against e15.md

| protein | contig | coords | strand | len |
|---|---|---|---|---|
| `r2_phage\|cds00003113` | NC_048046.1 | 114589–115999 | − | 469 aa |
| `r3_phage\|cds00020343` | NC_019406.1 | 107046–108456 | − | 469 aa |

Both coordinates, strands and lengths match e15.md exactly. All seven E15 cores and
the 987-nt array length reproduce from our genome copy base-for-base.

## Result 1 — the protein family is much larger than two

| method | setting | distinct subjects |
|---|---|---:|
| blastp | E ≤ 1e-5 | 7 (2 known + 5 new) |
| jackhmmer | 5 iterations, E ≤ 1e-3 | **57** |

blastp's new hits: 36–38% identity (NC_070966.1, NC_070967.1, NC_049977.1) and
27–29% (NC_021328.1, NC_073474.1). The two known proteins recover each other at
94.9%, which is the positive control.

NC_070966.1 and NC_070967.1 are 74.9% identical over 211 kb — sister phages, so
they count as one genome cluster, not two independent discoveries.

## Result 2 — no new E15-like ARRAY was found

The array is the defining feature; protein similarity alone is insufficient (e15.md §1).

Positive control, blastn of the 987-nt E15 array:

| subject | identity | aligned | coords |
|---|---:|---:|---|
| NC_048046.1 (self) | 100.0% | 987 | 116154–117140 |
| NC_019406.1 (peer) | 95.3% | 992 | 108612–109601 |

The peer is recovered over its full length at exactly the documented envelope.

Against the five new blastp regions (±500 nt): **0 HSPs at E ≤ 1.**

Of the 57 jackhmmer proteins, 54 have some candidate within 12 kb and 42 have a
modelled array, but only 12 have a 12-nt array, and the highest-scoring of those
are **tRNA gene clusters**, not E15 analogs:

| protein | contig | bits | dist to CDS | consensus | what it is |
|---|---|---:|---:|---|---|
| `r2_phage\|cds00003113` | NC_048046.1 | 660 | 43 | CTACGCCACCCT | **E15 (self)** |
| `r3_phage\|cds00020343` | NC_019406.1 | 653 | 42 | CTACGCCACCCT | **known peer** |
| `r3_phage\|cds00017789` | NC_048730.1 | 514 | 3053 | GGTTCGAGTCCC | tRNA cluster (40 tRNAs) |
| `r3_phage\|cds00002841` | NC_048724.1 | 504 | 3018 | GGTTCGAATCCC | tRNA cluster (42 tRNAs) |
| `db10_phage\|cds00030535` | NC_048728.1 | 494 | 704 | GGGTTCGAATCC | tRNA cluster (33 tRNAs) |
| `r3_phage\|cds00004091` | NC_049977.1 | 478 | 274 | TTGAAGTAACTC | non-tRNA array, 5 copies — best Tier-B lead |

Only the two known loci carry a `CTACGCCACC..` core at ~42 nt from the CDS.

## Tier assignment (e15.md §7)

- **Tier A, close sequence homolog:** 1, the known peer — recovered as the positive
  control, so **zero new** Tier-A loci.
- **Tier B, architectural analog:** no confident call. `NC_049977.1` is the strongest
  lead (478 bits, array 274 nt from the CDS, 5 × 12-nt copies, non-tRNA consensus),
  but it has no nucleotide homology to the E15 array and its array is on the
  opposite side of the CDS, which e15.md §4 says to flag as a different repeat slot
  rather than merge.
- **Tier C, exploratory:** the remaining 12-nt-array loci, mostly 3 copies.

## A caution worth recording

The two strongest new protein hits, NC_070966.1 (38%) and NC_070967.1 (37%), each
have a high-scoring detector candidate immediately 5' of the gene — 139133–140309
(score 3.73, 10 peaks) and 138306–139442 (score 3.53, 9 peaks) — in the same
relative position as the E15 array. But the explainer models neither, their 12-mer
content is at chance (p = 1.0), and they carry a tRNA-Leu overlap. Our own earlier
review had recorded the first as `no_repeat_found / not_worth_followup`.

These are either diverged arrays our ungapped explainer cannot represent, or they
are not arrays at all. **The data here does not decide between those**, and they
should be inspected manually rather than counted either way.

## Scope limits

- Searched 232 Mb of our own panels only. A public-database search was not run and
  would very likely add family members.
- jackhmmer at 5 iterations can drift; the 57 members are a retrieval set, not a
  validated family. No TROVE/Ro60 domain verification was performed, so "Ro60-like"
  here means "similar to the E15 query protein", not a domain-level assignment.
- Route A (full-array and inter-core interval search across all genomes) and route C
  (de novo array-first discovery) were not run.

## Result 3 — the negative is now genome-wide, not just the five blastp regions

The 987-nt E15 array was searched by blastn against the 10-kb context of **every one
of the 57 family proteins** (~1.2 Mb of sequence), at a deliberately permissive
E ≤ 10:

| subject | identity | aligned | evalue |
|---|---:|---:|---:|
| E15 self context | 100.0% | 987 | 0.0 |
| known peer context | 95.3% | 992 | 0.0 |

**Two HSPs, both expected. Nothing else, anywhere.** Within the 232 Mb searched, the
E15 array sequence occurs at exactly two loci, and every other member of the
Ro60-like protein family lacks it.

## Folder contents

| path | contents |
|---|---|
| `queries/E15_Ro60_like_proteins.faa` | the two 469-aa query proteins |
| `queries/E15_NC048046_array.fna` | the original 987-nt E15 array |
| `queries/E15_arrays_genomic_forward.fna` | both known array envelopes, genomic forward |
| `queries/E15_arrays_reverse_complement.fna` | reverse complements; **not** extra biological queries |
| `queries/E15_repeat_cores.fna` | all 14 cores with locus and coordinate provenance |
| `queries/E15_unique_core_seeds.fna` | the 5 distinct 12-nt seeds |
| `queries/E15_intercore_intervals.fna` | 12 intervals, both bounding cores excluded |
| `queries/E15_array_plus_3kb_flanks.fna` | each array with 3-kb flanks |
| `family/ro60_like_proteins_all57.faa` | all 57 jackhmmer family proteins, sorted by bitscore |
| `loci/best_candidate_regions.fna` | the highest-scoring detector candidate within 12 kb of each protein (54 loci) |
| `loci/protein_plus_10kb_context.fna` | each protein with 10 kb either side (57 loci) |
| `tables/ro60_loci.tsv` | master table: protein, coordinates, bitscore, nearest array, architecture, consensus, spacing |
| `tables/blastp_hits.tsv`, `tables/jackhmmer.tbl` | raw search output |
| `tables/blastn_E15array_vs_all_family_contexts.tsv` | the genome-wide array search |

Coordinates are 0-based half-open throughout. Protein IDs are `<panel>|<cds_id>`
and join to `protein_bundle_v1/protein_index.tsv`.

## Reproduction check against e15.md

Regenerated from our genome copies and verified against the document:

- all 7 E15 cores, base for base; array length 987 nt
- all 7 peer cores; array length 990 nt
- inter-core lengths — E15 163/147/118/185/158/132, peer 163/147/118/187/160/131
- 5 distinct 12-nt seeds across both loci
- both 469-aa proteins at the documented coordinates and strand

One correction made during this work: deriving the peer core positions from the
document's inter-core lengths put R5 at 109275, one base off. The true position is
109274, and the arithmetic version produced a frameshifted core (TACGCCACCCAC) and a
spurious sixth seed. Cores here were taken from the sequence, not from arithmetic.

`e15.md` section 2 refers to an `E15_search_queries/` folder that does not exist on
this system; the `queries/` directory here regenerates that file set from the
coordinates.
