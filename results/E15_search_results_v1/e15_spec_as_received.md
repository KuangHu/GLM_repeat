# Finding E15-like repeat arrays with RNA-related protein context

Prepared 20 October 2026. This is a search specification and query-sequence handoff, not a claim that the proposed expanded search has already been run.

## 1. What we are trying to recover

The target is a **multi-core repeat array with informative intervening sequences, supported by a conserved RNA-related protein neighborhood**. Retrieve both close sequence homologs and more divergent architectural analogs. Protein-context similarity alone is insufficient.

The original E15 locus has seven related 12-nt cores and six intervening intervals next to a predicted 469-aa Ro60-like protein. Three inter-core regions have suggestive YrlA-like RNA-model matches. All six matches across the two known loci are below the recorded Rfam YrlA gathering threshold. Use the name **“Ro60-associated repeat array with YrlA-like segments”**, rather than treating YrlA identity, transcription or Ro60 binding as established.

| Reference | Original E15 | Known related locus |
|---|---|---|
| Candidate | C00938 | C03482 |
| Accession | NC_048046.1 | NC_019406.1 |
| Detector candidate, 0-based half-open | 116042–117330 | 108498–109813 |
| Array envelope, first core through last core | **116153–117140** | **108611–109601** |
| Array length | **987 nt** | **990 nt** |
| Repeat cores | 7 × 12 nt | 7 × 12 nt |
| Inter-core lengths | 163, 147, 118, 185, 158, 132 nt | 163, 147, 118, 187, 160, 131 nt |
| Ro60-like CDS | 114589–115999, negative strand | 107046–108456, negative strand |
| Protein length | 469 aa | 469 aa |
| Current architecture IDs | V01737 / M03143 | V01737 / M03143 |

The array envelope is not an experimentally defined transcript boundary. The DNA FASTA is supplied in the original genomic forward orientation. Separate reverse-complement and RNA-alphabet files are explicitly labeled. The Ro60-like protein and all six YrlA-like model hits are on the negative strand.

## 2. Query files supplied

All links below are relative to this Markdown file. Copy the `E15_search_queries` folder together with this file when handing the search to another analyst.

| File | Contents and intended use |
|---|---|
| [E15_NC048046_array.fna](E15_search_queries/E15_NC048046_array.fna) | The **original 987-nt E15 array**, one DNA record; the primary nucleotide query. |
| [E15_arrays_genomic_forward.fna](E15_search_queries/E15_arrays_genomic_forward.fna) | Both known complete array envelopes; query both to reduce reference bias. |
| [E15_arrays_reverse_complement.fna](E15_search_queries/E15_arrays_reverse_complement.fna) | Reverse complements of those arrays, for orientation-aware inspection. Do not count these as additional biological queries. |
| [E15_candidate_envelopes.fna](E15_search_queries/E15_candidate_envelopes.fna) | The wider detector candidates; useful for boundary sensitivity checks. |
| [E15_array_plus_3kb_flanks.fna](E15_search_queries/E15_array_plus_3kb_flanks.fna) | Each array with 3-kb flanks, including the Ro60-like CDS and surrounding context. |
| [E15_repeat_cores.fna](E15_search_queries/E15_repeat_cores.fna) | All 14 observed cores, retaining locus and coordinate provenance. |
| [E15_unique_core_seeds.fna](E15_search_queries/E15_unique_core_seeds.fna) | Five distinct observed 12-nt seeds; search both strands. |
| [E15_intercore_intervals.fna](E15_search_queries/E15_intercore_intervals.fna) | Twelve intervals, with both bounding cores excluded. Useful for segment-level homology searches. |
| [E15_Ro60_like_proteins.faa](E15_search_queries/E15_Ro60_like_proteins.faa) | Two full-length protein anchors with global IDs and genomic coordinates. |
| [E15_YrlA_like_segments.fna](E15_search_queries/E15_YrlA_like_segments.fna) | Six below-threshold RNA-model hit sequences in negative-strand model orientation, written as DNA. |
| [E15_YrlA_like_segments_RNA.fasta](E15_search_queries/E15_YrlA_like_segments_RNA.fasta) | The same six sequences with U replacing T. These are **not measured transcripts**. |
| [self_hit_exclusions.tsv](E15_search_queries/self_hit_exclusions.tsv) | Source records, known peer and previously identified accession aliases. |
| [sequence_metadata.json](E15_search_queries/sequence_metadata.json) | Coordinates, sequence lengths, interpretation and SHA256 for every record. |
| [MANIFEST.sha256](E15_search_queries/MANIFEST.sha256) | File checksums. |

The seven forward-strand cores in the original E15 array are:

```text
Core  Start0  End0    Sequence
R1    116153  116165  CTACGCCACCCT
R2    116328  116340  CTGCGCCACCCT
R3    116487  116499  CTACGCCACCCT
R4    116617  116629  CTACGCCACCAT
R5    116814  116826  CTACGCCACCCA
R6    116984  116996  CTACGCCACCTG
R7    117128  117140  CTACGCCACCCA
```

These short cores may be parts of larger repeat units. A match to one 12-mer is only a seed for investigation.

## 3. Search design: three complementary entry points

### A. Full-array and interval sequence search

Search both full array envelopes against the chosen bacterial and phage genome collections using a nucleotide local-alignment method such as BLASTN. Also search the twelve inter-core intervals independently. Preserve all useful HSPs, both strands, query and target coordinates, aligned lengths, identities, scores and database-size provenance.

For the complete array, aggregate **non-overlapping query coverage** from collinear HSPs at the same locus; do not sum overlapping HSP lengths. Inspect whether aligned blocks preserve the core/interval order. A high-identity 20–40 nt match should not be called a homolog of the entire 987-nt array. Report full-array homology separately from individual RNA-like segment matches and core-only matches.

Use the original masking settings as the main analysis and retain any relaxed-masking search as a labeled sensitivity check. The sequences are repeat-rich, so relaxed searches require stronger contextual review. Search both complete arrays and intervals because the whole-array query can miss loci that retain only some RNA-like modules.

This route finds sequence relatives; it does not establish RNA structure or function.

### B. Protein-first search, followed by local repeat discovery

1. Search the two Ro60-like proteins against a **complete target proteome**, with globally unique IDs and a protein-to-CDS-coordinate map. The earlier exported candidate-neighborhood proteome is insufficient for an unbiased whole-genome search.
2. Retain full-length and partial hits as separate classes. Inspect TROVE/Ro60-related domain support; a generic vWA-domain match alone is insufficient. Preserve alternatives rather than declaring function from one best hit.
3. Extract the protein locus with, initially, 10 kb of sequence on each side. Adjust to record boundaries; wrap circular records correctly and deduplicate wrapped hits. Search both sides even though the known E15 orientation suggests one side.
4. Search these windows for recurring short cores and array architecture, including loci that never passed the earlier information-score detector. The old candidate list must not define the search universe.
5. Recover neighboring proteins and compare family membership, order, strand and repeat-relative location. Normalize a global reverse complement before comparing architectures.

For example, after preparing `target_proteins.faa` and its coordinate map, a local DIAMOND starting command is:

```bash
diamond makedb --in target_proteins.faa --db target_proteins
diamond blastp \
  --query E15_search_queries/E15_Ro60_like_proteins.faa \
  --db target_proteins --very-sensitive --evalue 1e-5 \
  --max-target-seqs 0 --threads 8 \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore \
  --out Ro60_similarity_hits.tsv
```

This is a retrieval starting point, not a validated Ro60 classification threshold. Inspect alignment coverage and domain architecture downstream. Record tool/database versions; verify supported options with the installed executable. The example has not been executed as a new database search in this handoff. See the [official DIAMOND command-line documentation](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options).

For more distant candidates, profile/domain-based searches can supplement pairwise similarity. Build new profiles only after assembling a curated, nonredundant alignment; the two closely related E15 proteins alone provide limited family diversity.

### C. Array-first discovery, followed by protein-context comparison

This route allows both the repeat sequence and protein sequence to diverge.

1. Scan genome windows for repeated, sufficiently complex short seeds. Include the five observed E15 seeds, but also allow de novo core discovery so the search is not locked to `CTACGCCACCCT`.
2. For an initial E15-seeded pass, allow up to two substitutions in a 12-nt seed. Keep this as a labeled retrieval setting, not evidence of homology. Search direct repeats on both strands; label inverted-repeat structures separately.
3. Collapse duplicate and overlapping seed hits that describe the same physical core. Merge shifted hits to one core event using alignment quality, rather than inflating copy count with several overlapping k-mers.
4. Group physical core events into arrays and extract every intervening sequence. Inspect whole repeat neighborhoods, not just the core consensus.
5. Extend the aligned sequence around each core to determine whether a longer unit is conserved. Do not equate the discovery seed length with the true repeat length.
6. Annotate nearby proteins and compare architectures. Preserve strong arrays even when the nearest protein lacks an assigned function.

An initial **review-priority screen**, inherited from the current local analysis, is:

```yaml
minimum_nonoverlapping_cores: 4
minimum_core_length_nt: 8
median_leave_one_out_core_identity: 0.80
minimum_consensus_base_entropy_bits_per_base: 1.5
minimum_median_intercore_gap_nt: 20
maximum_intercore_gap_nt: 1000
minimum_fraction_of_gaps_at_least_10_nt: 0.70
maximum_array_CDS_union_fraction: 0.25
regular_spacing_priority:
  minimum_every_gap_nt: 20
  maximum_gap_CV: 0.60
retain_three_core_arrays_as_separate_leads: true
```

These are adjustable screening rules, not biologically universal cutoffs. Store near-threshold hits separately. Compute CDS coverage from the union of actual coding intervals, rather than summing overlapping CDS. Short predicted ORFs may be spurious; the original E15 includes a 41-aa ORF, so any-CDS-overlap exclusion would discard a known positive example.

## 4. Combine repeat evidence with a permissive gene architecture

Use two levels of context support:

- **One-anchor leads:** a credible Ro60-like protein plus a strong nearby array. Keep these discoverable even if the other flanking proteins are absent or unrecognizable.
- **Multi-anchor architecture support:** two or more sequence families retain their relative order and direction around the repeat. Allow peripheral gene insertion, deletion and replacement. Treat dUTPase-like or other neighboring genes as context evidence, not automatically as required components.

When a shared single-copy gene changes from one side of the repeat to the other, flag a different repeat slot instead of merging automatically. This rule keeps the adjacent C00939 locus separate from E15. Avoid merging by transitive chains where A matches B and B matches C but A and C lack direct architectural support.

Do not require 50% protein identity and 80% coverage as the only route to discovery. Those settings defined the earlier conservative sequence clusters; more remote candidates need separate profile and architecture evidence.

Rank by **repeat quality first**, then protein-context support and evolutionary diversity. A large protein-context group with no coherent array, such as M00033, should not outrank a small, well-supported array family. Do not count multiple anchor-pair descriptions of the same loci as independent systems.

## 5. RNA-family annotation and structure checks

Annotate candidate windows with a broad Rfam/Infernal screen, including YrlA, classical Y RNA and tRNA, and use independent tRNA-specific annotation when practical. Preserve known CRISPR, ART, VIPR and known-RNA matches as labeled controls rather than deleting them before comparison.

Run family-threshold and exploratory searches separately. For a downloaded YrlA model and extracted candidate windows, examples are:

```bash
# Curated family-threshold annotation
cmsearch --cpu 8 --cut_ga --tblout YrlA_GA.tbl RF02553.cm candidate_windows.fna > YrlA_GA.out

# Exploratory below-threshold search on focused windows only
cmsearch --cpu 8 --max -T 10 --tblout YrlA_exploratory.tbl RF02553.cm candidate_windows.fna > YrlA_exploratory.out
```

Do not reuse an arbitrary database-size setting from a previous run. Record the actual searched collection and any `-Z` value, especially for partitioned searches. Official workflow guidance is available in the [Rfam genome-annotation documentation](https://docs.rfam.org/en/latest/genome-annotation.html).

The 19 October model snapshot used YrlA GA=69 bits. All six E15 hits scored only 29.9–42.3 bits. **Requiring YrlA GA as a mandatory discovery filter would reject the original E15 examples.** Keep below-threshold resemblance as provisional evidence, while keeping family assignment distinct.

For RNA-like intervals:

- Inspect model coverage, boundary placement, strand and overlap with short cores/CDS.
- Compare model-guided pairing with unconstrained folds; compatibility is not independent validation.
- Align corresponding intervals from more divergent genomes. Two related genomes are insufficient for a strong covariation claim.
- Check whether similarity is confined to generic tRNA-like substructures.
- Keep transcription, processing, binding and biological activity as untested until directly demonstrated.

The original source tables had no tRNA-overlap flag, but a later focused search found a 167-nt tRNA-model hit in NC_019406.1 at 109216–109383 on the negative strand, 29.3 bits against GA=29. Preserve this ambiguity. The unflagged status in the original package is not proof of absence of known RNA.

## 6. Exclude self-hits and avoid overstating independent discoveries

Keep two distinct output sets:

1. **Positive-control recovery:** original E15 and the known peer, retained to verify the search recovers its starting examples.
2. **New discovery set:** excludes both source genomes, their known aliases and equivalent source records.

Resolve accession versions and assembly/record aliases. At minimum, exclude the source identities listed in `self_hit_exclusions.tsv`: `NC_048046` / `MH588545` and `NC_019406` / `JX100810`. Version changes do not create a new discovery. Compare source-record and locus-sequence hashes to catch renamed or repackaged identical data. A short core alone is not a reliable duplicate identifier.

Near-identical genomes should be grouped as close relatives, with one representative used for discovery counts. Report distinct loci, accessions and genome clusters separately. An identical locus in a demonstrably different genomic context can be recorded as an identical-sequence occurrence, but supplies no new sequence-divergence evidence. Multiple HSPs, two query orientations, or both known phage queries hitting the same bacterial locus are not independent discoveries.

For bacterial matches, report alignment length and both query/subject coverage. The earlier 42-nt tRNA-Lys match covered only 28.6% of I2; it did not establish a bacterial homolog of the whole interval or array, nor a tRNA-targeting mechanism.

## 7. Controls, outputs and acceptance criteria

Use matched controls: Ro60-like loci without an array; strong arrays lacking Ro60; known tRNA/Yrl/CRISPR/ART/VIPR loci; and genomic windows matched for length, base composition and coding context. Randomized sequence controls are useful but do not replace genomic backgrounds. If claiming enrichment, specify a valid genome-level denominator, account for related genomes and apply the same detector to controls.

Every retained locus should have:

```text
locus_id, accession.version, genome_cluster_id, source_database_version
array_start0, array_end0, orientation, core_coordinates, core_sequences
core_count, inferred_core_length, intercore_lengths, core_consistency
array_CDS_union_fraction, RNA_overlap_details, repeat_model_id
anchor_protein_ids, protein_alignment_identity_and_coverage, domain_hits
gene_order, gene_strands, repeat_relative_gene_roles
self_hit_status, duplicate_status, known_system_labels
full_array_query_coverage, interval_homology, RNA_model_scores_and_thresholds
support_tier, uncertainty_notes, source_files_and_SHA256
```

Assign support tiers rather than a single “new system” label:

| Tier | Required evidence |
|---|---|
| A: close sequence homolog | Full-array or substantial collinear interval homology plus an array; self-hits removed. |
| B: architectural analog | Strong array plus credible protein-family/order/strand/position support; limited nucleotide homology is explicitly reported. |
| C: exploratory lead | Three cores, partial/remote protein support, or uncertain RNA resemblance; retained for review, not pooled with stronger calls. |
| Control / known explanation | Known-RNA or known-repeat evidence retained with its annotation details. |

Deliver a table, per-locus FASTA, side-by-side architecture figures, core alignments and interval comparisons. Distinguish raw computational screening from loci actually inspected visually. Report null results only within the searched database and detection settings.

**Success means finding a repeat-supported, non-self locus whose sequence and architecture can be inspected—not merely increasing the number of protein-neighborhood matches.**

## 8. Evidence snapshot and limitations

The supplied sequences were extracted again from the accession-version references. All fourteen core coordinates and both 469-aa protein translations were checked. The RNA-alphabet file is a representation of predicted model-hit sequences only. The two known loci remain the only E15 members in the current v4 architecture group; no broader search was launched for this sequence handoff.

See the [English E15 evidence report](YRLA_ARRAY_REPORT_EN.html) for the six subthreshold YrlA hits, structural comparisons, the later tRNA-like ambiguity and source-file provenance. The report separates literature precedent for Ro60–RNA systems from evidence specific to E15.
