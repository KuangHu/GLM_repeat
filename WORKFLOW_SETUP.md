# GLM-repeat: workflow and agent implementation specification

**Version:** 0.1, 2026-10-14  
**Task:** raw DNA -> evidence of repetitive sequence organization.  
**Backbone:** `tattabio/gLM2_650M`.  
**Status:** setup specification plus a reference scorer, not a finished genome-wide detector.

## 0. Scope and first milestone

Build ONE DNA-only repeat-sensitive model/workflow. CRISPR, VIPR and ART are separate validation families, not three hard-coded search modes. The eventual discovery scan must not require Cas, RT, protein anchors, an ncRNA prediction, a known motif, a palindrome, a target match, or a specific repeat length.

Include simple repeats and homopolymers in the repeat-positive concept. Detecting repetitive organization does **not** establish functional RNA, a new molecular system, or genome-editing activity. Treat ART as an array-associated RT system; do not label its individual RNAs as experimentally established msr-msd units.

First milestone:

```text
Pinned gLM2 + correct nucleotide/token coordinates
    -> exact masked scoring on a small benchmark
    -> local predictability tracks and calibrated window ranking
    -> separate CRISPR / VIPR / ART performance report
    -> validated faster scoring on a 10-Mb pilot
```

Do not begin with end-to-end fine-tuning or the full genome collection. Train a small head only if the frozen-model scores are inadequate.

## 1. Integration boundaries

Suggested independent project location:

```text
/global/home/users/kh36969/tools/GLM_repeat
```

Use a separate environment, for example `$HOME/.conda/envs/glm_repeat`. Do not install Torch into `opfi`, modify the existing Minerva environment, or change the Protein_ncRNA caller, decoys, thresholds, or outputs. Existing genome indexes may be read without modification. All new output fields use `glmrep_`.

Run model inference and substantial preprocessing on Slurm compute nodes, not the login node. Reuse the site's working scratch conventions; test scratch read/write responsiveness before a large download. One independent model process per allocated GPU is the default. Do not assume a partition/account name or eight free A40 GPUs.

## 2. What is supplied versus what the agent must implement

| Supplied and locally checked | Still to implement or validate on the cluster |
|---|---|
| Model download/pinning helper | Actual checkpoint download and source-code review |
| Exact and random-partition scoring primitive | GPU checkpoint-loading smoke test |
| Strict DNA/token mapping and four-base normalization | Real benchmark sequence acquisition and provenance |
| Synthetic fixture generator | Whole-genome windowing, aggregation and calling |
| Default configuration and Slurm reference job | Calibration and blind-test reports |
| 16 CPU tests using dummy models | Optional classification head and its training |

The helper does not consume `config/defaults.yaml` automatically. The agent must implement the configuration-driven orchestrator; until then, use the helper's explicit CLI flags. Do not advertise unimplemented commands as available.

## 3. Environment setup

Official gLM2 uses lowercase nucleotide tokens, mixed protein/DNA pretraining, a masked language-model objective and a 4,096-token context [S1]. This workflow deliberately tests DNA-only input, including coding DNA; that is an application choice, not an assurance that all genome regions are equally well modeled.

A conservative **candidate** bootstrap is shown below. It is not a GPU-validated lockfile. PyTorch 2.5.1 CUDA 12.1 wheels are documented by PyTorch [S4]; the chosen CUDA wheel still has to match the cluster driver.

```bash
export PROJECT=/global/home/users/kh36969/tools/GLM_repeat
export GLMREP_ENV="$HOME/.conda/envs/glm_repeat"
# Unpack this bundle into PROJECT before running the following commands.
cd "$PROJECT"
conda create -y -p "$GLMREP_ENV" python=3.11 pip
export GLMREP_PYTHON="$GLMREP_ENV/bin/python"

# Run setup/download on a site-approved setup or transfer node.
# Check the cluster driver before choosing this CUDA wheel.
"$GLMREP_PYTHON" -m pip install torch==2.5.1 \
  --index-url https://download.pytorch.org/whl/cu121
"$GLMREP_PYTHON" -m pip install -r requirements-bootstrap.txt
"$GLMREP_PYTHON" -m pip check
"$GLMREP_PYTHON" -m pytest -q tests
```

After the first real GPU smoke test passes, record `pip freeze`, Python version, driver, CUDA runtime, GPU model, model revision and source hashes in `results/environment/`. Freeze that working environment. An installation failure is not permission to change existing project environments.

### Download, pin, inspect

```bash
"$GLMREP_PYTHON" scripts/download_model.py \
  --model tattabio/gLM2_650M --out models/glm2_650M
```

This resolves `main` to a full commit SHA and writes a file-hash manifest. It downloads weights/configuration/code but does not execute the downloaded Python. Review the model code before granting `trust_remote_code`. Use that local snapshot offline during jobs; never allow a running campaign to float with `main` [S5].

For likelihood scoring use **AutoModelForMaskedLM**, not the README's embedding-only AutoModel example. Verify that the pretrained MLM head actually loaded; missing or newly initialized head weights are a fatal error [S2]. A successful encoder load alone does not establish a working scorer.

## 4. Inputs and benchmark separation

Genome input: FASTA plus a manifest of `genome_id`, path, accession/version, assembly identity and source checksum. Benchmark input: true genomic DNA sequences plus an annotation manifest:

```text
record_id  family  locus_group  accession_version  contig
start0  end0  reference_source  split  sequence_sha256
```

`family` is CRISPR, VIPR, ART, simple_repeat, background or synthetic. It is metadata only and must never enter model input. All windows, both orientations, flanks, perturbations and related copies of one locus stay in the same split. Split by locus family/lineage before window extraction, not by random windows. Use separate development, calibration and locked-test groups. Report limited independent sample sizes honestly.

Use authentic sequences with traceable accessions or author data for the real benchmark. Do not reconstruct biological positives from a PDF logo, guess accessions, or substitute this bundle's toys for real data. The uploaded ART manuscript is a source of methods and locus leads, not a sequence file. Its 28 array-bearing loci are distinct from its larger RT collection [S3]. If a sequence cannot be obtained, record `missing_reference` and do not claim that family was tested.

Include ordinary coding and intergenic DNA, multiple GC backgrounds, and bacterial/phage/plasmid contexts. Random genomic windows can contain undiscovered repeats: call them **background**, not guaranteed negatives. A smaller manually assessed non-repeat subset can support a true classification-error estimate without making exhaustive negative curation a prerequisite.

### Engineering fixtures

```bash
"$GLMREP_PYTHON" scripts/make_toys.py --out data/toy/toys.fa
```

These are synthetic architecture fixtures, not CRISPR/VIPR/ART biological benchmarks. Use a short record first for the model-loading smoke test; the complete toy set can then exercise longer contexts.

## 5. DNA normalization and windowing

Use lowercase DNA and a `<+>` prefix. Validate every nucleotide's actual token position, and obtain token IDs from the tokenizer rather than hard-coding them. Count special tokens in the 4,096-token limit. Do not truncate silently [S1, S2].

Default genomic window: 3,000 nt; stride: 1,500 nt. Retain the final shorter window and short contigs. Split at every non-ACGT base for v1, preserving the original genomic coordinates; do not delete ambiguous bases and join unrelated stretches. Score valid segments of at least 40 nt, and report the unscored-base fraction and short fragments. Do not filter coding or low-complexity sequence.

Score the forward sequence and its actual reverse complement, each with `<+>`. Map reverse positions back to forward genomic coordinates. `scoring_orientation` is not an inferred transcription strand. In linear contigs, never join across a contig boundary. Circular-origin coverage is an explicit later option and must use circular metadata.

A >3-kb array may be covered by multiple windows. This can detect parts of it; it is not equivalent to observing the whole array in one model context. A window containing no second copy cannot demonstrate cross-copy predictability. Record these limitations rather than silently extending the model context.

## 6. Exact reference scoring: reproduce the ART measurement

For each nucleotide i, construct one variant with only that nucleotide replaced by the mask token. Batch variants for efficiency, but keep one hidden nucleotide per variant. Feed both left and right context. From the masked-position logits, keep only the lowercase a/c/g/t logits and normalize over these four values:

```text
log_p4(i) = log_softmax(logits(i, [a,c,g,t]))[true_base(i)]
v(i)      = log_p4(i) - log(0.25)
```

This is the nucleotide-level, natural-log "above uniform chance" measurement described for gLM2 in the ART manuscript [S3]. `v(i) <= ln(4)`. Negative values are valid. Do not normalize across protein/special tokens, shift targets as in an autoregressive model, or extract true-base probabilities from an unmasked input. The last mistake leaks the answer into the model.

Save the unsmoothed track and, for comparison with the manuscript, a 16-nt running mean. Do not use the smoothed track as the only stored result: it can obscure very short periodic patterns. These values measure masked predictability, not probabilities that a locus is repetitive or functional.

The supplied helper implements this primitive for **already windowed, ACGT-only** records:

```bash
# Run this inside an allocated GPU job, after reviewing the model code.
"$GLMREP_PYTHON" scripts/reference_score.py \
  --snapshot models/glm2_650M --fasta data/bench/windows.fa \
  --out results/reference_exact --mode exact --batch-size 1 \
  --device cuda:0 --both-strands --trust-reviewed-code
```

Output: two mapped NPZ tracks per record when both strands are requested, and `run.json`. The helper verifies snapshot hashes and rejects incomplete checkpoint loading. It does not yet emit genomic BED calls or a repeat classifier score.

## 7. Convert local evidence into a window ranking

Implement a small feature extractor, without a known-repeat finder or family labels. Start with raw-score quantiles and standard deviation, plus local block statistics. For each supported block size in 64, 128, 256, 512, 1024 and 2048 nt, slide by one quarter of the block length and retain maximum block mean and maximum mean of the top 20% of positions. For very short records use the available record length. Record valid block counts.

The top-fraction statistic is a proposed way to retain sparse conserved cores in long variable regions; it is not a published universal repeat classifier. Initially evaluate both maximum-block-mean and top-fraction ranking. Select the score definition on development data, then freeze it before calibration. Do not keep tuning on the locked real-family tests.

Take the maximum orientation-level window score across the two physical orientations, and calibrate that **same maximum**. Do not merge nucleotide-wise orientation maxima and describe the synthetic track as one RNA strand. Preserve both original tracks.

A local-GC/Markov-adjusted likelihood may be added as a diagnostic, but must not be the only decision score: subtracting compositional predictability can remove the very homopolymers the user asked to include. Predictable coding or structured single-copy DNA can also score highly; calibration must determine whether the frozen signal is discriminative enough.

## 8. Calibration and biological acceptance

Use calibration backgrounds processed with exactly the same lengths, orientations, scales, scoring mode and merging rules as the discovery scan. Determine operating thresholds at 1, 10 and 100 merged background hits/Mb. The provisional pilot workload target is 10 hits/Mb, not a known biological constant or an asserted false-positive rate.

On the locked test set, report CRISPR, VIPR and ART separately: number of independent locus groups, locus recovery, recovery of annotated repeat cores, boundary/coverage errors, and sensitivity across the workload thresholds. Also report simple-repeat recovery and background hit burden. Count a locus as loosely recovered when a call overlaps its annotated repeat architecture; report the overlap fraction too, so a one-base intersection cannot masquerade as full array reconstruction.

Use paired perturbations for mechanistic checks: preserve the rest of a real window while independently disrupting each repeat copy; preserve NN while disrupting the conserved phase in VIPR-style synthetic controls; change NN while preserving the conserved phase as a **positive-preserving** perturbation. Independent shuffles of a homopolymer cannot remove its repeat signal, so that is not a valid negative. A score decline on a finite sample is an empirical result to measure, not an engineering assertion to hard-code.

Do not call background hits/Mb an FDR or a calibrated biological probability. Some background hits are real repeats. Always inspect a small random sample from high, middle and low score bins, not just attractive known positives. An optional low-complexity or ordinary repeat baseline is for evaluation only, not a mandatory discovery filter.

## 9. Faster scoring: validate before using at scale

Exact scoring needs L masked variants for an L-nt window. It is a reference measurement, not the default for millions of windows.

The supplied `random_partition` mode is an **experimental approximation**:

```text
For each of two independent rounds:
  randomly permute all nucleotide positions
  divide into 20 balanced groups
  mask one group per variant
  score only the hidden positions
Average the two log probabilities obtained for each nucleotide.
```

For a 3,000-nt record this uses 40 variants rather than 3,000 per orientation. This is a variant-count reduction, NOT a measured 75-fold wall-clock speedup. Each variant masks about 5% of positions; the conditioning differs from exact single-site masking. It must have its own calibration.

Never mask every fifth/tenth/twentieth position as a fixed phase: that can hide corresponding repeat positions together. Seeds must be deterministic. Use matching partitions for a sequence and its paired controls of the same length. The supplied helper uses the same explicit seed; a production worker may derive it from a stable parent-example ID, never from Python's randomized `hash()`.

Validate approximate versus exact tracks on a balanced subset of all real families, simple repeats and background windows. Evaluate both score ranking and downstream locus recovery, not correlation alone. A provisional acceptance gate is no more than 5 percentage points of recall loss for any adequately sampled family at matched background workload. With tiny family counts, inspect misses rather than claiming precision unsupported by the data. If the gate fails, increase groups/rounds or keep exact mode for the pilot; do not claim equivalence.

```bash
"$GLMREP_PYTHON" scripts/reference_score.py \
  --snapshot models/glm2_650M --fasta data/bench/windows.fa \
  --out results/reference_fast --mode random_partition \
  --groups 20 --repeats 2 --batch-size 1 \
  --device cuda:0 --both-strands --trust-reviewed-code
```

## 10. Blind pilot and production orchestration

Only after biological and runtime checks, scan at most 10 Mb of genomes not used for calibration. Start with one GPU and microbatch 1; increase only after measuring peak memory and numerical stability. Allocate independent FASTA/window shards across Slurm GPU jobs, with one checkpoint resident per process. The helper itself is not multi-GPU.

Implement per-shard atomic output and resumability keyed by input checksum, window coordinates, model SHA, score mode, mask seed and configuration hash. A missing/failed shard must not appear as a negative call. Measure seconds per window, scored nt/s, GPU peak memory, bytes written, and extrapolated cost including both strands and overlaps before a full-corpus run.

For overlapping windows, preserve a deterministic center-preferred per-base track and separately retain any above-threshold local block call. Merge overlapping blocks on the same contig; do not join arbitrary low-score gaps. Report coarse model-supported intervals, not exact repeat boundaries. Calibration must use this same postprocessing. Long repeat-copy spacings may require a later grouping step and should not be hidden by aggressive merging.

Audit 100 randomly selected blind windows with the exact scorer, **including fast-score negatives**, to detect recall loss. Do not evaluate the exact scorer only on selected fast positives. Both scored and skipped bases enter coverage summaries.

Expected production outputs:

```text
results/environment/                    locked environment and source hashes
results/benchmark/metrics_by_family.tsv  counts, recall, workload curves
results/benchmark/report.md             misses and interpretation
results/pilot/glmrep_windows.tsv         window scores and QC
results/pilot/glmrep_candidates.bed      coarse 0-based half-open intervals
results/pilot/glmrep_tracks/             original per-base tracks
results/pilot/glmrep_run.json            provenance, failures and coverage
```

Suggested window fields: `genome_id, contig, start0, end0, scoring_orientation, mode, model_sha, score_raw, score_name, threshold_id, valid_nt, unknown_nt, edge_flag, window_sha256`. Leave `repeat_probability` absent. RNA function and novelty are outside this release.

## 11. Optional single-model classifier, only after the reference result

If masked-predictability ranking cannot separate repetitive organization from background at a useful workload, train a small head on a **frozen gLM2 encoder**. This remains one DNA-only model, not a protein-anchored pipeline or three family-specific detectors.

Use contextual per-nucleotide embeddings, a small projection/1D-convolution layer and masked local/attention pooling with a window-level binary head. Preserve position-level information; do not start with only a whole-window mean embedding. Using unmasked input for this supervised encoder task is legitimate; using it for the masked-likelihood reference task is not.

Generate varied synthetic direct, partially conserved and spaced repeats with random motifs, lengths, spacing, degeneration and copy counts in matched real-DNA backgrounds. Include simple repeats as positives. Split source backgrounds, motif families and paired examples before augmentation. Keep real CRISPR/VIPR/ART test families untouched, while acknowledging possible foundation-model pretraining exposure. Do not interpret held-out task training as proof of pretraining novelty.

Cache compressed/local features only after estimating disk needs; full L-by-hidden-dimension tensors can become the bottleneck. Establish a small pilot training set first. Calibrate the head on untouched data and compare to the exact reference and simple baselines. Its sigmoid alone is not a real-world biological posterior. Whole-genome adoption requires the same per-family and background-workload report.

## 12. Acceptance checklist for the coding agent

Engineering: correct lowercase IDs, mask position, four-base softmax, token mapping, tail-window handling, reverse mapping, deterministic seeds, complete mask coverage, no silent truncation, no gaps joined, no missing MLM weights, finite logits, independent environment, successful Slurm run.

Scientific: authentic provenance; separate biological families; controls and related windows kept in one split; simple repeats retained; matched background evaluation; explicit discovery-versus-calibration separation; no claimed recall before real data testing; approximate scoring separately validated.

Delivery: commit code/configuration; retain CPU test results; add actual GPU smoke logs and environment lock; provide family-wise metrics and failure examples; mark unavailable benchmarks as missing. First delivery should answer whether the frozen model provides a useful repeat-sensitive signal, not merely whether a neural network loaded.

## Sources and verification boundary

[S1] TattaBio official gLM2 repository/model card: modality, lowercase DNA, context and inference interface.

```text
https://github.com/TattaBio/gLM2
https://huggingface.co/tattabio/gLM2_650M
```

[S2] Official custom model/configuration/tokenizer sources: separate encoder and pretrained MLM classes, output logits and lowercase token vocabulary.

```text
https://huggingface.co/tattabio/gLM2_650M/blob/main/modeling_glm2.py
https://huggingface.co/tattabio/gLM2_650M/blob/main/configuration_glm2.py
https://huggingface.co/tattabio/gLM2_650M/blob/main/tokenizer.json
```

[S3] User-provided `Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf`: Methods p.39, genomic language models; Results p.6, array-bearing loci. The 2,900-nt L0050 measurement uses per-base masking, four-nucleotide normalization, log probability minus ln(0.25), and 16-nt smoothing. Generalization to all requested families is a proposed benchmark, not established by that example.

[S4] Official PyTorch previous-version installation instructions.

```text
https://pytorch.org/get-started/previous-versions/
```

[S5] Hugging Face custom-code loading documentation. Pin and review remote code; retain hashes.

```text
https://huggingface.co/docs/transformers/custom_models
```

**Local verification for this bundle:** 16 CPU engineering tests passed using dummy models. No pretrained gLM2 weights were loaded, no GPU inference was run, and no biological recall was measured here. The bootstrap dependency combination still needs its real cluster smoke test. The pooling features, random-partition approximation, calibration choices and optional head are proposed workflow components, not claims from the source papers.
