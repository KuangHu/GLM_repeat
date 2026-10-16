# Source and review bundles

## Standalone explainer v0.3 review packet

[`GLM_repeat_explainer_v0.3_review.zip`](GLM_repeat_explainer_v0.3_review.zip)
contains all eight public original FASTA/GenBank records, frozen candidates before
and after CDS annotation, all v0.3 outputs, 17 semantic tests, the unchanged
random10/review selection and all regeneration dependencies. No repository
checkout, model, GPU, genome download or detector inference is required after
installing the declared CPU dependencies.

Download the ZIP and `SHA256SUMS`. Verify the ZIP's line in SHA256SUMS, extract it
to an empty directory and run there:

```bash
python -m pip install -r environment/requirements-tested.txt
python scripts/verify_explainer_packet.py --run-tests --rerun rerun
```

This verifies 57 exact-byte manifest entries, checks all source DNA against
GenBank, rebuilds CDS annotations, all 15 controls and 467 heldout candidates,
engineering validation and the review set, and compares the evidence tables.
The verifier reports exactly which provenance-only fields are allowed to change.
Frozen output tables are not overwritten. See the ZIP's README for individual
stage commands, source inventory and limitations. A complete isolated rerun is
recorded in [`../results/packaging/explainer_v03_standalone.json`](../results/packaging/explainer_v03_standalone.json).

The manifest hashes archived bytes without LF/CRLF normalization. ZIP extraction
preserves those bytes; do not resave files with newline conversion before checking.
The v0.2 packet's updated manifest instead matches its Git blobs, which include
both LF and CRLF files; the Windows checkout can still require Git-aware checking.
The v0.2 explainer code, output and findings are unchanged by that packaging fix.

## Cluster source bundle

`GLM_repeat_region_scan_v0.1.tar.gz` contains the CDS-guided scan workflow,
explainer v0.3, CPU tests, frozen pilot threshold metadata and Chinese runbook.
It excludes model weights, raw genomes and GPU score tracks.

Download the archive and `SHA256SUMS`, then on the cluster:

```bash
sha256sum --ignore-missing -c SHA256SUMS
tar -xzf GLM_repeat_region_scan_v0.1.tar.gz
cd GLM_repeat_region_scan_v0.1
```

Follow `REGION_SCAN_zh.md` to configure input, the existing GPU environment and
the reviewed model snapshot. The archive also contains `MANIFEST.sha256` for
verifying every packaged file. CPU validation passed; new GPU inference and
Slurm execution have not been run for this version.
