# Cluster source bundle

`GLM_repeat_region_scan_v0.1.tar.gz` contains the CDS-guided scan workflow,
explainer v0.3, CPU tests, frozen pilot threshold metadata and Chinese runbook.
It excludes model weights, raw genomes and GPU score tracks.

Download the archive and `SHA256SUMS`, then on the cluster:

```bash
sha256sum -c SHA256SUMS
tar -xzf GLM_repeat_region_scan_v0.1.tar.gz
cd GLM_repeat_region_scan_v0.1
```

Follow `REGION_SCAN_zh.md` to configure input, the existing GPU environment and
the reviewed model snapshot. The archive also contains `MANIFEST.sha256` for
verifying every packaged file. CPU validation passed; new GPU inference and
Slurm execution have not been run for this version.
