#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 05:00:00
#SBATCH -J glmrep_region
#SBATCH -o logs/%x_%A_%a.out
# Submit from repository root. See REGION_SCAN_zh.md for environment and arguments.
set -euo pipefail
: "${GLMREP_PYTHON:?Set absolute path to the existing GPU Python}"
: "${GLMREP_SNAPSHOT:?Set absolute path to the reviewed gLM2 snapshot}"
PREPARED=$1
OUT=$2
N=$3
MODE=${4:-random_partition}
# Slurm copies the batch script into its spool directory, so ${BASH_SOURCE[0]}
# does NOT point at the repository. Prefer the submission directory (the guide
# says to submit from the repository root); allow an explicit override; fall back
# to the script path only for direct non-Slurm execution.
ROOT=${GLMREP_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}
if [ ! -f "$ROOT/scripts/score_region_chunk.py" ]; then
  echo "scripts/score_region_chunk.py not found under ROOT=$ROOT" >&2
  echo "Submit from the repository root, or set GLMREP_REPO_ROOT." >&2
  exit 2
fi
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
"$GLMREP_PYTHON" "$ROOT/scripts/score_region_chunk.py" \
  --prepared "$PREPARED" --out "$OUT" --snapshot "$GLMREP_SNAPSHOT" \
  --chunks "$N" --task "${SLURM_ARRAY_TASK_ID:?Use --array=1-N}" --mode "$MODE" \
  --trust-reviewed-code
