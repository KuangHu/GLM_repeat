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
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
"$GLMREP_PYTHON" "$ROOT/scripts/score_region_chunk.py" \
  --prepared "$PREPARED" --out "$OUT" --snapshot "$GLMREP_SNAPSHOT" \
  --chunks "$N" --task "${SLURM_ARRAY_TASK_ID:?Use --array=1-N}" --mode "$MODE" \
  --trust-reviewed-code
