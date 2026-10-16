#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH -J glmrep_arr
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%A_%a.out
# One FASTA record per array task (both orientations).
# usage: sbatch --array=1-N slurm/score_array_lrc.sbatch <fasta> <record_list> <out_root> <policy> [mode args...]
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
FASTA=$1; LIST=$2; ROOT=$3; POLICY=$4; shift 4
REC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -n "$REC" ] || { echo "no record for task $SLURM_ARRAY_TASK_ID"; exit 1; }
hostname; nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
"$GLMREP_PYTHON" scripts/timing_score.py --snapshot models/glm2_650M \
  --fasta "$FASTA" --records "$REC" --out "$ROOT/$REC" --batch-size 1 \
  --both-strands --attention-mask-policy "$POLICY" --device cuda:0 \
  --trust-reviewed-code "$@"
