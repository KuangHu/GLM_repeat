#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 05:00:00
#SBATCH -J glmrep_chunk
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%A_%a.out
# usage: sbatch --array=1-N slurm/score_chunk_lrc.sbatch <fasta> <record_list> <out_root> <nchunks> <mode args...>
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
FASTA=$1; LIST=$2; ROOT=$3; N=$4; shift 4
mapfile -t ALL < "$LIST"
RECS=()
for ((i=SLURM_ARRAY_TASK_ID-1; i<${#ALL[@]}; i+=N)); do RECS+=("${ALL[$i]}"); done
echo "task $SLURM_ARRAY_TASK_ID: ${#RECS[@]} records"
hostname; nvidia-smi --query-gpu=name --format=csv,noheader
"$GLMREP_PYTHON" scripts/timing_score.py --snapshot models/glm2_650M \
  --fasta "$FASTA" --records "${RECS[@]}" --out "$ROOT/part$SLURM_ARRAY_TASK_ID" \
  --batch-size 1 --both-strands --device cuda:0 --trust-reviewed-code "$@"
