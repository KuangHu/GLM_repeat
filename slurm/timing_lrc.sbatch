#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -J glmrep_timing
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%j.out
# usage: sbatch slurm/timing_lrc.sbatch <fasta> <out_dir> <batch_size> [extra timing_score.py args]
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
FASTA=$1; OUT=$2; BS=$3; shift 3
hostname; nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
"$GLMREP_PYTHON" scripts/timing_score.py --snapshot models/glm2_650M \
  --fasta "$FASTA" --out "$OUT" --mode exact --batch-size "$BS" \
  --device cuda:0 --trust-reviewed-code "$@"
