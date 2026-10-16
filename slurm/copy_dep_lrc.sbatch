#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:30:00
#SBATCH -J glmrep_copydep
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%j.out
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
hostname; nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
"$GLMREP_PYTHON" scripts/copy_dependence.py score --out "$1" \
  --snapshot models/glm2_650M --device cuda:0 --batch-size 8 \
  --attention-mask-policy always --trust-reviewed-code
