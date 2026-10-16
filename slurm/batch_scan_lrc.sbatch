#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_normal
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -J glmrep_bscan
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%j.out
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
FASTA=data/bench/GLM_REPEAT_REAL_SEEDS_V1/data/timing6_3kb.fa
for BS in 1 4 16 32; do
  "$GLMREP_PYTHON" scripts/timing_score.py --snapshot models/glm2_650M --fasta $FASTA \
    --records w0001 --out results/benchmark/bscan_${SLURM_JOB_ID}/bs$BS --mode exact \
    --batch-size $BS --device cuda:0 --trust-reviewed-code || echo "bs=$BS FAILED"
done
