#!/bin/bash
#SBATCH --job-name=glmrep_ref
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=glmrep_%j.log
set -euo pipefail
: "${PROJECT:?Set PROJECT to the unpacked setup/project directory}"
: "${GLMREP_PYTHON:?Set GLMREP_PYTHON to the isolated environment Python}"
: "${GLMREP_SNAPSHOT:?Set GLMREP_SNAPSHOT to the reviewed pinned snapshot}"
: "${GLMREP_FASTA:?Set GLMREP_FASTA to a pre-windowed benchmark FASTA}"
: "${GLMREP_OUT:?Set GLMREP_OUT to a new output directory}"
: "${SLURM_JOB_ID:?GPU scoring must run in a Slurm job}"
cd "$PROJECT"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
"$GLMREP_PYTHON" scripts/reference_score.py \
  --snapshot "$GLMREP_SNAPSHOT" --fasta "$GLMREP_FASTA" \
  --out "$GLMREP_OUT" --device cuda:0 --mode exact \
  --batch-size 1 --both-strands --trust-reviewed-code
