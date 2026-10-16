#!/bin/bash
#SBATCH -A pc_rubinlab
#SBATCH -p es1
#SBATCH -q es_debug
#SBATCH --gres=gpu:A40:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH -J glmrep_smoke
#SBATCH -o /global/scratch/users/kh36969/GLM_repeat/logs/%x_%j.out
set -euo pipefail
source /global/home/users/kh36969/tools/GLM_repeat/env.sh
cd "$PROJECT"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
FASTA=${1:-data/toy/smoke_short.fa}
OUT=results/smoke_${SLURM_JOB_ID}
hostname; nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
"$GLMREP_PYTHON" -c "import torch,transformers;print('torch',torch.__version__,'cuda',torch.version.cuda,'avail',torch.cuda.is_available(),'bf16',torch.cuda.is_bf16_supported());print('transformers',transformers.__version__)"
"$GLMREP_PYTHON" scripts/reference_score.py \
  --snapshot models/glm2_650M --fasta "$FASTA" --out "$OUT" \
  --mode exact --batch-size 1 --device cuda:0 --both-strands --trust-reviewed-code
"$GLMREP_PYTHON" - "$OUT" <<'PY'
import json,sys,numpy as np,pathlib
o=pathlib.Path(sys.argv[1]); m=json.loads((o/'run.json').read_text())
for r in m['records']:
    v=np.load(o/r['track'])['above_uniform']
    print(r['id'],r['scoring_orientation'],r['length'],f"{r['seconds']:.1f}s",
          'mean %.3f min %.3f max %.3f'%(v.mean(),v.min(),v.max()))
print('peak GPU MB', __import__('torch').cuda.max_memory_allocated()//2**20)
PY
