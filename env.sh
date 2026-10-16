# source this in shells and sbatch jobs before using the glm_repeat env
export PROJECT=/global/home/users/kh36969/tools/GLM_repeat
export GLMREP_ENV=$HOME/.conda/envs/glm_repeat
export GLMREP_PYTHON=$GLMREP_ENV/bin/python
export PYTHONNOUSERSITE=1            # keep ~/.local site-packages out of this env
export HF_HOME=/global/scratch/users/kh36969/GLM_repeat/hf_cache
export TMPDIR=/global/scratch/users/kh36969/GLM_repeat/tmp
