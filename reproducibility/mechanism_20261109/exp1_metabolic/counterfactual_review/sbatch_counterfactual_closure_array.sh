#!/bin/bash
#SBATCH --job-name=exp1_cf2_fva_arr
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2/logs/exp1_cf2_fva_arr_%A_%a.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2/logs/exp1_cf2_fva_arr_%A_%a.err

set -euo pipefail

CLOSURE_DIR="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2"
mkdir -p "${CLOSURE_DIR}/tmp/task_${SLURM_ARRAY_TASK_ID}" "${CLOSURE_DIR}/logs" "${CLOSURE_DIR}/shards"
export TMPDIR="${CLOSURE_DIR}/tmp/task_${SLURM_ARRAY_TASK_ID}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

python3 "${CLOSURE_DIR}/code/run_targeted_candidate_fva_closure.py" \
  --mode shard \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --num-tasks 8 \
  --workers 16 \
  --results-dir /scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results \
  --code-v2-dir /scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/code_v2 \
  --out-dir "${CLOSURE_DIR}"
