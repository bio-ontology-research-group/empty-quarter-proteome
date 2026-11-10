#!/bin/bash
#SBATCH --job-name=exp1_cf2_pilot
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2/logs/exp1_cf2_pilot_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2/logs/exp1_cf2_pilot_%j.err

set -euo pipefail

CLOSURE_DIR="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure_v2"
mkdir -p "${CLOSURE_DIR}/tmp" "${CLOSURE_DIR}/logs"
export TMPDIR="${CLOSURE_DIR}/tmp"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

python3 "${CLOSURE_DIR}/code/run_targeted_candidate_fva_closure.py" \
  --mode pilot \
  --results-dir /scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results \
  --code-v2-dir /scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/code_v2 \
  --out-dir "${CLOSURE_DIR}"
