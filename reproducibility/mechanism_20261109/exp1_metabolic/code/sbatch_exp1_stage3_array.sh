#!/bin/bash
#SBATCH --job-name=exp1_st3_arr
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --array=0-7
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_arr_%A_%a.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_arr_%A_%a.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
CODE_DIR="${BASE}/code_v2"
export TMPDIR="${BASE}/tmp/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results/shards" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 3 Array Shard ${SLURM_ARRAY_TASK_ID}/8 on $(hostname)"

python3 "${CODE_DIR}/stage3_full484_worker.py" \
  --inventory-484-tsv "${BASE}/results/exp1_model_inventory_and_qc_484.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-dir "${BASE}/results" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --num-tasks 8 \
  --n-target-solutions 20 \
  --milp-time-limit-s 6.0 \
  --workers 32

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 3 Array Shard ${SLURM_ARRAY_TASK_ID}/8"
