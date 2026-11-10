#!/bin/bash
#SBATCH --job-name=exp1_st3_col
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
CODE_DIR="${BASE}/code_v2"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 3 Collector on $(hostname)"

python3 "${CODE_DIR}/stage3_collect_and_summarize.py" \
  --inventory-484-tsv "${BASE}/results/exp1_model_inventory_and_qc_484.tsv" \
  --shards-dir "${BASE}/results/shards" \
  --out-dir "${BASE}/results" \
  --num-tasks 8

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 3 Collector"
