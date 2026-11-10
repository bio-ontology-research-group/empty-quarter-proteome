#!/bin/bash
#SBATCH --job-name=exp1_st3_col_v4
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_v4_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_v4_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"

python3 "${BASE}/code_v5/stage3_collect_and_summarize_v4.py" \
  --inventory-484-tsv "${BASE}/results/exp1_model_inventory_and_qc_484.tsv" \
  --shards-dir "${BASE}/results/shards" \
  --masking-validation-receipt "${BASE}/results/exp1_masking_controls_7hq_v4_validation_receipt.json" \
  --out-dir "${BASE}/results" \
  --num-tasks 8

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed v4 full484 collection with separately validated control usability"
