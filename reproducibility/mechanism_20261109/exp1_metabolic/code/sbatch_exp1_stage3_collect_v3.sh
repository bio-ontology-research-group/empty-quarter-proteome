#!/bin/bash
#SBATCH --job-name=exp1_st3_col_v3
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_v3_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st3_col_v3_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"

echo "[$(date -Is)] Starting exp1_metabolic Stage 3 Collector (v3) on $(hostname)"

python3 "${BASE}/code_v3/stage3_collect_and_summarize_v3.py" \
  --inventory-484-tsv "${BASE}/results/exp1_inventory_484_deduplicated.tsv" \
  --shards-dir "${BASE}/results/shards" \
  --masking-controls-tsv "${BASE}/results/exp1_masking_controls_7hq_v4.tsv" \
  --out-dir "${BASE}/results" \
  --num-tasks 8

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 3 Collector (v3 PASSED)"
