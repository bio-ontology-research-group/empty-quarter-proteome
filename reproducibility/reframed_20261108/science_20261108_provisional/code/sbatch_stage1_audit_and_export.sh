#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --job-name=eq_sci_stg1
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
cd "$EQ_SCIENCE_OUT"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

echo "[Stage 1] Running unit test gate on compute node $(hostname) at $(date -Is)..."
python3 -m pytest -o cache_dir="$EQ_SCIENCE_OUT/logs/.pytest_cache" "$EQ_SCIENCE_OUT/tests/test_science_20261108.py" -v

echo "[Stage 1] Launching stage1_audit_and_export.py at $(date -Is)..."
python3 "$EQ_SCIENCE_OUT/code/stage1_audit_and_export.py"
echo "[Stage 1] Done at $(date -Is)."
