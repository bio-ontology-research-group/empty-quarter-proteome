#!/bin/bash
#SBATCH --job-name=s4_v3stat
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v3stat_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v3stat_%j.err

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 4 v3 Statistical Audit Job ${SLURM_JOB_ID:-local} on $(hostname)"
echo "Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/scratch_tmp"
mkdir -p "${TMPDIR}" "${EXP3_ROOT}/logs" "${EXP3_ROOT}/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

cd "${EXP3_ROOT}"

echo "[Step 1] Running stage4_statistical_audit_v3.py..."
python "${EXP3_ROOT}/code/stage4_statistical_audit_v3.py" --exp3-root "${EXP3_ROOT}" --n-perms 9999 --seed 20261109

echo "========================================================================"
echo "Stage 4 v3 Statistical Audit Completed at: $(date)"
echo "========================================================================"
