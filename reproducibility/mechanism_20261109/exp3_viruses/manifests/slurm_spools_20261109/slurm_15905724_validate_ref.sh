#!/bin/bash
#SBATCH --job-name=exp3_val_ref
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --time=04:00:00
#SBATCH --mem=350G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/validate_ref_%j.log

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 4 Validate Reference Job $SLURM_JOB_ID on $(hostname)"
echo "Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/stage4_validate_competitive_reference.py"

echo "Completed at: $(date)"
