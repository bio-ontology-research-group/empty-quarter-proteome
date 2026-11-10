#!/bin/bash
#SBATCH --job-name=exp3_st4_bldref
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=350G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_build_ref_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_build_ref_%j.err

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 4 Build Full Unified Competitive Reference Job $SLURM_JOB_ID on $(hostname)"
echo "Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/stage4_build_combined_competitive_reference.py" --threads 64

echo "Completed at: $(date)"
