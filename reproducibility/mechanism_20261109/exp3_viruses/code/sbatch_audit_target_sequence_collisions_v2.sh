#!/bin/bash
#SBATCH --job-name=exp3_col_v2
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/audit_collision_v2_%j.log

set -euo pipefail

echo "========================================================================"
echo "EXP3 Target Sequence Collision Audit v2 Job $SLURM_JOB_ID on $(hostname)"
echo "Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/audit_target_sequence_collisions_v2.py"

echo "Completed at: $(date)"
