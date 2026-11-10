#!/bin/bash
#SBATCH --job-name=exp3_s1_col
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64000M
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_collect_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_collect_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s1_collect_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 1 CRISPR collect on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 "$EXP3_ROOT/code/stage1_crispr_recovery.py" --mode collect
rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE1_COLLECT_DONE"
echo "[$(date -Is)] Completed Stage 1 CRISPR collect."
