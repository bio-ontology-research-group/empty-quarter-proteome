#!/bin/bash
#SBATCH --job-name=exp3_s1_arr
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --array=1-13
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=02:30:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_arr_%A_%a.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_arr_%A_%a.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export CCTYPER_DB=$EXP3_ROOT/inputs/cctyper_db_v1.8.0
export TMPDIR=$EXP3_ROOT/tmp/s1_arr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity/4.3.2 || true
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 1 CRISPR chunk ${SLURM_ARRAY_TASK_ID} on $(hostname) (Job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID})"
python3 "$EXP3_ROOT/code/stage1_crispr_recovery.py" --mode chunk --chunk-id "${SLURM_ARRAY_TASK_ID}" --max-workers 16 --cctyper-threads 8
rm -rf "$TMPDIR"
echo "[$(date -Is)] Completed Stage 1 CRISPR chunk ${SLURM_ARRAY_TASK_ID}."
