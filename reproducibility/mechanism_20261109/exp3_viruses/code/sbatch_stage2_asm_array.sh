#!/bin/bash
#SBATCH --job-name=exp3_s2_arr
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --array=0-7%4
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_arr_%A_%a.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_arr_%A_%a.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s2_arr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity/4.3.2 || true
export PYTHONNOUSERSITE=1

TOTAL_SHARDS=8
SHARD_ID=${SLURM_ARRAY_TASK_ID}

echo "[$(date -Is)] Starting Stage 2 assembly shard ${SHARD_ID}/${TOTAL_SHARDS} on $(hostname) (JOB=${SLURM_ARRAY_JOB_ID}_${SHARD_ID})"
python3 "$EXP3_ROOT/code/stage2_viral_identity_and_spacers_v2.py" \
  --mode asm-shard \
  --shard-id "${SHARD_ID}" \
  --num-shards "${TOTAL_SHARDS}" \
  --max-workers 24 \
  --blast-threads 8 \
  --raw-dir-name stage2_raw_outputs_v2

rm -rf "$TMPDIR"
echo "[$(date -Is)] Completed Stage 2 assembly shard ${SHARD_ID}/${TOTAL_SHARDS}."
