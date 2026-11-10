#!/bin/bash
#SBATCH --job-name=exp3_s2_col
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_col_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_col_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s2_col_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity/4.3.2 || true
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Running authoritative standalone audit on allocated node..."
python3 "$EXP3_ROOT/code/audit_stage2_receipts_with_corrected_validator.py" \
  --raw-dir "$EXP3_ROOT/results/stage2_raw_outputs_v2" \
  --max-workers 24

echo "[$(date -Is)] Starting Stage 2 collect & annotation on $(hostname) (JOB=${SLURM_JOB_ID})"
python3 "$EXP3_ROOT/code/stage2_viral_identity_and_spacers_v3.py" \
  --mode collect \
  --max-workers 24 \
  --blast-threads 8 \
  --raw-dir-name stage2_raw_outputs_v2

rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE2_DONE"
echo "[$(date -Is)] Completed Stage 2 collect & annotation."
