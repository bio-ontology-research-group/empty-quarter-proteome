#!/bin/bash
#SBATCH --job-name=exp3_s3_ledg
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s3_ledger_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s3_ledger_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s3_ledger_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 3 591-archaeal-MAG host-link ledger and historical C1 audit on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 "$EXP3_ROOT/code/stage3_host_link_ledger.py"
rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE3_DONE"
echo "[$(date -Is)] Completed Stage 3 host-link ledger."
