#!/bin/bash
#SBATCH --job-name=exp3_s0_prep
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s0_prep_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s0_prep_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results" "$EXP3_ROOT/provenance" "$EXP3_ROOT/inputs"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 0 prep and verification on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 -m pytest -v "$EXP3_ROOT/tests/test_exp3_viruses.py"
python3 "$EXP3_ROOT/code/stage0_prep_and_verify.py"
touch "$EXP3_ROOT/results/STAGE0_DONE"
echo "[$(date -Is)] Completed Stage 0 prep and verification."
