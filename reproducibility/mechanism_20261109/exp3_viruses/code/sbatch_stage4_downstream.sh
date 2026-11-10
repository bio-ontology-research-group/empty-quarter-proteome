#!/bin/bash
#SBATCH --job-name=exp3_s4_mech
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_mech_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_mech_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s4_mech_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 4 downstream mechanism analyses on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 "$EXP3_ROOT/code/stage4_downstream_mechanism.py"
rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE4_DONE"
echo "[$(date -Is)] Completed Stage 4 downstream mechanism analyses."
