#!/bin/bash
#SBATCH --job-name=exp3_s1_pilot
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_pilot_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s1_pilot_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export CCTYPER_DB=$EXP3_ROOT/inputs/cctyper_db_v1.8.0
export TMPDIR=$EXP3_ROOT/tmp/s1_pilot_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity/4.3.2 || true
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 1 CRISPR pilot (30 archaeal MAGs + tool fixtures) on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 -m pytest -v "$EXP3_ROOT/tests/test_exp3_viruses.py" "$EXP3_ROOT/tests/test_stage1_runtime_repair.py"
python3 "$EXP3_ROOT/code/stage1_crispr_recovery.py" --mode pilot --max-workers 16 --cctyper-threads 8
rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE1_PILOT_DONE"
echo "[$(date -Is)] Completed Stage 1 CRISPR pilot."
