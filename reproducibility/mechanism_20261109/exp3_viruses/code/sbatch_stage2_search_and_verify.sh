#!/bin/bash
#SBATCH --job-name=exp3_s2_vir
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_verify_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s2_verify_%j.err
set -euo pipefail

export EXP3_ROOT=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
export TMPDIR=$EXP3_ROOT/tmp/s2_verify_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity/4.3.2 || true
export PYTHONNOUSERSITE=1

echo "[$(date -Is)] Starting Stage 2 viral identity, raw assembly search, fresh CheckV/geNomad, and spacer scoring on $(hostname) (SLURM_JOB_ID=${SLURM_JOB_ID})"
python3 -m pytest -v "$EXP3_ROOT/tests/test_stage2_to_4_exp3.py"
python3 "$EXP3_ROOT/code/stage2_viral_identity_and_spacers_v2.py" --mode all --max-workers 24 --blast-threads 8 --raw-dir-name stage2_raw_outputs_v2
rm -rf "$TMPDIR"
touch "$EXP3_ROOT/results/STAGE2_DONE"
echo "[$(date -Is)] Completed Stage 2 viral identity and spacer scoring."
