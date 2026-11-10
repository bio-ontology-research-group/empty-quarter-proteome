#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --job-name=eq_sci_stg3
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
export TMPDIR="/dev/shm/${USER}_eq_sci_stg3_${SLURM_JOB_ID:-$$}"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR"
trap 'rm -rf "$TMPDIR"' EXIT

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity >/dev/null 2>&1 || true

echo "[Stage 3] Verifying singularity binary and container SIFs on $(hostname) at $(date -Is)..."
command -v singularity >/dev/null 2>&1 || { echo "FATAL: singularity not found in PATH" >&2; exit 1; }
test -r /project/k10018/sw/images/eq-annot/mafft.sif || { echo "FATAL: mafft.sif not readable" >&2; exit 1; }
test -r /project/k10018/sw/images/eq-annot/iqtree.sif || { echo "FATAL: iqtree.sif not readable" >&2; exit 1; }

echo "[Stage 3] Pre-flight unit test gate on $(hostname) at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_SCIENCE_OUT/tests/test_stage3_candidate_functional.py"

echo "[Stage 3] Launching stage3_candidate_functional.py on $(hostname) at $(date -Is)..."
python3 "$EQ_SCIENCE_OUT/code/stage3_candidate_functional.py"

echo "[Stage 3] Post-run real-output verification gate at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_SCIENCE_OUT/tests/test_stage3_candidate_functional.py"
echo "[Stage 3] Done at $(date -Is)."
