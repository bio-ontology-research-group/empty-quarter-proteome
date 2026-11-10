#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --job-name=eq_can_1_1b_4
#SBATCH --output=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/logs/%x_%j.err

set -euo pipefail
export EQ_ORIG_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
export EQ_CANONICAL_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy
export TMPDIR="/dev/shm/${USER}_eq_can_1_1b_4_${SLURM_JOB_ID:-$$}"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR="$TMPDIR/singularity_cache"
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR" "$EQ_CANONICAL_SCIENCE_OUT/logs" "$EQ_CANONICAL_SCIENCE_OUT/results"
trap 'rm -rf "$TMPDIR"' EXIT

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

echo "[Canonical 1/1b/4] Pre-flight unit test gate on $(hostname) at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_ORIG_SCIENCE_OUT/tests/test_canonical_taxonomy_sensitivity.py"

echo "[Canonical 1/1b/4] Running canonical_taxonomy_sensitivity.py --mode all_stage1_1b_4 on $(hostname) at $(date -Is)..."
python3 "$EQ_ORIG_SCIENCE_OUT/code/canonical_taxonomy_sensitivity.py" --mode all_stage1_1b_4

echo "[Canonical 1/1b/4] Post-run real-output verification gate at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_ORIG_SCIENCE_OUT/tests/test_canonical_taxonomy_sensitivity.py"
echo "[Canonical 1/1b/4] Done at $(date -Is)."
