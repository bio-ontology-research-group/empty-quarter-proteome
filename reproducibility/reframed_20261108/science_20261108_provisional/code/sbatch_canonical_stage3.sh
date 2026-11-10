#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --job-name=eq_can_stg3
#SBATCH --output=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/logs/%x_%j.err

set -euo pipefail
export EQ_ORIG_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
export EQ_CANONICAL_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261109-canonical-taxonomy
export EQ_SCIENCE_OUT="$EQ_CANONICAL_SCIENCE_OUT"
export EQ_ARCH_ROOT="$EQ_CANONICAL_SCIENCE_OUT/overlay_arch"
export EQ_SCIENCE_IN_ROOT="$EQ_CANONICAL_SCIENCE_OUT/inputs"
export TMPDIR="/dev/shm/${USER}_eq_can_stg3_${SLURM_JOB_ID:-$$}"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR" "$EQ_CANONICAL_SCIENCE_OUT/logs" "$EQ_CANONICAL_SCIENCE_OUT/results"
trap 'rm -rf "$TMPDIR"' EXIT

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity >/dev/null 2>&1 || true

echo "[Canonical Stage 3] Verifying singularity binary and container SIFs on $(hostname) at $(date -Is)..."
command -v singularity >/dev/null 2>&1 || { echo "FATAL: singularity not found in PATH" >&2; exit 1; }
test -r /project/k10018/sw/images/eq-annot/mafft.sif || { echo "FATAL: mafft.sif not readable" >&2; exit 1; }
test -r /project/k10018/sw/images/eq-annot/iqtree.sif || { echo "FATAL: iqtree.sif not readable" >&2; exit 1; }

echo "[Canonical Stage 3] Pre-flight unit test gate on $(hostname) at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_ORIG_SCIENCE_OUT/tests/test_stage3_candidate_functional.py" "$EQ_ORIG_SCIENCE_OUT/tests/test_canonical_taxonomy_sensitivity.py"

echo "[Canonical Stage 3] Running canonical_taxonomy_sensitivity.py --mode stage3_canonical on $(hostname) at $(date -Is)..."
python3 "$EQ_ORIG_SCIENCE_OUT/code/canonical_taxonomy_sensitivity.py" --mode stage3_canonical

echo "[Canonical Stage 3] Post-run real-output verification gate at $(date -Is)..."
pytest -o cache_dir="$TMPDIR/pytest_cache" -v "$EQ_ORIG_SCIENCE_OUT/tests/test_stage3_candidate_functional.py" "$EQ_ORIG_SCIENCE_OUT/tests/test_canonical_taxonomy_sensitivity.py"
echo "[Canonical Stage 3] Done at $(date -Is)."
