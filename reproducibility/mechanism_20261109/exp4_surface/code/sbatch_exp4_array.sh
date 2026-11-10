#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-16
#SBATCH --job-name=eq_exp4_arr
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

HISTORICAL_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"
export EQ_EXP4_ROOT="${EQ_EXP4_ROOT:-$HISTORICAL_ROOT}"
SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

CANON_HIST="$(realpath -m "$HISTORICAL_ROOT" 2>/dev/null || readlink -m "$HISTORICAL_ROOT" 2>/dev/null || echo "$HISTORICAL_ROOT")"
CANON_TARGET="$(realpath -m "$EQ_EXP4_ROOT" 2>/dev/null || readlink -m "$EQ_EXP4_ROOT" 2>/dev/null || echo "$EQ_EXP4_ROOT")"

# 1. Refuse historical root execution to protect frozen evidence (normalized path, symlink, trailing slash, or ..)
if [ "$CANON_TARGET" = "$CANON_HIST" ]; then
  echo "[exp4_surface Array ${SHARD_ID}] ERROR: Refusing execution in historical root ($HISTORICAL_ROOT)." >&2
  echo "[exp4_surface Array ${SHARD_ID}] Historical results are frozen. Reanalysis requires an explicit, fresh, empty EQ_EXP4_ROOT scratch directory." >&2
  exit 1
fi

SHARD_DONE="$EQ_EXP4_ROOT/results/shards/shard_${SHARD_ID}.v2.done"
SHARD_EVAL="$EQ_EXP4_ROOT/results/shards/shard_${SHARD_ID}_family_evaluation.tsv"

# 2. Refuse pre-existing outputs or .done markers; do not silently skip
if [ -f "$SHARD_DONE" ] || [ -f "$SHARD_EVAL" ]; then
  echo "[exp4_surface Array ${SHARD_ID}] ERROR: Pre-existing output ($SHARD_EVAL) or completion marker ($SHARD_DONE) detected in run root." >&2
  echo "[exp4_surface Array ${SHARD_ID}] Fresh-root reanalysis refuses pre-existing outputs and will not silently skip compute on .done markers." >&2
  exit 1
fi

# 3. Verify required input bundle exists
SHARD_PKL="$EQ_EXP4_ROOT/results/shard_inputs/shard_${SHARD_ID}.pkl"
if [ ! -f "$SHARD_PKL" ]; then
  echo "[exp4_surface Array ${SHARD_ID}] ERROR: Required input bundle missing: $SHARD_PKL" >&2
  echo "[exp4_surface Array ${SHARD_ID}] Stage required inputs from immutable original path before submitting shard array." >&2
  exit 1
fi

export TMPDIR="$EQ_EXP4_ROOT/tmp/shard_${SLURM_JOB_ID:-manual}_${SHARD_ID}"
mkdir -p "$TMPDIR" "$EQ_EXP4_ROOT/logs" "$EQ_EXP4_ROOT/results/shards"
trap 'rm -rf "$TMPDIR"' EXIT

cd "$EQ_EXP4_ROOT"
[ -f /project/k10018/sw/eq-env.sh ] && source /project/k10018/sw/eq-env.sh
command -v conda >/dev/null 2>&1 && conda activate /project/k10018/sw/envs/mapping 2>/dev/null || true
module load singularity >/dev/null 2>&1 || true
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C

echo "[exp4_surface Array ${SHARD_ID}] Launching run_shard on $(hostname) at $(date -Is)..."
python3 "$EQ_EXP4_ROOT/code/exp4_surface_pipeline.py" --mode run_shard --shard-id "$SHARD_ID" --threads 192

if [ ! -f "$SHARD_EVAL" ]; then
  echo "[exp4_surface Array ${SHARD_ID}] ERROR: Pipeline completed but output table missing: $SHARD_EVAL" >&2
  exit 1
fi

touch "$SHARD_DONE"
echo "[exp4_surface Array ${SHARD_ID}] Completed successfully at $(date -Is)."
