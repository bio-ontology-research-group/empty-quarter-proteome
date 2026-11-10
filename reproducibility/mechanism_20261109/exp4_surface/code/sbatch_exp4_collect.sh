#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=eq_exp4_col
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

HISTORICAL_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"
export EQ_EXP4_ROOT="${EQ_EXP4_ROOT:-$HISTORICAL_ROOT}"

CANON_HIST="$(realpath -m "$HISTORICAL_ROOT" 2>/dev/null || readlink -m "$HISTORICAL_ROOT" 2>/dev/null || echo "$HISTORICAL_ROOT")"
CANON_TARGET="$(realpath -m "$EQ_EXP4_ROOT" 2>/dev/null || readlink -m "$EQ_EXP4_ROOT" 2>/dev/null || echo "$EQ_EXP4_ROOT")"

# 1. Refuse historical root execution (normalized path, symlink, trailing slash, or ..)
if [ "$CANON_TARGET" = "$CANON_HIST" ]; then
  echo "[exp4_surface Collect] ERROR: Refusing execution in historical root ($HISTORICAL_ROOT)." >&2
  echo "[exp4_surface Collect] Reanalysis requires an explicit fresh EQ_EXP4_ROOT." >&2
  exit 1
fi

COLLECT_OUT="$EQ_EXP4_ROOT/results/exp4_surface_family_evaluation_v2.tsv"
# 2. Refuse pre-existing collector output
if [ -f "$COLLECT_OUT" ]; then
  echo "[exp4_surface Collect] ERROR: Pre-existing output ($COLLECT_OUT) detected in run root. Refusing execution: reanalysis will not silently overwrite existing results." >&2
  exit 1
fi

# 3. Verify all 17 shard evaluations exist
for s_id in $(seq 0 16); do
  SHARD_EVAL="$EQ_EXP4_ROOT/results/shards/shard_${s_id}_family_evaluation.tsv"
  if [ ! -f "$SHARD_EVAL" ]; then
    echo "[exp4_surface Collect] ERROR: Missing evaluation output for shard $s_id ($SHARD_EVAL)." >&2
    echo "[exp4_surface Collect] Array stage must complete for all 17 shards before running collect." >&2
    exit 1
  fi
done

export TMPDIR="$EQ_EXP4_ROOT/tmp/collect_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR" "$EQ_EXP4_ROOT/logs" "$EQ_EXP4_ROOT/results"
trap 'rm -rf "$TMPDIR"' EXIT

cd "$EQ_EXP4_ROOT"
[ -f /project/k10018/sw/eq-env.sh ] && source /project/k10018/sw/eq-env.sh
command -v conda >/dev/null 2>&1 && conda activate /project/k10018/sw/envs/mapping 2>/dev/null || true

echo "[exp4_surface Collect] Launching collect on $(hostname) at $(date -Is)..."
python3 "$EQ_EXP4_ROOT/code/exp4_surface_pipeline.py" --mode collect

if [ ! -f "$COLLECT_OUT" ]; then
  echo "[exp4_surface Collect] ERROR: Collector completed but output missing: $COLLECT_OUT" >&2
  exit 1
fi
echo "[exp4_surface Collect] Completed successfully at $(date -Is)."
