#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --job-name=eq_exp4_gspa
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

HISTORICAL_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"
export EQ_EXP4_ROOT="${EQ_EXP4_ROOT:-$HISTORICAL_ROOT}"

CANON_HIST="$(realpath -m "$HISTORICAL_ROOT" 2>/dev/null || readlink -m "$HISTORICAL_ROOT" 2>/dev/null || echo "$HISTORICAL_ROOT")"
CANON_TARGET="$(realpath -m "$EQ_EXP4_ROOT" 2>/dev/null || readlink -m "$EQ_EXP4_ROOT" 2>/dev/null || echo "$EQ_EXP4_ROOT")"

# Refuse historical root execution (normalized path, symlink, trailing slash, or ..)
if [ "$CANON_TARGET" = "$CANON_HIST" ]; then
  echo "[exp4_surface GSPA Eval] ERROR: Refusing execution in historical root ($HISTORICAL_ROOT)." >&2
  echo "[exp4_surface GSPA Eval] Reanalysis requires an explicit fresh EQ_EXP4_ROOT." >&2
  exit 1
fi

GSPA_OUT="$EQ_EXP4_ROOT/results/gspa_evaluation_summary.json"
if [ -f "$GSPA_OUT" ]; then
  echo "[exp4_surface GSPA Eval] ERROR: Pre-existing output ($GSPA_OUT) detected in run root. Refusing execution: reanalysis will not silently overwrite existing results." >&2
  exit 1
fi

export TMPDIR="$EQ_EXP4_ROOT/tmp/gspa_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR" "$EQ_EXP4_ROOT/logs" "$EQ_EXP4_ROOT/results"
trap 'rm -rf "$TMPDIR"' EXIT

cd "$EQ_EXP4_ROOT"
[ -f /project/k10018/sw/eq-env.sh ] && source /project/k10018/sw/eq-env.sh
command -v conda >/dev/null 2>&1 && conda activate /project/k10018/sw/envs/mapping 2>/dev/null || true

echo "[exp4_surface GSPA Eval] Launching gspa_eval on $(hostname) at $(date -Is)..."
python3 "$EQ_EXP4_ROOT/code/exp4_surface_pipeline.py" --mode gspa_eval --threads 64

if [ ! -f "$GSPA_OUT" ]; then
  echo "[exp4_surface GSPA Eval] ERROR: GSPA eval completed but output missing: $GSPA_OUT" >&2
  exit 1
fi
echo "[exp4_surface GSPA Eval] Completed successfully at $(date -Is)."
