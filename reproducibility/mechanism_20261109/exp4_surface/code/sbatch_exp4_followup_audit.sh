#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=370000M
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --job-name=eq_exp4_faudit
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

HISTORICAL_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"
export EQ_EXP4_ROOT="${EQ_EXP4_ROOT:-$HISTORICAL_ROOT}"

CANON_HIST="$(realpath -m "$HISTORICAL_ROOT" 2>/dev/null || readlink -m "$HISTORICAL_ROOT" 2>/dev/null || echo "$HISTORICAL_ROOT")"
CANON_TARGET="$(realpath -m "$EQ_EXP4_ROOT" 2>/dev/null || readlink -m "$EQ_EXP4_ROOT" 2>/dev/null || echo "$EQ_EXP4_ROOT")"

# 1. Refuse historical root execution (normalized path, symlink, trailing slash, or ..)
if [ "$CANON_TARGET" = "$CANON_HIST" ]; then
  echo "[exp4_surface Followup Audit] ERROR: Refusing execution in historical root ($HISTORICAL_ROOT)." >&2
  echo "[exp4_surface Followup Audit] Reanalysis requires an explicit fresh EQ_EXP4_ROOT." >&2
  exit 1
fi

# 2. Verify collector output exists before running audit
COLLECT_EVAL="$EQ_EXP4_ROOT/results/exp4_surface_family_evaluation_v2.tsv"
if [ ! -f "$COLLECT_EVAL" ]; then
  echo "[exp4_surface Followup Audit] ERROR: Required collector evaluation table missing in $EQ_EXP4_ROOT/results ($COLLECT_EVAL)." >&2
  echo "[exp4_surface Followup Audit] Ensure collector stage finishes before launching followup audit." >&2
  exit 1
fi

AUDIT_OUT="$EQ_EXP4_ROOT/results/followup_v2_review/exp4_surface_family_evaluation_v2_audited.tsv"
# 3. Refuse pre-existing audit output
if [ -f "$AUDIT_OUT" ]; then
  echo "[exp4_surface Followup Audit] ERROR: Pre-existing output ($AUDIT_OUT) detected in run root. Refusing execution: reanalysis will not silently overwrite existing results." >&2
  exit 1
fi

export TMPDIR="$EQ_EXP4_ROOT/tmp/followup_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR" "$EQ_EXP4_ROOT/logs" "$EQ_EXP4_ROOT/results" "$EQ_EXP4_ROOT/results/followup_v2_review"
trap 'rm -rf "$TMPDIR"' EXIT

cd "$EQ_EXP4_ROOT"
[ -f /project/k10018/sw/eq-env.sh ] && source /project/k10018/sw/eq-env.sh
command -v conda >/dev/null 2>&1 && conda activate /project/k10018/sw/envs/mapping 2>/dev/null || true
module load singularity >/dev/null 2>&1 || true
export PATH="/sw/ex113genoa/singularity/ce/4.3.2/bin:$PATH"
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C

echo "[exp4_surface Followup Audit] Launching targeted audit on $(hostname) at $(date -Is)..."
python3 "$EQ_EXP4_ROOT/code/exp4_followup_targeted_audit.py"

echo "[exp4_surface Independent Replay] Launching independent replay across all 17 shards on $(hostname) at $(date -Is)..."
python3 "$EQ_EXP4_ROOT/code/exp4_independent_polymer_and_msa_replay.py"

if [ ! -f "$AUDIT_OUT" ]; then
  echo "[exp4_surface Followup Audit] ERROR: Followup audit completed but output missing: $AUDIT_OUT" >&2
  exit 1
fi
echo "[exp4_surface Followup Audit & Replay] Completed successfully at $(date -Is)."
