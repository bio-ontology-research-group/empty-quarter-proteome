#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000M
#SBATCH --time=00:30:00
#SBATCH --job-name=eq_exp4_stage
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

HISTORICAL_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"
export EQ_EXP4_ROOT="${EQ_EXP4_ROOT:-$HISTORICAL_ROOT}"

CANON_HIST="$(realpath -m "$HISTORICAL_ROOT" 2>/dev/null || readlink -m "$HISTORICAL_ROOT" 2>/dev/null || echo "$HISTORICAL_ROOT")"
CANON_TARGET="$(realpath -m "$EQ_EXP4_ROOT" 2>/dev/null || readlink -m "$EQ_EXP4_ROOT" 2>/dev/null || echo "$EQ_EXP4_ROOT")"

# Refuse historical root execution (normalized path, symlink, trailing slash, or ..)
if [ "$CANON_TARGET" = "$CANON_HIST" ]; then
  echo "[exp4_surface Stage] ERROR: Refusing staging into historical root ($HISTORICAL_ROOT)." >&2
  echo "[exp4_surface Stage] Historical results are frozen. Staging requires an explicit, fresh, empty EQ_EXP4_ROOT scratch directory." >&2
  exit 1
fi

mkdir -p "$EQ_EXP4_ROOT"/{code,inputs/pdb_chains,results/shard_inputs,logs,provenance}

echo "[exp4_surface Stage] Staging immutable input assets from $HISTORICAL_ROOT on $(hostname) at $(date -Is)..."
cp "$HISTORICAL_ROOT/inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz" "$EQ_EXP4_ROOT/inputs/"
cp "$HISTORICAL_ROOT/inputs/curated_role_reference_panel."* "$EQ_EXP4_ROOT/inputs/"
cp "$HISTORICAL_ROOT/inputs/pdb_chains/"*.pdb "$EQ_EXP4_ROOT/inputs/pdb_chains/"
cp "$HISTORICAL_ROOT/results/shard_inputs/"shard_*.pkl "$EQ_EXP4_ROOT/results/shard_inputs/"
cp "$HISTORICAL_ROOT/provenance/frozen_family_universe.tsv" "$EQ_EXP4_ROOT/provenance/"

echo "[exp4_surface Stage] Verifying 17 staged shard bundles against relative manifest..."
cd "$EQ_EXP4_ROOT"
sha256sum -c "$EQ_EXP4_ROOT/manifests/upstream_shard_bundles.sha256"
echo "[exp4_surface Stage] Input staging and SHA256 preflight completed successfully at $(date -Is)."
