#!/bin/bash
#SBATCH --job-name=prostt5-smoke
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --gpus=nvidia_h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/smoke_%j_%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/smoke_%j_%r.err

set -euo pipefail

ROOT=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu

# Enforce Slurm worker execution guard
if [[ ${PROSTT5_WORKER_VERIFIED:-0} != 1 ]]; then
  export PROSTT5_WORKER_VERIFIED=1
  exec "$ROOT/require_slurm_worker.sh" /bin/bash "$ROOT/smoke_test.sbatch"
fi

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
PYTHON=/home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python

mkdir -p "$ROOT/logs" "$ROOT/smoke_results"

FINAL="$ROOT/smoke_results"
if [[ -f "$FINAL/COMPLETE" && ${FORCE_RERUN:-0} != 1 ]]; then
  echo "[$(date -Iseconds)] [Job $SLURM_JOB_ID] Final results already exist in $FINAL. Verifying checksums..."
  (cd "$FINAL" && sha256sum -c checksums.sha256 --quiet)
  echo "[$(date -Iseconds)] [Job $SLURM_JOB_ID] Verified up-to-date. Skipping recomputation."
  exit 0
fi

WORK="$ROOT/smoke_results/.work_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$WORK"

test ! -e "$WORK/STARTED"
touch "$WORK/STARTED"

echo "[$(date -Iseconds)] [Job $SLURM_JOB_ID] Starting ProstT5 smoke test on $(hostname)"

# Monitor GPU details and dynamic utilization
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv > "$WORK/gpu_info.csv"
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv -l 1 > "$WORK/gpu_utilization.csv" &
GPU_MONITOR=$!
trap 'kill "$GPU_MONITOR" 2>/dev/null || true' EXIT

# Model cache hash computation as short provenance step on worker
echo "[$(date -Iseconds)] Computing model cache checksums as provenance step on worker..."
MODEL_DIR="$ROOT/model"
if [[ ! -f "$MODEL_DIR/model_cache.sha256" ]]; then
  (cd "$MODEL_DIR" && sha256sum config.json model.safetensors spiece.model tokenizer_config.json > model_cache.sha256)
fi
cp "$MODEL_DIR/model_cache.sha256" "$WORK/model_cache.sha256"

# Run embedding with single-vs-batched consistency verification
echo "[$(date -Iseconds)] Executing ProstT5 embedding runner..."
"$PYTHON" "$ROOT/prostt5_embed.py" \
  --fasta "$ROOT/smoke_references.fasta" \
  --output-dir "$WORK" \
  --model-dir "$MODEL_DIR" \
  --precision half \
  --batch-tokens 4096 \
  --max-seq-len 2048 \
  --verify-consistency

# Stop GPU monitoring
kill "$GPU_MONITOR" 2>/dev/null || true
wait "$GPU_MONITOR" 2>/dev/null || true

# Run independent verification and emit report
echo "[$(date -Iseconds)] Running independent audit and verification..."
"$PYTHON" "$ROOT/verify_embeddings.py" \
  --output-dir "$WORK" \
  --fasta "$ROOT/smoke_references.fasta" \
  --report "$WORK/verification_report.json"

printf 'job=%s node=%s restart=%s date=%s\n' "$SLURM_JOB_ID" "$(hostname -s)" "${SLURM_RESTART_COUNT:-0}" "$(date -Iseconds)" > "$WORK/slurm_execution.txt"
cp "$ROOT/smoke_test.sbatch" "$WORK/executed_smoke.sbatch"
cp "$ROOT/require_slurm_worker.sh" "$WORK/executed_guard.sh"
cp "$ROOT/smoke_references_manifest.tsv" "$WORK/smoke_references_manifest.tsv"

# Compute bundle checksum
(cd "$WORK" && sha256sum embeddings.npz embeddings.npy sequence_manifest.tsv provenance.json verification_report.json gpu_info.csv gpu_utilization.csv slurm_execution.txt > final_bundle.sha256)

touch "$WORK/COMPLETE"

# Atomic sync to FINAL
echo "[$(date -Iseconds)] Publishing final results to $FINAL..."
for f in "$WORK"/*; do
  fname=$(basename "$f")
  cp -a "$f" "$FINAL/$fname"
done

echo "[$(date -Iseconds)] Smoke test job $SLURM_JOB_ID finished successfully."
