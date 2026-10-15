#!/bin/bash
#SBATCH --job-name=prostt5-smoke-v3
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --gpus=nvidia_h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/smoke_v3_%j_%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/smoke_v3_%j_%r.err

set -euo pipefail

ROOT=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu

# Enforce Slurm worker execution guard
if [[ ${PROSTT5_WORKER_VERIFIED:-0} != 1 ]]; then
  export PROSTT5_WORKER_VERIFIED=1
  exec "$ROOT/require_slurm_worker.sh" /bin/bash "$ROOT/smoke_test_v3.sbatch"
fi

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
PYTHON=/home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python

FINAL="$ROOT/smoke_results_v3"
mkdir -p "$ROOT/logs" "$FINAL"

echo "[$(date -Iseconds)] [Job $SLURM_JOB_ID] Starting ProstT5 smoke v3 test on $(hostname)"

# Monitor GPU details and dynamic utilization
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv > "$FINAL/gpu_info.csv"
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv -l 1 > "$FINAL/gpu_utilization.csv" &
GPU_MONITOR=$!
trap 'kill "$GPU_MONITOR" 2>/dev/null || true' EXIT

# Model cache hash computation as short provenance step on worker
echo "[$(date -Iseconds)] Computing model cache checksums as provenance step on worker..."
MODEL_DIR="$ROOT/model"
if [[ ! -f "$MODEL_DIR/model_cache.sha256" ]]; then
  (cd "$MODEL_DIR" && sha256sum config.json model.safetensors spiece.model tokenizer_config.json special_tokens_map.json added_tokens.json > model_cache.sha256)
fi
cp "$MODEL_DIR/model_cache.sha256" "$FINAL/model_cache.sha256"

# Always invoke runner directly to evaluate comprehensive signature and perform atomic publishing / resume
echo "[$(date -Iseconds)] Executing ProstT5 embedding runner v3 (evaluating comprehensive signature)..."
"$PYTHON" "$ROOT/prostt5_embed.py" \
  --fasta "$ROOT/smoke_references.fasta" \
  --output-dir "$FINAL" \
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
  --output-dir "$FINAL" \
  --fasta "$ROOT/smoke_references.fasta" \
  --report "$FINAL/verification_report.json"

printf 'job=%s node=%s restart=%s date=%s\n' "$SLURM_JOB_ID" "$(hostname -s)" "${SLURM_RESTART_COUNT:-0}" "$(date -Iseconds)" > "$FINAL/slurm_execution.txt"
cp "$ROOT/smoke_test_v3.sbatch" "$FINAL/executed_smoke_v3.sbatch"
cp "$ROOT/require_slurm_worker.sh" "$FINAL/executed_guard.sh"
cp "$ROOT/smoke_references_manifest.tsv" "$FINAL/smoke_references_manifest.tsv"

# Compute bundle checksum
(cd "$FINAL" && sha256sum embeddings.npz embeddings.npy sequence_manifest.tsv provenance.json verification_report.json gpu_info.csv gpu_utilization.csv slurm_execution.txt > final_bundle.sha256)

touch "$FINAL/COMPLETE"

echo "[$(date -Iseconds)] Smoke test v3 job $SLURM_JOB_ID finished successfully."
