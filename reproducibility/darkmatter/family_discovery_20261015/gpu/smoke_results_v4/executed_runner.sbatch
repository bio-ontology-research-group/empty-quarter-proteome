#!/bin/bash
#SBATCH --job-name=prostt5-embed
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --gpus=nvidia_h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/embed_%j_%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/logs/embed_%j_%r.err

#
# Reusable ProstT5 AA2fold GPU Embedding Runner Slurm Script.
#
# Usage:
#   sbatch --export=ALL,INPUT_FASTA=/path/to/proteins.faa,OUTPUT_DIR=/path/to/output run_embed.sbatch
#
# Optional environment overrides:
#   BATCH_TOKENS (default: 4096)
#   MAX_SEQ_LEN (default: 2048)
#   PRECISION (default: half)
#   VERIFY_CONSISTENCY (default: 0)
#

set -euo pipefail

ROOT=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu

# Enforce Slurm worker execution guard
if [[ ${PROSTT5_WORKER_VERIFIED:-0} != 1 ]]; then
  export PROSTT5_WORKER_VERIFIED=1
  exec "$ROOT/require_slurm_worker.sh" /bin/bash "$ROOT/run_embed.sbatch"
fi

if [[ -z "${INPUT_FASTA:-}" || -z "${OUTPUT_DIR:-}" ]]; then
  echo "FATAL: INPUT_FASTA and OUTPUT_DIR must be specified." >&2
  exit 64
fi

test -f "$INPUT_FASTA" || { echo "Input FASTA $INPUT_FASTA not found." >&2; exit 66; }

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
PYTHON=/home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python
MODEL_DIR="$ROOT/model"

BATCH_TOKENS="${BATCH_TOKENS:-4096}"
MAX_SEQ_LEN="${MAX_SEQ_LEN:-2048}"
PRECISION="${PRECISION:-half}"
VERIFY_CONSISTENCY="${VERIFY_CONSISTENCY:-0}"

mkdir -p "$ROOT/logs" "$OUTPUT_DIR"

echo "[$(date -Iseconds)] [Job $SLURM_JOB_ID] Starting ProstT5 embedding workflow for $INPUT_FASTA on $(hostname)"

# Monitor GPU
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv > "$OUTPUT_DIR/gpu_info.csv"
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv -l 2 > "$OUTPUT_DIR/gpu_utilization.csv" &
GPU_MONITOR=$!
trap 'kill "$GPU_MONITOR" 2>/dev/null || true' EXIT

# Model cache hash computation as short provenance step on worker
if [[ ! -f "$MODEL_DIR/model_cache.sha256" ]]; then
  echo "[$(date -Iseconds)] Computing model cache checksums as provenance step on worker..."
  (cd "$MODEL_DIR" && sha256sum config.json model.safetensors spiece.model tokenizer_config.json special_tokens_map.json added_tokens.json > model_cache.sha256)
fi
cp "$MODEL_DIR/model_cache.sha256" "$OUTPUT_DIR/model_cache.sha256"

EXTRA_ARGS=()
if [[ "$VERIFY_CONSISTENCY" == "1" ]]; then
  EXTRA_ARGS+=("--verify-consistency")
fi

# Always invoke runner directly to evaluate comprehensive signature and handle atomic publication / resume
echo "[$(date -Iseconds)] Invoking ProstT5 runner (evaluates comprehensive signature and performs atomic execution)..."
"$PYTHON" "$ROOT/prostt5_embed.py" \
  --fasta "$INPUT_FASTA" \
  --output-dir "$OUTPUT_DIR" \
  --model-dir "$MODEL_DIR" \
  --precision "$PRECISION" \
  --batch-tokens "$BATCH_TOKENS" \
  --max-seq-len "$MAX_SEQ_LEN" \
  "${EXTRA_ARGS[@]}"

kill "$GPU_MONITOR" 2>/dev/null || true
wait "$GPU_MONITOR" 2>/dev/null || true

# Independent verification
echo "[$(date -Iseconds)] Running independent verification..."
"$PYTHON" "$ROOT/verify_embeddings.py" \
  --output-dir "$OUTPUT_DIR" \
  --fasta "$INPUT_FASTA" \
  --report "$OUTPUT_DIR/verification_report.json"

printf 'job=%s node=%s restart=%s input=%s date=%s\n' "$SLURM_JOB_ID" "$(hostname -s)" "${SLURM_RESTART_COUNT:-0}" "$INPUT_FASTA" "$(date -Iseconds)" > "$OUTPUT_DIR/slurm_execution.txt"
cp "$ROOT/run_embed.sbatch" "$OUTPUT_DIR/executed_runner.sbatch"
cp "$ROOT/require_slurm_worker.sh" "$OUTPUT_DIR/executed_guard.sh"

(cd "$OUTPUT_DIR" && sha256sum embeddings.npz embeddings.npy sequence_manifest.tsv provenance.json verification_report.json gpu_info.csv gpu_utilization.csv slurm_execution.txt > final_bundle.sha256)

touch "$OUTPUT_DIR/COMPLETE"

echo "[$(date -Iseconds)] Job $SLURM_JOB_ID completed and verified successfully."
