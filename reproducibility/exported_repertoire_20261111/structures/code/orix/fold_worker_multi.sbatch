#!/bin/bash
# One GPU, NPROC concurrent colabfold_batch queue loops on that GPU (each uses a few GB of GPU memory), so the
# GPU stays busy while another process compiles (CPU) or writes results. Introduced after ORIX flagged the
# one-process-per-GPU layout for low utilisation (51-68 % mean). Each loop claims shards atomically
# (mkdir lock, owner = jobid.slot); colabfold_batch is resume-safe (per-query .done.txt).
# submit: sbatch --partition=pi-hohndor --qos=pi-hohndor --gres=gpu:1 --export=ALL,NPROC=4,SETS="..." fold_worker_multi.sbatch
#SBATCH --job-name=eqdark_wm
#SBATCH --output=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/wm_%j.%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/wm_%j.%r.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=2-00:00:00
set -uo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
BIN=$HOME/pixi-envs/colabfold/.pixi/envs/default/bin
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
export PATH=$BIN:$PATH
export XDG_CACHE_HOME=$D/cache MPLCONFIGDIR=$D/cache/mpl
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
NPROC=${NPROC:-4}
mkdir -p "$MPLCONFIGDIR" "$D/gpu_logs" "$D/locks"
echo "[$(date -Is)] multi-worker job=$SLURM_JOB_ID node=$(hostname) nproc=$NPROC sets=$SETS restart=${SLURM_RESTART_COUNT:-0}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv,noheader -l 30 \
    >> "$D/gpu_logs/wm_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}.csv" &
MON=$!
# Without MPS the processes only time-slice the GPU (measured: 4 processes -> ~17x slower per model, i.e. high
# utilisation but lower throughput). With the CUDA MPS daemon the contexts run concurrently.
if [ "${MPS:-1}" = 1 ] && command -v nvidia-cuda-mps-control >/dev/null; then
  export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/mps_${SLURM_JOB_ID}/pipe CUDA_MPS_LOG_DIRECTORY=$D/logs/mps_${SLURM_JOB_ID}
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  nvidia-cuda-mps-control -d && echo "MPS daemon started" || echo "MPS daemon failed to start"
  sleep 2
fi
trap 'kill $MON 2>/dev/null || true; [ -n "${CUDA_MPS_PIPE_DIRECTORY:-}" ] && echo quit | nvidia-cuda-mps-control 2>/dev/null; rm -rf /dev/shm/mps_${SLURM_JOB_ID}' EXIT
loop() {
  local SLOT=$1 ME="$SLURM_JOB_ID.$1" SET IN SHARD OUT LOCK OWNER MODEL RANK MM MAXMSA
  sleep $((SLOT * 20))
  for SET in $SETS; do
    case "$SET" in monomer*) MODEL=alphafold2_ptm; RANK=plddt; MM=512:1024;; *) MODEL=alphafold2_multimer_v3; RANK=multimer; MM=508:1024;; esac
    case "$SET" in *_shallow) MAXMSA="--max-msa $MM";; *) MAXMSA="";; esac
    for IN in "$D"/inputs/"$SET"/shard_*; do
      SHARD=$(basename "$IN"); OUT=$D/predictions/$SET/$SHARD; LOCK=$D/locks/${SET}_${SHARD}
      [ -e "$OUT/SHARD_COMPLETE" ] && continue
      if ! mkdir "$LOCK" 2>/dev/null; then
        OWNER=$(cat "$LOCK/owner" 2>/dev/null || echo none)
        if [ "$OWNER" != "$ME" ]; then
          case "${OWNER%%.*}" in ''|none|resharded|*[!0-9]*) continue;; esac
          [ "${OWNER%%.*}" != "$SLURM_JOB_ID" ] || continue
          [ -z "$(squeue -h -t RUNNING,PENDING,COMPLETING,SUSPENDED -j "${OWNER%%.*}" -o %i 2>/dev/null)" ] || continue
          echo "[$(date -Is)] slot $SLOT takes over $SET $SHARD from inactive $OWNER"
        fi
      fi
      echo "$ME" > "$LOCK/owner"; sleep 2
      [ "$(cat "$LOCK/owner")" = "$ME" ] || continue
      mkdir -p "$OUT"
      echo "[$(date -Is)] slot $SLOT START $SET $SHARD n_inputs=$(ls "$IN" | wc -l) maxmsa=${MAXMSA:-default}"
      if colabfold_batch --model-type "$MODEL" --num-models 5 --num-recycle 3 --num-seeds 1 --random-seed 20261111 \
          --rank "$RANK" --data "$D/colabfold_data" --sort-queries-by length --recompile-padding 30 $MAXMSA \
          "$IN" "$OUT" > "$D/logs/colabfold_${SET}_${SHARD}_${SLURM_JOB_ID}.${SLURM_RESTART_COUNT:-0}.s$SLOT.log" 2>&1 \
         && python "$D/scripts/collect_metrics.py" "$SET" "$SHARD" >> "$D/logs/collect_${SLURM_JOB_ID}.log" 2>&1; then
        touch "$OUT/SHARD_COMPLETE"
        echo "[$(date -Is)] slot $SLOT FOLD_SHARD_DONE $SET $SHARD done=$(ls "$OUT"/*.done.txt 2>/dev/null | wc -l)"
      else
        echo "[$(date -Is)] slot $SLOT FOLD_SHARD_FAILED $SET $SHARD"; return 1
      fi
    done
  done
  echo "[$(date -Is)] slot $SLOT queue empty"
}
PIDS=()
for s in $(seq 1 "$NPROC"); do loop "$s" & PIDS+=($!); done
if [ -n "${EXTRA_CMD:-}" ]; then      # side task sharing this GPU (e.g. the Foldseek search)
  ( sleep 120; bash "$EXTRA_CMD" > "$D/logs/extra_${SLURM_JOB_ID}.log" 2>&1; echo "[$(date -Is)] EXTRA_CMD exit=$?" ) &
  PIDS+=($!)
fi
rc=0
for p in "${PIDS[@]}"; do wait "$p" || rc=1; done
echo "[$(date -Is)] MULTI_WORKER_DONE rc=$rc"
exit $rc
