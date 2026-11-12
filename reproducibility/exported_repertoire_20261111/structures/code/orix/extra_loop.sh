#!/bin/bash
# Extra colabfold queue loop added to an already running multi-worker job (same GPU, same MPS daemon):
#   srun --jobid=JOB --overlap --ntasks=1 bash extra_loop.sh SLOT "SETS"
# Used to raise a three-process job to four/five processes without restarting it.
set -uo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
export PATH=$HOME/pixi-envs/colabfold/.pixi/envs/default/bin:$PATH
export XDG_CACHE_HOME=$D/cache MPLCONFIGDIR=$D/cache/mpl
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
[ -d /dev/shm/mps_${SLURM_JOB_ID}/pipe ] && export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/mps_${SLURM_JOB_ID}/pipe CUDA_MPS_LOG_DIRECTORY=$D/logs/mps_${SLURM_JOB_ID}
SETS="$2"
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
loop "$1"
