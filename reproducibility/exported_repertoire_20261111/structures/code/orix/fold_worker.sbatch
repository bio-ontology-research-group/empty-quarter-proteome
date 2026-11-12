#!/bin/bash
# Queue worker: one GPU, drains unclaimed shards of the listed sets in priority order.
# A shard is claimed atomically (mkdir lock); colabfold_batch is resume-safe, so a requeued worker continues
# its own claimed shard. Any number of identical workers can be submitted on any partition.
# submit: sbatch [partition/qos/gres options] --export=ALL,SETS="monomer_tail complex_real complex_control" fold_worker.sbatch
#SBATCH --job-name=eqdark_worker
#SBATCH --output=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/worker_%j.%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/worker_%j.%r.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
set -uo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
BIN=$HOME/pixi-envs/colabfold/.pixi/envs/default/bin
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
export PATH=$BIN:$PATH
export XDG_CACHE_HOME=$D/cache MPLCONFIGDIR=$D/cache/mpl
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
mkdir -p "$MPLCONFIGDIR" "$D/gpu_logs" "$D/locks"
echo "[$(date -Is)] worker job=$SLURM_JOB_ID node=$(hostname) sets=$SETS restart=${SLURM_RESTART_COUNT:-0}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv,noheader -l 30 \
    >> "$D/gpu_logs/worker_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}.csv" &
MON=$!
trap 'kill $MON 2>/dev/null || true' EXIT
rc=0
for SET in $SETS; do
  case "$SET" in monomer*) MODEL=alphafold2_ptm; RANK=plddt; MM=512:1024;; *) MODEL=alphafold2_multimer_v3; RANK=multimer; MM=508:1024;; esac
  # *_shallow sets hold alignments with <= 1500 rows: reduced MSA padding, no sequence is dropped
  case "$SET" in *_shallow) MAXMSA="--max-msa $MM";; *) MAXMSA="";; esac
  for IN in "$D"/inputs/"$SET"/shard_*; do
    SHARD=$(basename "$IN"); OUT=$D/predictions/$SET/$SHARD; LOCK=$D/locks/${SET}_${SHARD}
    [ -e "$OUT/SHARD_COMPLETE" ] && continue
    if ! mkdir "$LOCK" 2>/dev/null; then
      OWNER=$(cat "$LOCK/owner" 2>/dev/null || echo none)
      if [ "$OWNER" != "$SLURM_JOB_ID" ]; then
        # take over only when the owning job no longer exists in the queue (cancelled or failed worker)
        case "$OWNER" in ''|none|resharded|*[!0-9]*) continue;; esac
        [ -z "$(squeue -h -j "$OWNER" -o %i 2>/dev/null)" ] || continue
        echo "[$(date -Is)] taking over $SET $SHARD from inactive job $OWNER"
      fi
    fi
    echo "$SLURM_JOB_ID" > "$LOCK/owner"
    mkdir -p "$OUT"
    echo "[$(date -Is)] START $SET $SHARD n_inputs=$(ls "$IN" | wc -l) maxmsa=${MAXMSA:-default}"
    if colabfold_batch --model-type "$MODEL" --num-models 5 --num-recycle 3 --num-seeds 1 --random-seed 20261111 \
        --rank "$RANK" --data "$D/colabfold_data" --sort-queries-by length --recompile-padding 10 $MAXMSA \
        "$IN" "$OUT" > "$D/logs/colabfold_${SET}_${SHARD}_${SLURM_JOB_ID}.${SLURM_RESTART_COUNT:-0}.log" 2>&1 \
       && python "$D/scripts/collect_metrics.py" "$SET" "$SHARD"; then
      touch "$OUT/SHARD_COMPLETE"
      echo "[$(date -Is)] FOLD_SHARD_DONE $SET $SHARD done=$(ls "$OUT"/*.done.txt 2>/dev/null | wc -l)"
    else
      # a failure here is usually pre-emption (SIGTERM): stop at once, the requeued job resumes this shard
      echo "[$(date -Is)] FOLD_SHARD_FAILED $SET $SHARD"; exit 1
    fi
  done
done
echo "[$(date -Is)] WORKER_DONE rc=$rc"
exit $rc
