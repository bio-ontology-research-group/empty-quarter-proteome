#!/bin/bash
# ColabFold/AlphaFold2 prediction of several shards in one job, one colabfold_batch process per GPU
# (the account may hold at most 16 submitted jobs, so shards are grouped per job).
# submit: sbatch --gres=gpu:nvidia_h200:4 --cpus-per-task=32 --mem=256G \
#   --export=ALL,SET=complex_real,MODEL=alphafold2_multimer_v3,RANK=multimer,SHARDS="0 1 2 3" fold_multi.sbatch
#SBATCH --job-name=eqdark_foldm
#SBATCH --output=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/foldm_%x_%j.%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/foldm_%x_%j.%r.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
set -euo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
BIN=$HOME/pixi-envs/colabfold/.pixi/envs/default/bin
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
export PATH=$BIN:$PATH
export XDG_CACHE_HOME=$D/cache MPLCONFIGDIR=$D/cache/mpl
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
mkdir -p "$MPLCONFIGDIR" "$D/gpu_logs"
read -ra SH <<< "$SHARDS"
NG=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
echo "[$(date -Is)] job=$SLURM_JOB_ID node=$(hostname) set=$SET model=$MODEL shards=$SHARDS gpus=$NG restart=${SLURM_RESTART_COUNT:-0}"
[ "${#SH[@]}" -le "$NG" ] || { echo "more shards than GPUs"; exit 2; }
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
nvidia-smi --query-gpu=timestamp,index,utilization.gpu,memory.used --format=csv,noheader -l 30 \
    > "$D/gpu_logs/${SET}_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}.csv" &
MON=$!
trap 'kill $MON 2>/dev/null || true' EXIT
PIDS=()
g=0
for s in "${SH[@]}"; do
  SHARD=$(printf 'shard_%02d' "$s"); IN=$D/inputs/$SET/$SHARD; OUT=$D/predictions/$SET/$SHARD
  test -d "$IN"; mkdir -p "$OUT"
  ( export CUDA_VISIBLE_DEVICES=$g
    colabfold_batch --model-type "$MODEL" --num-models 5 --num-recycle 3 --num-seeds 1 --random-seed 20261111 \
        --rank "$RANK" --data "$D/colabfold_data" --sort-queries-by length --recompile-padding 10 \
        "$IN" "$OUT" > "$D/logs/foldm_${SET}_${SHARD}_${SLURM_JOB_ID}.${SLURM_RESTART_COUNT:-0}.log" 2>&1
    python "$D/scripts/collect_metrics.py" "$SET" "$SHARD"
    touch "$OUT/SHARD_COMPLETE"
    echo "[$(date -Is)] FOLD_SHARD_DONE $SET $SHARD done=$(ls $OUT/*.done.txt 2>/dev/null | wc -l)" ) &
  PIDS+=($!); g=$((g + 1))
done
rc=0
for p in "${PIDS[@]}"; do wait "$p" || rc=1; done
echo "[$(date -Is)] FOLD_MULTI_DONE rc=$rc"
exit $rc
