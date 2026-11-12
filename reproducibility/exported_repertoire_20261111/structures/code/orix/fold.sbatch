#!/bin/bash
# ColabFold/AlphaFold2 prediction of one shard of precomputed a3m alignments (GPU only; no MSA server).
# submit: sbatch --array=0-(N-1) --export=ALL,SET=monomer,MODEL=alphafold2_ptm,RANK=plddt fold.sbatch
#SBATCH --job-name=eqdark_fold
#SBATCH --output=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/fold_%x_%A_%a.%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-exported-dark-20261111/logs/fold_%x_%A_%a.%r.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
set -euo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
SHARD=$(printf 'shard_%02d' "$SLURM_ARRAY_TASK_ID")
IN=$D/inputs/$SET/$SHARD
OUT=$D/predictions/$SET/$SHARD
BIN=$HOME/pixi-envs/colabfold/.pixi/envs/default/bin
case "$(hostname -s)" in *login*) echo "refusing to run on a login node"; exit 64;; esac
test -d "$IN"; mkdir -p "$OUT" "$D/gpu_logs"
export PATH=$BIN:$PATH
export XDG_CACHE_HOME=$D/cache MPLCONFIGDIR=$D/cache/mpl
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
mkdir -p "$MPLCONFIGDIR"
echo "[$(date -Is)] job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID node=$(hostname) set=$SET model=$MODEL restart=${SLURM_RESTART_COUNT:-0} n_inputs=$(ls $IN | wc -l)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv,noheader -l 30 \
    > "$D/gpu_logs/${SET}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURM_RESTART_COUNT:-0}.csv" &
MON=$!
trap 'kill $MON 2>/dev/null || true' EXIT
colabfold_batch --model-type "$MODEL" --num-models 5 --num-recycle 3 --num-seeds 1 --random-seed 20261111 \
    --rank "$RANK" --data "$D/colabfold_data" --sort-queries-by length --recompile-padding 10 \
    "$IN" "$OUT"
python "$D/scripts/collect_metrics.py" "$SET" "$SHARD"   # seconds: per-shard metrics table and crops
echo "[$(date -Is)] FOLD_SHARD_DONE $SET $SHARD done=$(ls $OUT/*.done.txt 2>/dev/null | wc -l)"
touch "$OUT/SHARD_COMPLETE"
