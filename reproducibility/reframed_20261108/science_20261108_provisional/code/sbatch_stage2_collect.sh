#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=eq_sci_s2_col
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
cd "$EQ_SCIENCE_OUT"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

echo "[Stage 2 Collect] Running unit test gate on compute node $(hostname) at $(date -Is)..."
python3 -m pytest -o cache_dir="$EQ_SCIENCE_OUT/logs/.pytest_cache" "$EQ_SCIENCE_OUT/tests/test_science_20261108.py" -v

SHARD_DIR="$EQ_SCIENCE_OUT/results/stage2_shards"
M8_OUT="$EQ_SCIENCE_OUT/results/thermo538_vs_uniref50_s75_e1e3.m8.gz"

if [ ! -s "$M8_OUT" ]; then
  echo "[Stage 2 Collect] Verifying all 16 shards in $SHARD_DIR at $(date -Is)..."
  for i in $(seq 0 15); do
    if [ ! -f "$SHARD_DIR/shard_${i}.done" ] || [ ! -s "$SHARD_DIR/shard_${i}.m8.gz" ]; then
      echo "ERROR: Missing shard $i in $SHARD_DIR" >&2
      exit 1
    fi
  done
  echo "[Stage 2 Collect] Combining 16 shards into $M8_OUT..."
  zcat "$SHARD_DIR"/shard_{0..15}.m8.gz | gzip -1 > "${M8_OUT}.tmp"
  mv "${M8_OUT}.tmp" "$M8_OUT"
fi

echo "[Stage 2 Collect] Running stage2_uniref50_harmonised.py at $(date -Is)..."
python3 "$EQ_SCIENCE_OUT/code/stage2_uniref50_harmonised.py" "$M8_OUT"
rm -rf "$SHARD_DIR"
echo "[Stage 2 Collect] Done at $(date -Is)."
