#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-15
#SBATCH --job-name=eq_sci_s2_arr
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/%x_%A_%a.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/%x_%A_%a.err

set -euo pipefail
export EQ_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
cd "$EQ_SCIENCE_OUT"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity >/dev/null 2>&1 || true
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C

SHARD_ID="${SLURM_ARRAY_TASK_ID}"
N_SHARDS=16
SHARD_DIR="$EQ_SCIENCE_OUT/results/stage2_shards"
mkdir -p "$SHARD_DIR"

SHARD_M8="$SHARD_DIR/shard_${SHARD_ID}.m8.gz"
SHARD_DONE="$SHARD_DIR/shard_${SHARD_ID}.done"

if [ -s "$SHARD_M8" ] && [ -f "$SHARD_DONE" ]; then
  echo "[Stage 2 Array ${SHARD_ID}] Shard already completed ($SHARD_M8). Skipping."
  exit 0
fi

T="/dev/shm/eq_stage2_${SLURM_JOB_ID}_${SHARD_ID}"
MM_TMP="$EQ_SCIENCE_OUT/results/mmseqs_tmp_${SLURM_JOB_ID}_${SHARD_ID}"
mkdir -p "$T" "$MM_TMP"
trap 'rm -rf "$T" "$MM_TMP"' EXIT

echo "[Stage 2 Array ${SHARD_ID}/${N_SHARDS}] Extracting shard sequences to $T/shard.faa at $(date -Is)..."
python3 - "$T/shard.faa" "$SHARD_ID" "$N_SHARDS" << 'EOF'
import sys, tarfile
from pathlib import Path

out_faa = sys.argv[1]
shard_id = int(sys.argv[2])
n_shards = int(sys.argv[3])
tar_p = "/project/k10018/work/eq-annot-20261030/results/c7_thermo_evol/thermo_proteomes.tar.gz"

seq_idx = -1
write_cur = False
n_written = 0
with tarfile.open(tar_p, "r:gz") as tf, open(out_faa, "w") as out:
    for m in tf:
        if not m.isfile() or not m.name.endswith(".faa"):
            continue
        gid = Path(m.name).stem
        f = tf.extractfile(m)
        if f is None:
            continue
        for raw in f:
            line = raw.decode("utf-8", errors="replace")
            if line.startswith(">"):
                seq_idx += 1
                write_cur = (seq_idx % n_shards) == shard_id
                if write_cur:
                    tok = line[1:].strip().split()[0]
                    out.write(f">{gid}|{tok}\n")
                    n_written += 1
            else:
                if write_cur:
                    out.write(line)
print(f"Shard {shard_id}/{n_shards}: wrote {n_written} sequences (out of {seq_idx + 1} total)", flush=True)
EOF

echo "[Stage 2 Array ${SHARD_ID}] Running MMseqs2 easy-search (-s 7.5 -e 1e-3 --max-accept 10) vs UniRef50 2026_01 at $(date -Is)..."
singularity exec -B /project,/scratch,/dev/shm /project/k10018/sw/images/eq-annot/mmseqs2.sif \
  mmseqs easy-search "$T/shard.faa" /project/k10018/sw/db/eq-c2/uniref50/uniref50 \
  "$T/shard.m8" "$MM_TMP" \
  -s 7.5 -e 1e-3 --max-accept 10 --threads "${SLURM_CPUS_ON_NODE:-192}" --split-memory-limit 250G \
  --format-output query,target,pident,alnlen,evalue,bits >/dev/null

gzip -1 -c "$T/shard.m8" > "${SHARD_M8}.tmp"
mv "${SHARD_M8}.tmp" "$SHARD_M8"
date -Is > "$SHARD_DONE"
rm -rf "$MM_TMP" "$T"
echo "[Stage 2 Array ${SHARD_ID}] Completed successfully at $(date -Is)."
