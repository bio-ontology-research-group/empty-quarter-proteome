#!/bin/bash
# 02_prostt5 — ProstT5 3Di (createdb ONLY) for the ~1.4M dark reps that lack a 3Di.
# createdb from the proven scripts/06d_prostt5_array.sh; the broad-DB SEARCH (AFDB v6 + PDB + MGnify,
# mirroring the PI's protein_novelty) is the separate 03_search step. createdb needs NO target DB,
# so this runs NOW without AFDB access. Keeps the queryDB (3Di) for 03_search.
#
# Submit after 01_prep: sbatch .../00_structure_clean/scripts/02_prostt5.sh
#
#SBATCH --job-name=sc_prostt5
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc_prostt5_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc_prostt5_%A_%a.err
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=06:00:00
#SBATCH --partition=batch
#SBATCH --exclude=gpu201-23-r

set -euo pipefail
module load rl9-gpustack cuda/12.4.1

FOLDSEEK=/ibex/scratch/projects/c2014/rob/metagenomics/tools/foldseek/bin/foldseek   # official gpu binary (conda/module foldseek crash on GPU)
[ -x "$FOLDSEEK" ] || { echo "ERROR: foldseek-linux-gpu not found"; exit 2; }

STRUCT=/ibex/scratch/projects/c2014/rob/metagenomics/04_structure
PROSTT5_MODEL=$STRUCT/prostt5_weights/prostt5-f16.gguf
AFDB=$STRUCT/afdb_pad/afdb_pad   # GPU-padded AFDB; data symlinks -> /ibex/ai/reference/KSL (needs read access to that reference)

CHUNK=$(printf 'chunk_%03d' "$SLURM_ARRAY_TASK_ID")
CHUNK_FA=$STRUCT/00_structure_clean/01_prep/missing_chunks/$CHUNK.fa
OUTDIR=$STRUCT/00_structure_clean/02_prostt5/$CHUNK
TMPDIR=$STRUCT/00_structure_clean/tmp/prostt5_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUTDIR" "$TMPDIR"

[ -f "$CHUNK_FA" ]      || { echo "[$CHUNK] missing chunk fasta"; exit 2; }
[ -f "$PROSTT5_MODEL" ] || { echo "ProstT5 weights missing at $PROSTT5_MODEL"; exit 2; }
# AFDB not needed here (createdb-only); the broad-DB search is 03_search.

NSEQ=$(grep -c '^>' "$CHUNK_FA")
echo "[$(date -Iseconds)] [$CHUNK] ProstT5 + FoldSeek vs AFDB on $NSEQ proteins"
STAGGER=$(( SLURM_ARRAY_TASK_ID % 4 * 30 + (RANDOM % 15) )); echo "[$CHUNK] stagger ${STAGGER}s"; sleep "$STAGGER"

QUERY_DB=$OUTDIR/queryDB
# Step 1: ProstT5 (AA -> 3Di). Idempotent + 3 retries (GGML cuBLAS init can race).
if [ -f "${QUERY_DB}_ss.dbtype" ]; then
    echo "[$(date -Iseconds)] [$CHUNK] queryDB exists — skip createdb"
else
    for attempt in 1 2 3; do
        if "$FOLDSEEK" createdb "$CHUNK_FA" "$QUERY_DB" --prostt5-model "$PROSTT5_MODEL" --gpu 1 --threads 8; then
            echo "[$(date -Iseconds)] [$CHUNK] createdb ok (attempt $attempt)"; break
        fi
        echo "[$CHUNK] createdb attempt $attempt failed; retry"; rm -rf "$QUERY_DB"*; sleep $((attempt*60))
        [ "$attempt" -eq 3 ] && { echo "[$CHUNK] createdb failed x3"; exit 4; }
    done
fi

# (search vs AFDB v6 + PDB + MGnify = the separate 03_search step, mirroring the PI's protein_novelty.
#  createdb needs no target DB, so this step runs without AFDB access.)
echo "[$(date -Iseconds)] [$CHUNK] 3Di createdb done: $NSEQ proteins. queryDB kept for 03_search."
rm -rf "$TMPDIR"
