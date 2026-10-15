#!/bin/bash
# darkclean_3a — split 02_uncovered_reps.faa.gz into N_CHUNKS chunks for the UniRef50 array.
# Mirrors 27a_split_delta_reps_for_array.sh. Default 120 chunks (~315k reps/chunk for ~37.8M).
#
# Submit (after darkclean_2): sbatch .../darkmatter_clean/scripts/darkclean_3a_chunk.sh
#
#SBATCH --job-name=velazqal_darkclean3a
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3a_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3a_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate pyrodigal

: "${N_CHUNKS:=120}"
M=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$M/04_protein_catalog/darkmatter_clean
INFA=$OUT/02_uncovered_reps.faa.gz
CHUNK_DIR=$OUT/03_uncovered_vs_uniref50/chunks
mkdir -p "$CHUNK_DIR"

[ -s "$INFA" ] || { echo "ERROR: $INFA missing — run darkclean_2 first"; exit 1; }

echo "[$(date -Iseconds)] counting reps"
N_TOTAL=$(zcat "$INFA" | grep -c '^>')
PER_CHUNK=$(( (N_TOTAL + N_CHUNKS - 1) / N_CHUNKS ))
echo "[$(date -Iseconds)] splitting $N_TOTAL reps into $N_CHUNKS chunks (~$PER_CHUNK each)"

zcat "$INFA" | awk -v per="$PER_CHUNK" -v n="$N_CHUNKS" -v dir="$CHUNK_DIR" '
BEGIN { chunk=0; seen=0; fn=sprintf("%s/chunk_%03d.faa", dir, chunk) }
/^>/  { if (seen>0 && seen%per==0 && chunk<n-1) { close(fn); chunk++; fn=sprintf("%s/chunk_%03d.faa", dir, chunk) } seen++ }
      { print > fn }
END   { close(fn); print "wrote "seen" records across "(chunk+1)" chunks" > "/dev/stderr" }'

echo "[$(date -Iseconds)] compressing chunks"
for f in "$CHUNK_DIR"/chunk_*.faa; do pigz -p 2 "$f" & done
wait

echo "$N_CHUNKS" > "$CHUNK_DIR/COUNT"
echo "[$(date -Iseconds)] done — COUNT=$(cat "$CHUNK_DIR/COUNT"), files=$(ls "$CHUNK_DIR"/chunk_*.faa.gz | wc -l)"
echo "DONE_DARKCLEAN3A"
