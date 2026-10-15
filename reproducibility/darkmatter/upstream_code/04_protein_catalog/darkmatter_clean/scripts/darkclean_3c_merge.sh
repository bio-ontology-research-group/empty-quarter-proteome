#!/bin/bash
# darkclean_3c — merge the uncovered-rep UniRef50 chunks into hit/nohit id lists.
# Mirrors 27c_merge_uniref50_results.sh.
#
# Submit (after darkclean_3b array): sbatch .../darkmatter_clean/scripts/darkclean_3c_merge.sh
#
#SBATCH --job-name=velazqal_darkclean3c
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3c_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3c_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate pyrodigal

M=/ibex/scratch/projects/c2014/rob/metagenomics
OUTDIR=$M/04_protein_catalog/darkmatter_clean/03_uncovered_vs_uniref50
PER_CHUNK=$OUTDIR/per_chunk
IDS=$M/04_protein_catalog/darkmatter_clean/01_uncovered_rep_ids.txt

HITS=$OUTDIR/hits.tsv
HIT_IDS=$OUTDIR/hit_ids.txt
ALL_IDS=$OUTDIR/all_query_ids.txt
NOHIT_IDS=$OUTDIR/nohit_ids.txt

echo "[$(date -Iseconds)] concat per-chunk hits"
cat "$PER_CHUNK"/hits_chunk_*.tsv > "${HITS}.tmp" && mv "${HITS}.tmp" "$HITS"
echo "  $(wc -l < "$HITS") hit rows"

echo "[$(date -Iseconds)] hit ids (unique col1)"
cut -f1 "$HITS" | sort -u > "$HIT_IDS"

echo "[$(date -Iseconds)] all query ids = the uncovered reps we searched"
sort -u "$IDS" > "$ALL_IDS"

echo "[$(date -Iseconds)] nohit = all - hit"
comm -23 "$ALL_IDS" "$HIT_IDS" > "$NOHIT_IDS"

echo "[$(date -Iseconds)] done."
echo "  all:   $(wc -l < "$ALL_IDS")"
echo "  hit:   $(wc -l < "$HIT_IDS")"
echo "  nohit: $(wc -l < "$NOHIT_IDS")   (these are the NEW dark reps from the uncovered set)"
echo "DONE_DARKCLEAN3C"
