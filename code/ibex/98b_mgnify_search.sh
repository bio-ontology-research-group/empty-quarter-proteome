#!/bin/bash
# Phase D (heavy) .2 — DIAMOND blastp the 23,382 dark rep proteins vs MGnify90, sensitive, e-value 1e-3.
# A hit => the EQ dark family also occurs in global environmental metagenomes (globally shared); no hit
# => absent even from MGnify (the largest env protein catalogue) => desert/EQ-endemic.
#   sbatch --dependency=afterok:<98a> scripts/ibex/98b_mgnify_search.sh
#SBATCH --job-name=hohndor_mgnify_search
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_search_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_search_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DB=$WD/databases/mgnify90
OUT=$WD/06_darkmatter/darkfam_ecology/endemism
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
[ -f "$DB/MGNIFY_DB_DONE" ] || { echo "run 98a first"; exit 1; }

echo "[$(date -Iseconds)] diamond blastp $(grep -c '^>' "$OUT/query.faa") reps vs MGnify90"
diamond blastp -q "$OUT/query.faa" -d "$DB/mgy_clusters" \
    --sensitive --evalue 1e-3 -k 25 -p "$T" \
    -f 6 qseqid sseqid pident length evalue bitscore \
    -o "$OUT/mgnify_hits.m8" 2> "$OUT/diamond.log"
echo "  hit rows: $(wc -l < "$OUT/mgnify_hits.m8")"
echo "  reps with >=1 hit: $(cut -f1 "$OUT/mgnify_hits.m8" | sort -u | wc -l) / $(grep -c '^>' "$OUT/query.faa")"
date -Iseconds > "$OUT/MGNIFY_SEARCH_DONE"
echo "[$(date -Iseconds)] search done -> $OUT/mgnify_hits.m8"
