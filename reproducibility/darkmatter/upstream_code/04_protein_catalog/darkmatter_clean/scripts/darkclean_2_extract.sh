#!/bin/bash
# darkclean_2 — extract the uncovered rep sequences from all_reps90 (input to the UniRef50 search).
#
# Submit (after darkclean_1): sbatch .../darkmatter_clean/scripts/darkclean_2_extract.sh
#
#SBATCH --job-name=velazqal_darkclean2
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean2_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean2_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate pyrodigal   # has pigz

M=/ibex/scratch/projects/c2014/rob/metagenomics
CAT=$M/04_protein_catalog
OUT=$CAT/darkmatter_clean
IDS=$OUT/01_uncovered_rep_ids.txt
REPS=$CAT/all_linclust90/all_reps90.faa.gz
OUTFA=$OUT/02_uncovered_reps.faa.gz

[ -s "$IDS" ] || { echo "ERROR: $IDS missing/empty — run darkclean_1 first"; exit 1; }

echo "[$(date -Iseconds)] extracting $(wc -l < "$IDS") uncovered reps from all_reps90"
pigz -dc "$REPS" | awk -v ids="$IDS" '
BEGIN { while ((getline l < ids) > 0) keep[l]=1; close(ids) }
/^>/  { id=substr($1,2); pr=(id in keep) }
      { if (pr) print }
' | pigz -p 8 > "${OUTFA}.tmp"
mv "${OUTFA}.tmp" "$OUTFA"

N_OUT=$(zcat "$OUTFA" | grep -c '^>')
echo "[$(date -Iseconds)] extracted $N_OUT / target $(wc -l < "$IDS")"
[ "$N_OUT" -eq "$(wc -l < "$IDS")" ] || { echo "WARN: count mismatch"; }
echo "DONE_DARKCLEAN2"
