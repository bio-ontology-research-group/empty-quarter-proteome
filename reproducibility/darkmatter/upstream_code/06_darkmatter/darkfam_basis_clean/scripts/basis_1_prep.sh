#!/bin/bash
# darkfam_basis_clean step 1 — extract query sequences + chunk for the strict re-search.
# Query = 23,382 current dark-family reps (chunk_1) + 667,566 never-screened >=3-member
# orthogroup reps (chunks 2-9, round-robin). IDs prebuilt in results/.
#
# Submit: sbatch scripts/basis_1_prep.sh
#
#SBATCH --job-name=ale_basis_prep
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/prep_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
FAA=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/orthogroups/all_rep_proteins.faa
R=$B/results
NGAP=8   # gap chunks 2..9

for f in "$R/family_rep_ids.txt" "$R/gap_rep_ids.txt"; do
    [ -s "$f" ] || { echo "ERROR: $f missing"; exit 1; }
done
mkdir -p "$R/chunks"

echo "[$(date -Iseconds)] chunk_1 = family reps ($(wc -l < "$R/family_rep_ids.txt"))"
awk 'NR==FNR{ids[$1];next} /^>/{h=substr($1,2); p=(h in ids)} p' \
    "$R/family_rep_ids.txt" "$FAA" > "$R/chunks/chunk_1.faa"
n1=$(grep -c '^>' "$R/chunks/chunk_1.faa")
echo "[$(date -Iseconds)] chunk_1 seqs: $n1 (expect 23382)"
[ "$n1" -eq 23382 ] || { echo "ERROR: chunk_1 count mismatch"; exit 1; }

echo "[$(date -Iseconds)] gap chunks 2..$((NGAP+1)) ($(wc -l < "$R/gap_rep_ids.txt") reps, round-robin)"
awk -v R="$R" -v NG=$NGAP 'NR==FNR{ids[$1];next}
    /^>/{h=substr($1,2); p=(h in ids); if(p){n++; f=sprintf("%s/chunks/chunk_%d.faa", R, ((n-1)%NG)+2)}}
    p{print > f}' "$R/gap_rep_ids.txt" "$FAA"
tot=0
for i in $(seq 2 $((NGAP+1))); do
    c=$(grep -c '^>' "$R/chunks/chunk_${i}.faa"); tot=$((tot+c)); echo "  chunk_$i: $c"
done
echo "[$(date -Iseconds)] gap total: $tot (expect $(wc -l < "$R/gap_rep_ids.txt"))"
[ "$tot" -eq "$(wc -l < "$R/gap_rep_ids.txt")" ] || { echo "ERROR: gap count mismatch"; exit 1; }
touch "$R/CHUNKS_DONE"
echo "[$(date -Iseconds)] done"
