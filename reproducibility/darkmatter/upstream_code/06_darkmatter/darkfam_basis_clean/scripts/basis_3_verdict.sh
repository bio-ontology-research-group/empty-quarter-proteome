#!/bin/bash
# darkfam_basis_clean step 3 — verdict (local, seconds; run after basis_2 chunks land).
# Family verdict needs only chunk_1; gap verdict needs chunks 2-9.
#   SURVIVING_FAMILIES.txt — family reps with NO UniRef50 hit under the locked rule
#   DROPPED_FAMILIES.txt   — family reps with a hit (not dark under the locked rule)
#   NEWDARK_GAP.txt        — never-screened >=3-member reps that are strict-dark
#                            (candidate NEW families; have no June tables — decision item)
# Usage: bash scripts/basis_3_verdict.sh
set -euo pipefail
export LC_ALL=C
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
R=$B/results

if [ -s "$R/hits_chunk_1.tsv" ]; then
    cut -f1 "$R/hits_chunk_1.tsv" | sort -u > "$R/_fam_hit.txt"
    sort "$R/family_rep_ids.txt" > "$R/_fam_all.txt"
    comm -23 "$R/_fam_all.txt" "$R/_fam_hit.txt" > "$R/SURVIVING_FAMILIES.txt"
    comm -12 "$R/_fam_all.txt" "$R/_fam_hit.txt" > "$R/DROPPED_FAMILIES.txt"
    echo "FAMILIES: total $(wc -l < "$R/_fam_all.txt")  surviving $(wc -l < "$R/SURVIVING_FAMILIES.txt")  dropped $(wc -l < "$R/DROPPED_FAMILIES.txt")"
else
    echo "FAMILIES: chunk_1 not done yet"
fi

done_gap=0
for i in 2 3 4 5 6 7 8 9; do [ -s "$R/hits_chunk_${i}.tsv" ] && done_gap=$((done_gap+1)); done
if [ "$done_gap" -eq 8 ]; then
    cat "$R"/hits_chunk_{2,3,4,5,6,7,8,9}.tsv | cut -f1 | sort -u > "$R/_gap_hit.txt"
    sort "$R/gap_rep_ids.txt" > "$R/_gap_all.txt"
    comm -23 "$R/_gap_all.txt" "$R/_gap_hit.txt" > "$R/NEWDARK_GAP.txt"
    echo "GAP: total $(wc -l < "$R/_gap_all.txt")  strict-dark (candidate new families) $(wc -l < "$R/NEWDARK_GAP.txt")"
else
    echo "GAP: $done_gap/8 chunks done"
fi
