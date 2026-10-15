#!/bin/bash
# 03c_candidates — combine per-chunk hits, restrict to the >=100 aa MAG-dark targets, report the
# retrieval decomposition (AFDB-any / PDB-any / dual / single-reference / neither) and the reps with
# NO structural match in AFDB v6 or PDB = candidate novel-fold proteins (feed 04_cluster next).
# Hard gate: refuses to run unless every manifest chunk has banked a non-empty hit_ids.txt
# (03b writes hit_ids.txt only after BOTH search legs finish, so presence == complete chunk).
# Submit after all 03b array tasks finish.
#SBATCH --job-name=sc3_cand
#SBATCH --account=pi-hohndor
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_cand_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_cand_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$M/04_structure/00_structure_clean/03_search
SORT=(sort -S 48G --parallel "${SLURM_CPUS_PER_TASK:-8}" -T "$OUT/tmp")
mkdir -p "$OUT/tmp"

# --- gate: all chunks banked ---
NEXP=$(wc -l < "$OUT/manifest.txt"); NOK=0; MISS=""
while read -r q; do
  tag=$(basename "$q")
  if [ -s "$OUT/hits/$tag/hit_ids.txt" ]; then NOK=$((NOK+1)); else MISS="$MISS $tag"; fi
done < "$OUT/manifest.txt"
[ "$NOK" -eq "$NEXP" ] || { echo "ABORT: only $NOK/$NEXP chunks banked; missing:$MISS"; exit 2; }
echo "[$(date -Iseconds)] gate passed: $NOK/$NEXP chunks banked"

# --- unions (query IDs with >=1 hit, per reference and overall) ---
"${SORT[@]}" -u "$OUT/target_ids.txt" > "$OUT/target_sorted.txt"
cut -f1 "$OUT"/hits/*/afdb_hits.tsv | "${SORT[@]}" -u > "$OUT/afdb_any_raw.txt"
cut -f1 "$OUT"/hits/*/pdb_hits.tsv  | "${SORT[@]}" -u > "$OUT/pdb_any_raw.txt"
# restrict to targets (defensive; chunk queries are built from target_ids)
comm -12 "$OUT/target_sorted.txt" "$OUT/afdb_any_raw.txt" > "$OUT/afdb_hit_ids.txt"
comm -12 "$OUT/target_sorted.txt" "$OUT/pdb_any_raw.txt"  > "$OUT/pdb_hit_ids.txt"
rm -f "$OUT/afdb_any_raw.txt" "$OUT/pdb_any_raw.txt"
"${SORT[@]}" -u -m "$OUT/afdb_hit_ids.txt" "$OUT/pdb_hit_ids.txt" > "$OUT/all_hit_ids.txt"
comm -12 "$OUT/afdb_hit_ids.txt" "$OUT/pdb_hit_ids.txt" > "$OUT/dual_hit_ids.txt"

comm -12 "$OUT/target_sorted.txt" "$OUT/all_hit_ids.txt" > "$OUT/target_hits.txt"
comm -23 "$OUT/target_sorted.txt" "$OUT/all_hit_ids.txt" > "$OUT/structural_candidates.txt"

NT=$(wc -l < "$OUT/target_sorted.txt")
NA=$(wc -l < "$OUT/afdb_hit_ids.txt"); NP=$(wc -l < "$OUT/pdb_hit_ids.txt")
ND=$(wc -l < "$OUT/dual_hit_ids.txt")
NH=$(wc -l < "$OUT/target_hits.txt"); NC=$(wc -l < "$OUT/structural_candidates.txt")
{
  echo -e "targets_ge100_mag_dark\t$NT"
  echo -e "afdb_any_hit\t$NA"
  echo -e "pdb_any_hit\t$NP"
  echo -e "dual_hit_afdb_and_pdb\t$ND"
  echo -e "afdb_only_hit\t$((NA-ND))"
  echo -e "pdb_only_hit\t$((NP-ND))"
  echo -e "with_structural_hit_any\t$NH"
  echo -e "candidates_no_structural_hit\t$NC"
  awk -v c="$NC" -v t="$NT" 'BEGIN{printf "candidate_fraction\t%.4f\n", (t? c/t : 0)}'
} | tee "$OUT/SUMMARY.txt"
[ $((NH+NC)) -eq "$NT" ] || { echo "ABORT: hit+candidate != targets"; exit 2; }
echo "DONE_03c"
