#!/bin/bash
# 03a_prep_query — REUSE existing 3Di (no re-folding). Build the >=100 aa MAG-dark target list, then
# subset the tier1/tier2 merged 3Di DBs to those targets (the ~2.9M reuse portion). The ~1.4M newly
# folded reps are already clean DBs in 02_prostt5/chunk_* and are searched as-is by 03c. CPU only.
#SBATCH --job-name=sc3_prep
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_prep_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_prep_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
S=$M/04_structure; SC=$S/00_structure_clean
FS=$M/tools/foldseek/bin/foldseek
PC=$M/04_protein_catalog/darkmatter_clean
OUT=$SC/03_search; mkdir -p "$OUT/tier_sub"

# target ids = MAG-mapped dark  AND  >=100 aa (~1.47M)
zcat "$PC/05_mag_dark/mag_dark_rep_ids.txt.gz" | sort -u > "$OUT/mag_ids.txt"
zcat "$PC/06_dark_lengths/dark_rep_len.tsv.gz" | awk '$2>=100{print $1}' | sort -u > "$OUT/ge100_ids.txt"
comm -12 "$OUT/mag_ids.txt" "$OUT/ge100_ids.txt" > "$OUT/target_ids.txt"
echo "[03a] targets (>=100aa MAG dark): $(wc -l < "$OUT/target_ids.txt")"

# subset a tier's merged 3Di DB (base + _ss + _h + .lookup) to the target reps present in it
subset_tier () {   # $1 = source merged DB base   $2 = out base
  awk 'NR==FNR{t[$1]=1;next} ($2 in t){print $1}' "$OUT/target_ids.txt" "$1.lookup" > "$2.keys"
  for suf in "" _ss _h; do "$FS" createsubdb "$2.keys" "$1$suf" "$2$suf"; done
  awk 'NR==FNR{k[$1]=1;next} ($1 in k)' "$2.keys" "$1.lookup" > "$2.lookup"
  echo "[03a] $(basename "$2"): $(wc -l < "$2.keys") reps subset"
}
subset_tier "$S/tier1_cluster/merged" "$OUT/tier_sub/q_t1"
subset_tier "$S/tier2_cluster/merged" "$OUT/tier_sub/q_t2"
echo "DONE_03a"
