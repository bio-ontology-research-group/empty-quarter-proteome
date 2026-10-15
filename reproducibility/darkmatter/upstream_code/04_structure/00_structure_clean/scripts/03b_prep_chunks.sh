#!/bin/bash
# 03b_prep — re-chunk the >=100aa MAG-dark 3Di query into ~15k-seq sub-DBs for the search array.
# Filters to target_ids (>=100aa MAG dark = 1,466,943). Reuses existing 3Di (no re-folding).
# Sources: tier_sub/q_t1,q_t2 (already >=100aa) + 02_prostt5/chunk_001..008 (filter to target here).
#SBATCH --job-name=sc3_prep
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_prep_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_prep_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --partition=batch
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
SC=$M/04_structure/00_structure_clean; O=$SC/03_search
FS=$M/tools/foldseek/bin/foldseek
Q=$O/qchunks; rm -rf "$Q"; mkdir -p "$Q"; : > "$O/manifest.txt"
PER=7000
[ -s "$O/target_ids.txt" ] || { echo "missing target_ids.txt"; exit 2; }

split_one () {   # $1=source DB base   $2=tag   $3=keyfile(all keys to include)
  local base="$1" tag="$2" keyf="$3"
  [ -s "$keyf" ] || { echo "[prep] $tag: no keys, skip"; return 0; }
  split -l "$PER" -d -a 4 "$keyf" "$Q/${tag}_k_"
  for kf in "$Q/${tag}_k_"*; do
    local o="$Q/${tag}_$(basename "$kf" | sed "s/${tag}_k_//")"
    "$FS" createsubdb "$kf" "$base"      "$o"      >/dev/null
    "$FS" createsubdb "$kf" "${base}_ss" "${o}_ss" >/dev/null
    "$FS" createsubdb "$kf" "${base}_h"  "${o}_h"  >/dev/null
    echo "$o" >> "$O/manifest.txt"
  done
  echo "[prep] $tag: $(wc -l < "$keyf") reps -> $(ls "$Q/${tag}_k_"* | wc -l) chunks"
}

# q_t1 / q_t2 are already >=100aa target subsets -> use all their keys (col1 of .lookup)
cut -f1 "$O/tier_sub/q_t1.index" > "$Q/qt1.keys"; split_one "$O/tier_sub/q_t1" qt1 "$Q/qt1.keys"
cut -f1 "$O/tier_sub/q_t2.index" > "$Q/qt2.keys"; split_one "$O/tier_sub/q_t2" qt2 "$Q/qt2.keys"
# folded chunks: filter to target_ids (>=100aa) first, then split
for i in 1 2 3 4 5 6 7 8; do
  c=$SC/02_prostt5/chunk_00${i}/queryDB
  awk 'NR==FNR{t[$1]=1;next}($2 in t){print $1}' "$O/target_ids.txt" "${c}.lookup" > "$Q/c${i}.keys"
  split_one "$c" "c${i}" "$Q/c${i}.keys"
done

NCH=$(wc -l < "$O/manifest.txt")
NSEQ=$(cat "$Q"/*.keys | wc -l)
echo "[prep] total sub-chunks: $NCH   total reps: $NSEQ   (target = 1,466,943)"
echo "DONE_03b_prep"
