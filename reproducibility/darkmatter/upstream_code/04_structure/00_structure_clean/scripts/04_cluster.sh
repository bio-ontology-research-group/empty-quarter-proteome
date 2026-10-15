#!/bin/bash
# 04_cluster — structural clustering of the candidate novel-fold proteins (03c output) into
# fold families. Extracts the candidates' 3Di from the ten source DBs (tier_sub/q_t1, q_t2,
# 02_prostt5/chunk_001..008), concatenates, and clusters with the proven 62b parameters
# (--single-step-clustering --cluster-mode 2 -s 4). Census + admitted-family list at the end.
# ADMISSION RULE (decided 2026-08-15, session 6): a candidate fold family needs >=3 members.
# Each member is already the representative of a distinct 90%-identity sequence cluster, so >=3
# means independently diverged recurrence (consistent with the §1.4 admission logic; no
# sample-span requirement — §1.3 makes a structure claim, not an ecology claim).
# Singletons/pairs carry no claim and are not folded.
# Submit after 03c: sbatch scripts/04_cluster.sh
#SBATCH --job-name=sc4_cluster
#SBATCH --account=pi-hohndor
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc4_cluster_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc4_cluster_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
SC=$M/04_structure/00_structure_clean; O=$SC/03_search; W=$SC/04_cluster
FS=$M/tools/foldseek/bin/foldseek
CAND=$O/structural_candidates.txt
[ -s "$CAND" ] || { echo "ABORT: missing/empty $CAND (run 03c first)"; exit 2; }
NCAND=$(wc -l < "$CAND"); echo "[$(date -Iseconds)] candidates: $NCAND"
rm -rf "$W"; mkdir -p "$W/sub" "$W/tmp" "$W/concat"

# --- extract candidate keys per source DB (lookup: key \t name \t idx; index col1 = keys in DB) ---
subset_one () {  # $1=source DB base  $2=tag
  local base="$1" tag="$2"
  local keyf="$W/sub/keys_${tag}"
  awk -v IDX="${base}.index" 'BEGIN{while((getline l < IDX)>0){split(l,a,"\t");k[a[1]]=1}}
       NR==FNR{c[$1]=1;next} ($1 in k)&&($2 in c){print $1}' "$CAND" "${base}.lookup" > "$keyf"
  local n; n=$(wc -l < "$keyf")
  echo "[$(date -Iseconds)] $tag: $n candidate members"
  [ "$n" -gt 0 ] || return 0
  "$FS" createsubdb "$keyf" "$base"      "$W/sub/$tag"      >/dev/null
  "$FS" createsubdb "$keyf" "${base}_ss" "$W/sub/${tag}_ss" >/dev/null
  "$FS" createsubdb "$keyf" "${base}_h"  "$W/sub/${tag}_h"  >/dev/null
  # createsubdb copies .lookup/.source; our rebuilt lookups (file-idx 0) clash with the
  # multi-file .source and make concatdbs ERROR on the maxSetId check — and neither file
  # is needed downstream (cluster uses ''+_ss, createtsv uses _h). Drop them.
  rm -f "$W/sub/${tag}.lookup" "$W/sub/${tag}.source"
  echo "$W/sub/$tag" >> "$W/dblist.txt"
}
: > "$W/dblist.txt"
subset_one "$O/tier_sub/q_t1" qt1
subset_one "$O/tier_sub/q_t2" qt2
for i in 1 2 3 4 5 6 7 8; do subset_one "$SC/02_prostt5/chunk_00${i}/queryDB" "c${i}"; done

# --- concat (verbatim pattern from scripts/62b_tier1_foldcluster_fast.sh) ---
concat_list () {  # $1=listfile  $2=output base path
  local lst=$1 out=$2 first f nb x g
  first=$(head -1 "$lst")
  for f in ${first}*; do cp "$f" "$out${f#$first}"; done
  tail -n +2 "$lst" | while read -r nb; do
    for x in '' _ss _h; do "$FS" concatdbs "$out$x" "$nb$x" "$out.tmp$x" >/dev/null; done
    for g in "$out".tmp*; do mv -f "$g" "$out${g#$out.tmp}"; done
  done
}
concat_list "$W/dblist.txt" "$W/merged"
NM=$(wc -l < "$W/merged.index")
echo "[$(date -Iseconds)] merged candidate 3Di: $NM (expect $NCAND)"
[ "$NM" -eq "$NCAND" ] || { echo "ABORT: merged count != candidates"; exit 2; }

# --- structural clustering (proven 62b params) ---
"$FS" cluster "$W/merged" "$W/clu" "$W/tmp" --threads "${SLURM_CPUS_PER_TASK:-32}" \
  --split-memory-limit 200G --single-step-clustering --cluster-mode 2 -s 4
"$FS" createtsv "$W/merged" "$W/merged" "$W/clu" "$W/clusters.tsv"

# --- census + admission ---
cut -f1 "$W/clusters.tsv" | sort | uniq -c | sort -rn > "$W/cluster_sizes.txt"
awk '$1>=3{print $2}' "$W/cluster_sizes.txt" > "$W/families_ge3_reps.txt"
NP=$(wc -l < "$W/clusters.tsv"); NCL=$(wc -l < "$W/cluster_sizes.txt")
NF3=$(wc -l < "$W/families_ge3_reps.txt")
NM3=$(awk '$1>=3{s+=$1}END{print s+0}' "$W/cluster_sizes.txt")
N1=$(awk '$1==1' "$W/cluster_sizes.txt" | wc -l); N2=$(awk '$1==2' "$W/cluster_sizes.txt" | wc -l)
{
  echo -e "candidates_clustered\t$NP"
  echo -e "clusters_total\t$NCL"
  echo -e "families_ge3_members_ADMITTED\t$NF3"
  echo -e "members_in_admitted_families\t$NM3"
  echo -e "singleton_clusters_excluded\t$N1"
  echo -e "pair_clusters_excluded\t$N2"
} | tee "$W/SUMMARY.txt"
echo "top 20 cluster sizes:"; head -20 "$W/cluster_sizes.txt"
echo "DONE_04"
