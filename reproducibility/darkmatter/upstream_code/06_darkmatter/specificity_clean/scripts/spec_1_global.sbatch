#!/bin/bash
# §3.1 specificity recompute on the REBUILT 106,042-family universe (June 23,382 superseded).
# Adapted verbatim in rule from scripts/ibex/131_global_specificity.sh + 98b_mgnify_search.sh:
#   diamond blastp --sensitive --evalue 1e-3 -k 1 ; a query rep with >=1 hit occurs in that catalogue,
#   no hit => absent from it under the stated rule.
# CHANGE vs June: query is the 106,042 dark family reps PLUS the 106,042 length/size-matched ANNOTATED
#   control families (darkfam_basis_clean CONTROL_FAMILIES.txt) in one pass -> one DB scan serves both,
#   and every catalogue gets a known-protein floor. Databases are reused untouched (built June 2026).
#   sbatch spec_1_global.sbatch
#SBATCH --job-name=ale_spec_global
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_global_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_global_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DBROOT=$WD/databases
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
QUERY=$OUT/query_all.faa
source $WD/../miniconda3/etc/profile.d/conda.sh; conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
log(){ echo "[$(date -Iseconds)] $*"; }
log "diamond: $(diamond --version 2>&1 | head -1)"
log "query (dark + matched ctrl): $(grep -c '^>' "$QUERY")"

SUMM="$OUT/specificity_summary.tsv"
[ -f "$SUMM" ] || echo -e "catalog\tdb_seqs\tset\tn_with_hit\tn_absent\tn_total\tpct_absent" > "$SUMM"

score(){ # name  db_seqs  — split the hit table into dark vs ctrl and write both rows
  local name=$1 dbn=$2 m8="$OUT/${name}_hits.m8"
  cut -f1 "$m8" | sort -u > "$OUT/${name}_hit.ids"
  grep -v -P "^${name}\t" "$SUMM" > "$SUMM.tmp" && mv "$SUMM.tmp" "$SUMM"
  for set in dark ctrl; do
    local list; [ "$set" = dark ] && list="$B/FAMILIES_FINAL.txt" || list="$B/CONTROL_FAMILIES.txt"
    local tot hit abs pct
    tot=$(wc -l < "$list")
    hit=$(sort -u "$list" | comm -12 - "$OUT/${name}_hit.ids" | wc -l)
    abs=$((tot - hit)); pct=$(awk -v a=$abs -v n=$tot 'BEGIN{printf "%.2f",100*a/n}')
    log "[$name/$set] with-hit=$hit absent=$abs (${pct}%)"
    echo -e "${name}\t${dbn}\t${set}\t${hit}\t${abs}\t${tot}\t${pct}" >> "$SUMM"
  done
}

search(){ # name  db_basename  db_seqs
  local name=$1 db=$2 dbn=$3
  if [ -s "$OUT/${name}_hits.m8" ]; then log "[$name] hits exist — scoring only"; score "$name" "$dbn"; return; fi
  log "[$name] diamond blastp vs $db (--sensitive -e1e-3 -k1)"
  /usr/bin/time -v diamond blastp -q "$QUERY" -d "$db" --sensitive --evalue 1e-3 -k 1 -p "$T" \
     -f 6 qseqid sseqid pident evalue -o "$OUT/${name}_hits.m8" 2> "$OUT/${name}.log" \
     || { log "[$name] blastp FAILED"; return; }
  score "$name" "$dbn"
}

# depth figures are the DIAMOND makedb counts recorded in the June logs
search mgnify90 "$DBROOT/mgnify90/mgy_clusters" 717738164
search gem      "$DBROOT/gem/gem_db"            122978191

log "DONE_SPEC_GLOBAL"
cat "$SUMM"
