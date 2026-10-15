#!/bin/bash
# §3.1 desert-vs-desert specificity on the REBUILT 106,042 universe + matched annotated control.
# Rule identical to scripts/ibex/131b_desert_specificity.sh (diamond blastp --sensitive -e1e-3 -k 1).
# Databases reused untouched (built June 2026 from SPIRE per-study genecalls):
#   atacama=Uritskiy_2019_soil_Atacama · mojave=JGI_soil_desert_Mojave · arid=JGI_soil_arid
#   desert_autotroph=Bay_2021_autotrophy_desert · antarctic_colddesert=JGI_soil_Antarctica
# Writes to a separate summary file to avoid racing spec_1_global.
#   sbatch spec_2_deserts.sbatch
#SBATCH --job-name=ale_spec_deserts
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_deserts_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_deserts_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DESD=$WD/databases/deserts
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
QUERY=$OUT/query_all.faa
source $WD/../miniconda3/etc/profile.d/conda.sh; conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
log(){ echo "[$(date -Iseconds)] $*"; }
log "query (dark + matched ctrl): $(grep -c '^>' "$QUERY")"

SUMM="$OUT/desert_specificity_summary.tsv"
[ -f "$SUMM" ] || echo -e "catalog\tdb_seqs\tset\tn_with_hit\tn_absent\tn_total\tpct_absent" > "$SUMM"

score(){
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

search(){
  local name=$1 db=$2 dbn=$3
  if [ -s "$OUT/${name}_hits.m8" ]; then log "[$name] hits exist — scoring only"; score "$name" "$dbn"; return; fi
  log "[$name] diamond blastp vs $name"
  diamond blastp -q "$QUERY" -d "$db" --sensitive --evalue 1e-3 -k 1 -p "$T" \
     -f 6 qseqid sseqid pident evalue -o "$OUT/${name}_hits.m8" 2> "$OUT/${name}.log" \
     || { log "[$name] blastp FAILED"; return; }
  score "$name" "$dbn"
}

search atacama              "$DESD/atacama_db"              6193359
search desert_autotroph     "$DESD/desert_autotroph_db"     7200701
search arid                 "$DESD/arid_db"                 30569562
search antarctic_colddesert "$DESD/antarctic_colddesert_db" 49206550
search mojave               "$DESD/mojave_db"               54361006
search deserts_combined     "$DESD/deserts_combined_db"     147531178

log "DONE_SPEC_DESERTS"
cat "$SUMM"
