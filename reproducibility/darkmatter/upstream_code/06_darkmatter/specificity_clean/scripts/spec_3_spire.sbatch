#!/bin/bash
# §3.1 SPIRE leg — the comparator June never got. Closes the abstract's "pending SPIRE" hedge.
#
# WHAT THIS DB IS (verified on disk 2026-08-17, do not confuse with the desert leg):
#   databases/spire/rep.tar contains spire_representative_genomes/specI_v4_*.fa.gz — SPIRE's
#   DEREPLICATED species-REPRESENTATIVE GENOMES. Those were gene-called locally with prodigal -p meta
#   into spire_all.faa.parts/ (107,008 parts for 107,008 genomes; prodigal COMPLETED in June).
#   => genome-level content, all biomes, ~271M proteins (est. from 400 sampled parts).
#   The five desert DBs are a DIFFERENT slice of the same resource: SPIRE per-STUDY genecalls
#   (SAMN*.psa_megahit.prodigal.faa.gz) = every gene on every assembled contig, binned or not,
#   assembly-level, 5 studies, 147,531,178 proteins. Assembly-level is the MORE inclusive kind.
#
# WHY JUNE FAILED: scripts/ibex/131_global_specificity.sh:65 did `cat "$pd"/*.faa > "$out" 2>/dev/null`
#   — 107,008 args blew ARG_MAX and `2>/dev/null` swallowed the error, leaving spire_all.faa at 0 bytes.
#   Fixed here with find -print0 | xargs -0, and NO error-swallowing redirect.
#
# WHY A FIFO: staging the concatenated FASTA would cost ~118 GB and /ibex/project/c2014 is 99% full
#   (~1.6 TB of 145 TB shared). `diamond makedb --in <fifo>` was TESTED on 3 parts and reported
#   exactly the grep -c truth (12,999), so we stream instead of staging. The index itself (~75 GB est.)
#   is built under /ibex/user/velazqal (1.5 TB, 5% used), NOT on the shared project volume.
#
# Search rule identical to spec_1/spec_2 so the rows are comparable: --sensitive -e 1e-3 -k 1.
# Query is the same 212,084-seq dark+matched-control FASTA -> SPIRE gets a known-protein floor too.
#   sbatch spec_3_spire.sbatch
#SBATCH --job-name=ale_spec_spire
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_spire_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_spire_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
QUERY=$OUT/query_all.faa
PARTS=$WD/databases/spire/spire_all.faa.parts
BUILD=/ibex/user/velazqal/eqmags_31/spire      # own quota — see project_c2014_storage_pressure
DB=$BUILD/spire_rep_db
MKLOG=$BUILD/makedb.log
mkdir -p "$BUILD/tmp" "$S/logs"
source $WD/../miniconda3/etc/profile.d/conda.sh; conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
log(){ echo "[$(date -Iseconds)] $*"; }
log "diamond: $(diamond --version 2>&1 | head -1)"

# ---- guard: never fill a 99%-full volume, and check our own quota has room for the index ----
avail_gb(){ df -B1G --output=avail "$1" | tail -1 | tr -d ' '; }
A=$(avail_gb "$BUILD")
log "free on $BUILD: ${A} GB (index needs ~75 GB)"
[ "$A" -lt 200 ] && { log "ABORT: <200 GB free on the build volume"; exit 1; }

SUMM="$OUT/specificity_summary.tsv"
[ -f "$SUMM" ] || echo -e "catalog\tdb_seqs\tset\tn_with_hit\tn_absent\tn_total\tpct_absent" > "$SUMM"

score(){ # name  db_seqs
  # NOTE: these MUST be separate `local` statements. Bash expands the whole `local a=$1 b="${a}..."`
  # command line before assigning, so ${name} would resolve against the ENCLOSING scope; under
  # `set -u` that is an unbound-variable abort. spec_1/spec_2 survive the one-liner only because they
  # call score() from inside search(), which happens to have its own `local name`. Job 50613399 died
  # here after 3h of successful work — do not re-merge these lines.
  local name=$1
  local dbn=$2
  local m8="$OUT/${name}_hits.m8"
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

# ---------------- build: stream 107,008 parts through a FIFO into makedb ----------------
if [ ! -s "${DB}.dmnd" ]; then
  NP=$(find "$PARTS" -maxdepth 1 -type f -name '*.faa' | wc -l)
  log "parts found: $NP (expect 107008)"
  [ "$NP" -lt 107008 ] && log "WARNING: fewer parts than expected — proceeding, count is recorded"
  PIPE="$BUILD/tmp/spire.pipe"; rm -f "$PIPE"; mkfifo "$PIPE"
  log "makedb via FIFO (no ~118 GB intermediate)"
  find "$PARTS" -maxdepth 1 -type f -name '*.faa' -print0 | xargs -0 cat > "$PIPE" &
  WPID=$!
  diamond makedb --in "$PIPE" -d "$DB" -p "$T" --tmpdir "$BUILD/tmp" 2>&1 | tee "$MKLOG"
  MKRC=${PIPESTATUS[0]}
  wait $WPID; WRC=$?
  rm -f "$PIPE"
  log "makedb rc=$MKRC writer rc=$WRC"
  [ "$MKRC" -ne 0 ] && { log "ABORT: makedb failed"; exit 1; }
else
  log "index exists — skipping build"
fi

# depth figure comes from makedb's own log (a FIFO cannot be pre-counted)
NSEQ=$(grep -oP 'Database sequences\s+\K[0-9]+' "$MKLOG" 2>/dev/null | tail -1)
[ -z "${NSEQ:-}" ] && { log "ABORT: could not read sequence count from $MKLOG"; exit 1; }
log "SPIRE representative-genome proteins: $NSEQ"

# ---------------- search ----------------
if [ -s "$OUT/spire_hits.m8" ]; then
  log "[spire] hits exist — scoring only"
else
  log "[spire] diamond blastp vs $DB (--sensitive -e1e-3 -k1)"
  /usr/bin/time -v diamond blastp -q "$QUERY" -d "$DB" --sensitive --evalue 1e-3 -k 1 -p "$T" \
     -f 6 qseqid sseqid pident evalue -o "$OUT/spire_hits.m8" 2> "$OUT/spire.log" \
     || { log "[spire] blastp FAILED"; exit 1; }
fi
score spire "$NSEQ"

# ---------------- refresh the union rows now that SPIRE exists ----------------
log "recomputing union no-hit sets including SPIRE"
U="$OUT/union_specificity_summary.tsv"
tmp=$(mktemp -d)
sort -u "$B/FAMILIES_FINAL.txt"  > $tmp/dark.s
sort -u "$B/CONTROL_FAMILIES.txt" > $tmp/ctrl.s
sort -u "$OUT"/mgnify90_hit.ids "$OUT"/gem_hit.ids "$OUT"/deserts_combined_hit.ids                    > $tmp/u_prev.ids
sort -u "$OUT"/mgnify90_hit.ids "$OUT"/gem_hit.ids "$OUT"/deserts_combined_hit.ids "$OUT"/spire_hit.ids > $tmp/u_all.ids
sort -u "$OUT"/gem_hit.ids "$OUT"/spire_hit.ids                                                        > $tmp/u_genomes.ids
{
  echo -e "union\tset\tn_with_hit\tn_absent\tn_total\tpct_absent\tcomment"
  for set in dark ctrl; do
    tot=$(wc -l < $tmp/$set.s)
    for u in u_prev u_all u_genomes; do
      case $u in
        u_prev)    lbl="mgnify90+gem+deserts_combined"; cmt="pre-SPIRE union (superseded as headline)";;
        u_all)     lbl="mgnify90+gem+deserts+spire";    cmt="HEADLINE: all named references, SPIRE included";;
        u_genomes) lbl="gem+spire";                     cmt="known GENOME space only (both genome-resolved)";;
      esac
      h=$(comm -12 $tmp/$set.s $tmp/$u.ids | wc -l)
      awk -v l="$lbl" -v s=$set -v h=$h -v t=$tot -v c="$cmt" \
        'BEGIN{printf "%s\t%s\t%d\t%d\t%d\t%.2f\t%s\n",l,s,h,t-h,t,100*(t-h)/t,c}'
    done
  done
} > "$U"
comm -23 $tmp/dark.s $tmp/u_all.ids > "$OUT/dark_nohit_union_all.ids"
comm -23 $tmp/ctrl.s $tmp/u_all.ids > "$OUT/ctrl_nohit_union_all.ids"
log "union table -> $U"; cat "$U"
log "dark no-hit anywhere: $(wc -l < "$OUT/dark_nohit_union_all.ids")"
rm -rf $tmp
log "DONE. Index left at ${DB}.dmnd ($(du -sh "$BUILD" 2>/dev/null | cut -f1)) — keep-or-delete is Alejandra's call."
