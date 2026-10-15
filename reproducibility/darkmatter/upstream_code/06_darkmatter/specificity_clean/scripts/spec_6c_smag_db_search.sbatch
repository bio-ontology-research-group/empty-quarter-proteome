#!/bin/bash
# §3.1 SMAG leg, step c: protein DB, the one-pass search, scoring, size-only expectation, union.
#
# SMAG = 40,039 soil MAGs (Ma et al. 2023), genome-level like GEM and the SPIRE reps, gene-called here
# (spec_6b). The comparison has to be the same as every other row, so everything below is reused:
#   * DIAMOND v2.0.15 from the `eggnog` env, the binary behind every existing index and search
#   * rule: diamond blastp --sensitive --evalue 1e-3 -k 1, one pass of the same 212,084 queries
#     (106,042 dark reps + 106,042 matched annotated controls, results/query_all.faa)
#   * score() copied from spec_3_spire (separate `local` lines; see the note there)
#   * size-only expectation: raref_curve.py on the SAME permissive MGnify90 hit table
#     (spec_4, -e 1 -k 500), with SMAG's exact protein count added as an anchor. No interpolation.
# Storage follows the SPIRE leg: the index lives on own quota (/ibex/user/velazqal/eqmags_31/smag/),
# only the .dmnd + makedb.log + manifest are kept (`diamond getseq` recovers the proteins). The
# genomes (tmp/smag/), protein chunks (tmp/smag_proteins/) and the concatenated FASTA are deleted once
# makedb's sequence count equals the FASTA's header count.
# Outputs are NEW files; specificity_summary.tsv, union_specificity_summary.tsv, rarefaction_curve.tsv
# and content_vs_depth.tsv are not touched.
#   sbatch spec_6c_smag_db_search.sbatch
#SBATCH --job-name=ale_spec_smag_search
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_search_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_search_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
QUERY=$OUT/query_all.faa
LIST=$S/tmp/smag/genome_list.txt
CHUNKS=$S/tmp/smag_proteins
BUILD=/ibex/user/velazqal/eqmags_31/smag
DB=$BUILD/smag_db
MKLOG=$BUILD/makedb.log
MAN=$BUILD/smag_db_manifest.tsv
RAREF_M8=/ibex/user/velazqal/eqmags_31/raref/raref_mgnify90_e1_k500.m8
mkdir -p "$BUILD/tmp"
source $WD/../miniconda3/etc/profile.d/conda.sh; conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
log(){ echo "[$(date -Iseconds)] $*"; }
log "diamond: $(diamond --version 2>&1 | head -1)"

# ---------------- build (skipped once the manifest exists) ----------------
if [ ! -s "$MAN" ]; then
  NG=$(wc -l < "$LIST")
  NC_EXP=$(( (NG + 499) / 500 ))
  NC=$(find "$CHUNKS" -maxdepth 1 -name 'chunk_*.faa' -size +0 | wc -l)
  log "genomes $NG, protein chunks $NC / $NC_EXP"
  [ "$NC" -ne "$NC_EXP" ] && { log "ABORT: gene calling incomplete — rerun spec_6b"; exit 1; }
  FAA=$S/tmp/smag_all.faa
  find "$CHUNKS" -maxdepth 1 -name 'chunk_*.faa' -print0 | LC_ALL=C sort -z | xargs -0 cat > "$FAA" \
    || { log "ABORT: concat failed"; exit 1; }
  read -r NMAG NPROT NRES < <(awk '/^>/ {n++; split(substr($1,2),a,"::"); m[a[1]]=1; next} {r+=length($0)}
      END {printf "%d %d %d\n", length(m), n, r}' "$FAA")
  log "proteins $NPROT from $NMAG MAGs, $NRES residues"
  [ "$NMAG" -ne "$NG" ] && log "NOTE: $((NG - NMAG)) genomes yielded no protein"
  diamond makedb --in "$FAA" -d "$DB" -p "$T" 2>&1 | tee "$MKLOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { log "ABORT: makedb failed — inputs kept"; exit 1; }
  NSEQ=$(grep -oP 'Database sequences\s+\K[0-9]+' "$MKLOG" | tail -1)
  [ "${NSEQ:-0}" -ne "$NPROT" ] && { log "ABORT: makedb $NSEQ != FASTA $NPROT — inputs kept"; exit 1; }
  DMD5=$(md5sum "${DB}.dmnd" | cut -d' ' -f1)
  {
    echo -e "key\tvalue"
    echo -e "source\tZenodo 8223844 (Ma et al. 2023 Nat Commun 14:7318), CC-BY-4.0"
    echo -e "gene_caller\tpyrodigal 3.7.1 -p meta --no-stop-codon"
    echo -e "diamond\t$(diamond --version 2>&1 | head -1)"
    echo -e "genomes_listed\t$NG"
    echo -e "mags_with_proteins\t$NMAG"
    echo -e "proteins\t$NPROT"
    echo -e "residues\t$NRES"
    echo -e "db\t${DB}.dmnd"
    echo -e "db_bytes\t$(stat -c %s "${DB}.dmnd")"
    echo -e "db_md5\t$DMD5"
  } > "$MAN"
  cp "$MAN" "$OUT/smag_db_manifest.tsv"
  log "manifest -> $MAN"; cat "$MAN"
  log "deleting genomes, chunks and the concatenated FASTA ($(du -sh "$S/tmp" | cut -f1) under tmp/)"
  rm -rf "$S/tmp/smag" "$CHUNKS" "$FAA"
  rmdir "$S/tmp" 2>/dev/null
else
  log "manifest exists — index built, inputs already deleted"
fi
NSEQ=$(awk -F'\t' '$1=="proteins"{print $2}' "$MAN")
log "SMAG proteins: $NSEQ"

# ---------------- search, same rule as every other catalogue ----------------
SUMM="$OUT/smag_specificity_summary.tsv"
[ -f "$SUMM" ] || echo -e "catalog\tdb_seqs\tset\tn_with_hit\tn_absent\tn_total\tpct_absent" > "$SUMM"

score(){ # name  db_seqs   (as spec_3_spire: keep the `local` lines separate)
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

if [ -s "$OUT/smag_hits.m8" ]; then
  log "[smag] hits exist — scoring only"
else
  log "[smag] diamond blastp vs $DB (--sensitive -e1e-3 -k1)"
  /usr/bin/time -v diamond blastp -q "$QUERY" -d "$DB" --sensitive --evalue 1e-3 -k 1 -p "$T" \
     -f 6 qseqid sseqid pident evalue -o "$OUT/smag_hits.m8" 2> "$OUT/smag.log" \
     || { log "[smag] blastp FAILED"; exit 1; }
fi
score smag "$NSEQ"

# ---------------- size-only expectation at SMAG's exact size ----------------
CURVE=$OUT/rarefaction_curve_with_smag.tsv
if [ ! -s "$CURVE" ]; then
  python3 "$S/scripts/raref_curve.py" "$RAREF_M8" "$B/FAMILIES_FINAL.txt" "$B/CONTROL_FAMILIES.txt" \
     "$CURVE" "smag:$NSEQ" || { log "curve FAILED"; exit 1; }
fi
# self-check: every pre-existing grid row must reproduce rarefaction_curve.tsv exactly
DIFF=$(awk -F'\t' 'NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} FNR>1 && ($1"\t"$2 in a) && a[$1"\t"$2]!=$4"\t"$5 {n++} END{print n+0}' \
       "$OUT/rarefaction_curve.tsv" "$CURVE")
log "curve rows differing from rarefaction_curve.tsv: $DIFF (must be 0)"
[ "$DIFF" -ne 0 ] && { log "ABORT: curve does not reproduce"; exit 1; }
awk -F'\t' -v OFS='\t' '
  NR==FNR { if ($1=="smag" && $3=="dark") d=$7; if ($1=="smag" && $3=="ctrl") c=$7; n=$2; next }
  FNR>1 && $3 ~ /(^|,)smag(,|$)/ { de=$4; ce=$5 }
  END {
    print "catalogue","db_seqs","dark_obs","dark_exp_size_only","dark_gap_pts","ctrl_obs","ctrl_exp_size_only","ctrl_gap_pts"
    printf "smag\t%s\t%.2f\t%.2f\t%+.2f\t%.2f\t%.2f\t%+.2f\n", n, d, de, de-d, c, ce, ce-c
  }' "$SUMM" "$CURVE" > "$OUT/smag_content_vs_depth.tsv"
log "content vs depth:"; cat "$OUT/smag_content_vs_depth.tsv"

# ---------------- union with SMAG, a separate file ----------------
U="$OUT/union_specificity_with_smag_summary.tsv"
tmp=$(mktemp -d)
sort -u "$B/FAMILIES_FINAL.txt"  > $tmp/dark.s
sort -u "$B/CONTROL_FAMILIES.txt" > $tmp/ctrl.s
sort -u "$OUT"/mgnify90_hit.ids "$OUT"/gem_hit.ids "$OUT"/deserts_combined_hit.ids "$OUT"/spire_hit.ids > $tmp/u_all.ids
sort -u $tmp/u_all.ids "$OUT"/smag_hit.ids                                                             > $tmp/u_all_smag.ids
sort -u "$OUT"/gem_hit.ids "$OUT"/spire_hit.ids                                                        > $tmp/u_genomes.ids
sort -u $tmp/u_genomes.ids "$OUT"/smag_hit.ids                                                         > $tmp/u_genomes_smag.ids
{
  echo -e "union\tset\tn_with_hit\tn_absent\tn_total\tpct_absent\tcomment"
  for set in dark ctrl; do
    tot=$(wc -l < $tmp/$set.s)
    for u in u_all u_all_smag u_genomes u_genomes_smag; do
      case $u in
        u_all)          lbl="mgnify90+gem+deserts+spire";      cmt="existing headline, recomputed as a check";;
        u_all_smag)     lbl="mgnify90+gem+deserts+spire+smag"; cmt="all named references plus SMAG";;
        u_genomes)      lbl="gem+spire";                       cmt="existing known GENOME space, recomputed as a check";;
        u_genomes_smag) lbl="gem+spire+smag";                  cmt="known GENOME space plus SMAG (all genome-resolved)";;
      esac
      h=$(comm -12 $tmp/$set.s $tmp/$u.ids | wc -l)
      awk -v l="$lbl" -v s=$set -v h=$h -v t=$tot -v c="$cmt" \
        'BEGIN{printf "%s\t%s\t%d\t%d\t%d\t%.2f\t%s\n",l,s,h,t-h,t,100*(t-h)/t,c}'
    done
  done
} > "$U"
comm -23 $tmp/dark.s $tmp/u_all_smag.ids > "$OUT/dark_nohit_union_with_smag.ids"
comm -23 $tmp/ctrl.s $tmp/u_all_smag.ids > "$OUT/ctrl_nohit_union_with_smag.ids"
rm -rf $tmp
log "union -> $U"; cat "$U"
# the recomputed rows must match the shipped union file
for lbl in "mgnify90+gem+deserts+spire" "gem+spire"; do
  for set in dark ctrl; do
    a=$(awk -F'\t' -v l="$lbl" -v s=$set '$1==l && $2==s {print $4}' "$OUT/union_specificity_summary.tsv")
    b=$(awk -F'\t' -v l="$lbl" -v s=$set '$1==l && $2==s {print $4}' "$U")
    [ "$a" = "$b" ] && log "check $lbl/$set: $b absent, matches" || log "CHECK FAILED $lbl/$set: shipped $a vs $b"
  done
done
log "kept: $(du -sh "$BUILD" | cut -f1) in $BUILD"
log "DONE_SMAG_SEARCH"
