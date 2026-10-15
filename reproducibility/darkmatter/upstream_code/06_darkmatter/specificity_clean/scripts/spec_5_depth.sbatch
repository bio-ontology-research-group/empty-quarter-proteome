#!/bin/bash
# §3.1 depth-robustness leg — the reviewer-proofing question: "did you check the DEEPEST, NEWEST reference?"
#
# QUESTION. 30,479 dark families have no DIAMOND hit in any of MGnify90 (717.7M) + GEM (123.0M) +
# deserts_combined (147.5M) + SPIRE reps (269.5M) = 1.258 G sequences. Does a 2.3x deeper MGnify snapshot
# explain them? Search that no-hit set alone against MGnify Proteins **2026_07**, the complete MGnify90
# representative set: **1,661,142,385 sequences**, prebuilt MMseqs2 DB (~427 GB) at
#   rebuild_20260729/protein_novelty/references/mgnify_proteins_2026_07_mgy90/mmseqs/mgy90
# The 435 control families with no hit anywhere go in the same pass, so the control floor travels along.
#
# RECIPE IS THE PI's, DELIBERATELY. `-s 7.5 -e 1e-3 --max-seqs 10`, copied from his
# `rebuild_20260729/scripts/47_search_current_mgnify.sbatch`, so our answer is directly comparable with
# his protein-novelty rerun (`protein_novelty/mgnify_2026_07/integrated/mgnify_family_summary.json`:
# 299,927 UniRef50-dark queries, 155,609 with any hit -> **48.1% still no-hit** at this depth).
# His own summary records the framing we independently locked: "no hit means not detected in this
# versioned reference under these thresholds; it does not establish geographic endemicity."
#
# TOOL CHANGE IS A REAL CAVEAT. This leg is MMseqs2, every other leg is DIAMOND v2.0.15, and the release
# differs from the DIAMOND mgnify90 index. So this is a ROBUSTNESS check, NOT a row in tab:specificity,
# and its number must never be unioned with the DIAMOND absences. Report it as: "of the N families absent
# from every DIAMOND catalogue, M remain undetected in the complete 1.66-billion-sequence MGnify90
# 2026_07 representative set searched with MMseqs2 at -s 7.5, e<=1e-3."
#
#   sbatch spec_5_depth.sbatch
#SBATCH --job-name=ale_spec_depth
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_depth_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_depth_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=1400G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
# TAG selects which no-hit sets go in (added 2026-09-24). Unset = the original run: inputs
#   {dark,ctrl}_nohit_union_all.ids, outputs without suffix. TAG=with_smag = the union that also
#   includes SMAG (spec_6c): inputs {dark,ctrl}_nohit_union_with_smag.ids, every output suffixed
#   _with_smag, own work dir. Same search and resources either way.
#   sbatch --export=ALL,TAG=with_smag spec_5_depth.sbatch
TAG=${TAG:-}
SUF=${TAG:+_$TAG}
IDS=${TAG:-all}
DARK_IDS=$OUT/dark_nohit_union_${IDS}.ids
CTRL_IDS=$OUT/ctrl_nohit_union_${IDS}.ids
W=/ibex/user/velazqal/eqmags_31/depth$SUF      # own quota; project volume is 99% full
REF=$WD/rebuild_20260729/protein_novelty/references/mgnify_proteins_2026_07_mgy90
MMSEQS=$WD/../miniconda3/envs/mmseqs2/bin/mmseqs
mkdir -p "$W/tmp" "$S/logs"
log(){ echo "[$(date -Iseconds)] $*"; }
log "mmseqs: $("$MMSEQS" version 2>&1 | head -1)"
log "set: ${TAG:-original} — $DARK_IDS ($(wc -l < "$DARK_IDS")) + $CTRL_IDS ($(wc -l < "$CTRL_IDS"))"

# query = the no-hit dark + ctrl reps pulled from the one-pass query FASTA
if [ ! -s "$W/query_nohit.faa" ]; then
  cat "$DARK_IDS" "$CTRL_IDS" | sort -u > "$W/nohit_all.ids"
  awk 'NR==FNR {k[$1]=1; next} /^>/ {p = (substr($1,2) in k)} p' "$W/nohit_all.ids" "$OUT/query_all.faa" \
    > "$W/query_nohit.faa"
fi
NQ=$(grep -c '^>' "$W/query_nohit.faa")
log "query: $NQ seqs (dark + ctrl with no DIAMOND hit anywhere)"
[ "$NQ" -ne "$(wc -l < "$W/nohit_all.ids")" ] && { log "ABORT: query FASTA does not match the id list"; exit 1; }

A=$(df -B1G --output=avail "$W" | tail -1 | tr -d ' ')
log "free on $W: ${A} GB"
[ "$A" -lt 100 ] && { log "ABORT: <100 GB free"; exit 1; }
[ -e "${REF}/mmseqs/mgy90.dbtype" ] || { log "ABORT: reference DB missing"; exit 1; }

if [ ! -e "$W/qdb.dbtype" ]; then
  log "createdb"
  "$MMSEQS" createdb "$W/query_nohit.faa" "$W/qdb" || { log "ABORT: createdb failed"; exit 1; }
fi

if [ ! -e "$W/res.dbtype" ]; then
  log "search vs mgy90 (1,661,142,385 seqs) -s 7.5 -e 1e-3 --max-seqs 10"
  "$MMSEQS" search "$W/qdb" "$REF/mmseqs/mgy90" "$W/res" "$W/tmp" \
      -s 7.5 -e 1e-3 --max-seqs 10 \
      --threads "${SLURM_CPUS_PER_TASK:-48}" \
      --split-memory-limit 1200G \
    || { log "ABORT: search failed"; exit 1; }
fi

TSV=$OUT/depth_mgnify2026_07_hits$SUF.tsv
if [ ! -s "$TSV" ]; then
  log "convertalis"
  "$MMSEQS" convertalis "$W/qdb" "$REF/mmseqs/mgy90" "$W/res" "$TSV" \
      --format-output "query,target,fident,evalue,bits,qcov,tcov" \
      --threads "${SLURM_CPUS_PER_TASK:-48}" \
    || { log "ABORT: convertalis failed"; exit 1; }
fi
log "hit rows: $(wc -l < "$TSV")"

# ---------------- score: how much does 2.3x more depth explain? ----------------
HITS=$OUT/depth_hit$SUF.ids
cut -f1 "$TSV" | sort -u > "$HITS"
SUMM=$OUT/depth_robustness_summary$SUF.tsv
{
  printf "set\tn_queried_nohit\tn_now_hit_at_1.66G\tn_still_nohit\tpct_still_nohit\tpct_of_full_universe\n"
  for set in dark ctrl; do
    if [ "$set" = dark ]; then src="$DARK_IDS"; else src="$CTRL_IDS"; fi
    q=$(wc -l < "$src")
    h=$(sort -u "$src" | comm -12 - "$HITS" | wc -l)
    awk -v s="$set" -v q="$q" -v h="$h" \
      'BEGIN{printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\n", s, q, h, q-h, 100*(q-h)/q, 100*(q-h)/106042}'
  done
} > "$SUMM"
# the surviving set is the object any "truly unmatched" statement should rest on
sort -u "$DARK_IDS" | comm -23 - "$HITS" > "$OUT/dark_nohit_plus_deepmgnify$SUF.ids"
log "written: $SUMM"; column -t "$SUMM" 2>/dev/null || cat "$SUMM"
log "dark families unmatched even at 1.66 G: $(wc -l < "$OUT/dark_nohit_plus_deepmgnify$SUF.ids")"
log "DONE."
