#!/bin/bash
# Phase E.1 — for each ecology dark family (23,382 reps), pull its BEST C·V·T function assignment from
# candidate_scores_v3.tsv (11G), keyed on the family rep. Small output to bridge decoding<->ecology.
#   sbatch scripts/ibex/97_darkfam_decode_join.sh
#SBATCH --job-name=hohndor_darkfam_decode
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_decode_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_decode_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DM=$WD/06_darkmatter
OUT=$DM/darkfam_ecology
SC=$DM/scores_v3/candidate_scores_v3.tsv
META=$OUT/catalog/darkfam_meta.tsv

# best assignment per family rep (col2=rep, col28=priority); keep ec/rxn/pathways/confidence/secondary
awk -F'\t' '
  NR==FNR { if (FNR>1) fam[$1]=1; next }                       # family ids (col1 of meta)
  FNR==1 { next }
  ($1 in fam) {                                                  # col1 = candidate_protein = family rep
    pr = $28 + 0
    if (!($1 in best) || pr > best[$1]) {
      best[$1] = pr
      row[$1] = $5 "\t" $7 "\t" $8 "\t" $25 "\t" $19 "\t" $28   # missing_ec rxn_name pathways C_conf core_secondary priority
    }
  }
  END {
    print "family_id\tbest_ec\trxn_name\tpathways\tC_confidence\tcore_secondary\tpriority"
    for (f in row) print f "\t" row[f]
  }' "$META" "$SC" > "$OUT/family_ec.tsv"

echo "families with an assignment: $(tail -n +2 "$OUT/family_ec.tsv" | wc -l) / $(tail -n +2 "$META" | wc -l)"
echo "with a specific EC: $(tail -n +2 "$OUT/family_ec.tsv" | awk -F"\t" '$2!="" && $2!="NA"' | wc -l)"
date -Iseconds > "$OUT/DECODE_JOIN_DONE"
echo "[$(date -Iseconds)] decode join done -> $OUT/family_ec.tsv"
