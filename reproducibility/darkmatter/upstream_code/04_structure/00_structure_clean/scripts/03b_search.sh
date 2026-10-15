#!/bin/bash
# 03b_search v2 — FoldSeek 3Di search of each ~15k sub-chunk vs AFDB v6 + PDB.
# mmap (--db-load-mode 1), no split, largemem. One array task per manifest line.
# Submit AFTER 03b_prep, with the real count:  sbatch --array=1-<NCH>%24 03b_search.sh
#SBATCH --job-name=sc3_search
#SBATCH --account=pi-hohndor
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_search_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc3_search_%A_%a.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=768G
#SBATCH --time=12:00:00
#SBATCH --partition=largemem
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
SC=$M/04_structure/00_structure_clean; O=$SC/03_search
FS=$M/tools/foldseek/bin/foldseek
AFDB=$M/rebuild_20260729/structural_calibration/references/afdb_uniprot_v6_uniprot_2025_03/database/afdb
PDB=$M/04_structure/fsdb_pdb
[ -f "${AFDB}.dbtype" ] || { echo "AFDB v6 missing: $AFDB"; exit 2; }
[ -f "${PDB}.dbtype" ]  || { echo "PDB missing: $PDB"; exit 2; }

Q=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$O/manifest.txt")
[ -f "${Q}_ss.dbtype" ] || { echo "query 3Di missing: $Q"; exit 2; }
tag=$(basename "$Q"); W=$O/hits/$tag; mkdir -p "$W"
[ -s "$W/hit_ids.txt" ] && { echo "[$tag] hits exist — skip"; exit 0; }
TA=$O/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$TA/a" "$TA/p"

P=(--alignment-type 2 -s 9.5 -e 1e-3 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G --threads "${SLURM_CPUS_PER_TASK:-48}")
echo "[$(date -Iseconds)] [$tag] search vs AFDB v6 ($(wc -l < "${Q}.index") queries)"
"$FS" search "$Q" "$AFDB" "$W/a_res" "$TA/a" "${P[@]}"
"$FS" createtsv "$Q" "$AFDB" "$W/a_res" "$W/afdb_hits.tsv"
echo "[$(date -Iseconds)] [$tag] search vs PDB"
"$FS" search "$Q" "$PDB" "$W/p_res" "$TA/p" "${P[@]}"
"$FS" createtsv "$Q" "$PDB" "$W/p_res" "$W/pdb_hits.tsv"

cut -f1 "$W/afdb_hits.tsv" "$W/pdb_hits.tsv" 2>/dev/null | sort -u > "$W/hit_ids.txt"
rm -rf "$TA" "$W"/a_res* "$W"/p_res*
echo "[$(date -Iseconds)] DONE $tag  hits=$(wc -l < "$W/hit_ids.txt")"
