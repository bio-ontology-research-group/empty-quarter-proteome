#!/bin/bash
# 05c_afdb — AFDB v6 leg of the validation screen, split into NPART parallel tasks.
# The PDB leg is already banked (sc5c_parts/pdb_hits.raw.tsv, job 50603810 before it timed out);
# only AFDB (241,070,489 structures) is rerun. Single-job 05c_screen.sh TIMED OUT at 24h mid
# split 2/2 — splitting the 20,905 queries quarters the per-task prefilter work (~9h each).
# Params verbatim from rebuild_20260729/scripts/167_postfold (foldseek/10-941cd33).
# Per-part checkpointing: a task with a banked part exits immediately, so requeues are cheap.
# Submit: sbatch --array=1-4 scripts/05c_afdb.sh
#SBATCH --job-name=sc5c_afdb
#SBATCH --account=pi-hohndor
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=48
#SBATCH --mem=768G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5c_afdb_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5c_afdb_%A_%a.err
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
V=$M/04_structure/00_structure_clean/05_validate
AFDB=$M/rebuild_20260729/structural_calibration/references/afdb_uniprot_v6_uniprot_2025_03
P=/ibex/user/velazqal/sc5c_parts
FORMAT='query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,qcov,tcov,qlen,tlen,alnlen'
NPART=4
N=${SLURM_ARRAY_TASK_ID:?array task id required}
test -f "$AFDB/REFERENCE_COMPLETE"
[ -d "$V/query_pdb" ] || { echo "ABORT: missing $V/query_pdb"; exit 2; }

OUT=$P/afdb_hits.part${N}.tsv
[ -s "$OUT" ] && { echo "part $N already banked ($(wc -l < "$OUT") rows) — skip"; exit 0; }

# deterministic 1-of-NPART slice of the query models; symlinks resolved to real paths
QD=$P/q_part${N}; rm -rf "$QD"; mkdir -p "$QD" "$P/tmp_${N}"
ls "$V/query_pdb" | sort | awk -v n="$N" -v k="$NPART" 'NR % k == (n % k)' > "$P/part${N}.list"
while read -r f; do ln -s "$(readlink -f "$V/query_pdb/$f")" "$QD/$f"; done < "$P/part${N}.list"
NQ=$(ls "$QD" | wc -l)
echo "[$(date -Iseconds)] part $N: $NQ queries of 20,905"
[ "$NQ" -eq "$(wc -l < "$P/part${N}.list")" ] || { echo "ABORT: link count mismatch"; exit 2; }

module load foldseek/10-941cd33
foldseek easy-search "$QD" "$AFDB/database/afdb" "$OUT.tmp" "$P/tmp_${N}" \
  --alignment-type 1 -e 10 -s 9.5 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G \
  --threads "${SLURM_CPUS_PER_TASK:-48}" --format-output "$FORMAT"
mv "$OUT.tmp" "$OUT"
rm -rf "$P/tmp_${N}" "$QD"
echo "[$(date -Iseconds)] DONE part $N: $(wc -l < "$OUT") hit rows"
