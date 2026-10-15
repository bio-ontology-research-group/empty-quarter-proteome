#!/bin/bash
# bgc_clean — antiSMASH BGC census on the CURRENT catalogue (5,438 dRep reps), §2.4 recompute.
# Same design as the historical 86_antismash_array.sh (antismash/8.0.4 module, --minimal =
# core cluster detection only, prodigal gene calls) but on results/rep_index.tsv (mag \t fna).
# Census only — BGC novelty % is an upper bound (contig fragmentation), do NOT headline.
#
# Submit: sbatch --array=1-272%60 scripts/bgc_antismash_array.sh   # 5,438 / BATCH=20
#
#SBATCH --job-name=ale_bgc_census
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism/bgc_clean/logs/antismash_%A_%a.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --partition=batch

set -uo pipefail
module load antismash/8.0.4 2>/dev/null || { echo "antismash module load failed"; exit 2; }
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism/bgc_clean
INDEX=$B/results/rep_index.tsv
OUTBASE=$B/results/antismash
BATCH=20
mkdir -p "$OUTBASE"

total=$(wc -l < "$INDEX")
start=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH + 1 )); end=$(( SLURM_ARRAY_TASK_ID * BATCH ))
[ "$end" -gt "$total" ] && end=$total
[ "$start" -gt "$total" ] && { echo "nothing to do"; exit 0; }
echo "[$(date -Iseconds)] task $SLURM_ARRAY_TASK_ID rows $start..$end of $total"
antismash --version || true

sed -n "${start},${end}p" "$INDEX" | while IFS=$'\t' read -r mag fna; do
    ODIR=$OUTBASE/$mag
    if [ -s "$ODIR/${mag}.json" ]; then echo "[$mag] done — skip"; continue; fi
    [ -s "$fna" ] || { echo "[$mag] missing fna $fna"; continue; }
    rm -rf "$ODIR"
    antismash --minimal --genefinding-tool prodigal --cpus "$SLURM_CPUS_PER_TASK" \
        --output-dir "$ODIR" --output-basename "$mag" "$fna" 2>&1 | tail -2 || echo "[$mag] FAILED"
done
echo "[$(date -Iseconds)] task $SLURM_ARRAY_TASK_ID complete"
