#!/bin/bash
# resistome_clean — AMRFinderPlus (--plus: AMR + stress/virulence) on the CURRENT catalogue
# (5,438 dRep reps), §2.4 recompute. Same design as the historical 94_amrfinder_array.sh but
# on results/rep_faa_index.tsv (mag \t faa.gz; decompressed to $TMPDIR per MAG).
# Task 1 records tool + DB versions to results/VERSIONS.txt (Methods M.10).
#
# Submit: sbatch --array=1-37 scripts/amr_array.sh   # 5,438 / BATCH=150
#
#SBATCH --job-name=ale_amrfind
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism/resistome_clean/logs/amrfind_%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --partition=batch

set -uo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate amrfinder
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism/resistome_clean
INDEX=$B/results/rep_faa_index.tsv
OUTBASE=$B/results/amrfinder
BATCH=150
mkdir -p "$OUTBASE"

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ] && [ ! -s "$B/results/VERSIONS.txt" ]; then
    { amrfinder --version; amrfinder --database_version 2>&1 | tail -3; } > "$B/results/VERSIONS.txt" || true
fi

total=$(wc -l < "$INDEX")
start=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH + 1 )); end=$(( SLURM_ARRAY_TASK_ID * BATCH ))
[ "$end" -gt "$total" ] && end=$total
[ "$start" -gt "$total" ] && { echo "nothing to do"; exit 0; }
echo "[$(date -Iseconds)] task $SLURM_ARRAY_TASK_ID rows $start..$end of $total"

sed -n "${start},${end}p" "$INDEX" | while IFS=$'\t' read -r mag faagz; do
    OUT=$OUTBASE/${mag}.amrfinder.tsv
    [ -s "$OUT" ] && { echo "[$mag] done — skip"; continue; }
    [ -s "$faagz" ] || { echo "[$mag] missing faa $faagz"; continue; }
    TMPFAA=${TMPDIR:-/tmp}/${mag}.faa
    zcat "$faagz" > "$TMPFAA"
    amrfinder -p "$TMPFAA" --plus --name "$mag" --threads "$SLURM_CPUS_PER_TASK" -o "$OUT" \
        2>&1 | tail -1 || echo "[$mag] FAILED"
    rm -f "$TMPFAA"
done
echo "[$(date -Iseconds)] task $SLURM_ARRAY_TASK_ID complete"
