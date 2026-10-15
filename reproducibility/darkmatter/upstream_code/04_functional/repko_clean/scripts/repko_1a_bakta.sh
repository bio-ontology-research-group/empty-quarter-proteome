#!/bin/bash
# repko step 1a -- BAKTA only, for the 124 raw bins that were never annotated.
# Split from eggNOG (as in the originals 03_bakta_array.sh / 04_eggnog_array.sh): bakta is light and
# fast, eggNOG is the memory-hungry one, so they must not share a task's resource request.
# Writes into the per-sample tree at the canonical path:
#     <sample>/wgs/annotation/<mode>/bakta/<bin>/
# Identical stack to 03_bakta_array.sh: db 6.0, --skip-plot --force.
# Resources from the original ogb_bakta jobs: 8 cpu / 24 G, 4-10 min observed -> 1 h is ample.
#
#   N=$(tail -n +2 results/annotate_targets.tsv | wc -l)
#   sbatch --array=1-${N}%10 scripts/repko_1a_bakta.sh
#
#SBATCH --job-name=repko_bakta
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/bakta_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/bakta_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --partition=batch

set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
W=$M/04_functional/repko_clean
BAKTA_DB=/ibex/reference/KSL/bakta/6.0/db

LINE=$(tail -n +2 "$W/results/annotate_targets.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
FA=$(echo "$LINE" | cut -f1); SAMPLE=$(echo "$LINE" | cut -f2)
MODE=$(echo "$LINE" | cut -f3); BIN=$(echo "$LINE" | cut -f4)
[ -s "$FA" ] || { echo "ERROR: missing input FASTA $FA"; exit 2; }

WGS=$(readlink -f "$(dirname "$FA")/../../..")
BOUT="$WGS/annotation/$MODE/bakta/$BIN"
echo "[$(date -Iseconds)] $SAMPLE/$MODE/$BIN -> $BOUT"

# bakta writes <BIN>.tsv last-ish; require the .faa (what eggNOG consumes) as the completion marker
if [ -s "$BOUT/$BIN.faa" ]; then
    echo "  already done ($(grep -c '^>' "$BOUT/$BIN.faa") proteins) -- skip"; exit 0
fi
source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate bakta
mkdir -p "$BOUT"
bakta --db "$BAKTA_DB" --output "$BOUT" --prefix "$BIN" --threads "${SLURM_CPUS_PER_TASK:-8}" \
      --skip-plot --force "$FA"
[ -s "$BOUT/$BIN.faa" ] || { echo "ERROR: bakta produced no proteins for $BIN"; exit 3; }
echo "[$(date -Iseconds)] DONE $SAMPLE/$BIN  proteins=$(grep -c '^>' "$BOUT/$BIN.faa")"
