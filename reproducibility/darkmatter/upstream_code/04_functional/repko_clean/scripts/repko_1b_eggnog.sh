#!/bin/bash
# repko step 1b -- EGGNOG-MAPPER only, after step 1a (bakta) has produced the proteins.
# Writes into the per-sample tree at the canonical path:
#     <sample>/wgs/annotation/<mode>/eggnog/<bin>/
# Identical stack to 04_eggnog_array.sh: data_dir 2.1.12, --itype proteins --override.
#
# --dbmem loads the 41 GB eggnog.db into RAM once, then serves all annotation from memory.
# Without it, the annotation phase random-reads the shared file over WekaFS and, at >=2 tasks/node,
# blocks in the kernel (D state on wekafs_lock_helper, zero CPU). --dbmem is the standard way to run
# emapper on a cluster and eliminates both problems (speed AND the freeze).
# First task on a node pays ~1-3 min DB load; later tasks on the same node reuse the page cache.
#
# MEMORY: eggnog.db is 41 GB and the annotation phase random-reads it. A first attempt at
# 32 G with several tasks packed on the same node stalled in uninterruptible I/O wait
# (D state, wekafs_lock_helper, zero CPU) -- the file pages are charged to the cgroup, so a
# 32 G ceiling against a 41 GB DB thrashes. There is in-project precedent for eggNOG at 128 G
# (job eq_eggnog). We therefore request 128 G and keep concurrency LOW so tasks do not pack
# onto one node. Runtime reference: a complete run reports "## Rate: 3.66 q/s" -> ~20-25 min
# for a 3-5k protein bin; 6 h is a wide margin.
#
#   N=$(tail -n +2 results/annotate_targets.tsv | wc -l)
#   sbatch --array=1-${N}%2 scripts/repko_1b_eggnog.sh
#
#SBATCH --job-name=repko_eggnog
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/eggnog_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/eggnog_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=batch

set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
W=$M/04_functional/repko_clean
EGGNOG_DB=/ibex/reference/KSL/eggnog-mapper/2.1.12

LINE=$(tail -n +2 "$W/results/annotate_targets.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
FA=$(echo "$LINE" | cut -f1); SAMPLE=$(echo "$LINE" | cut -f2)
MODE=$(echo "$LINE" | cut -f3); BIN=$(echo "$LINE" | cut -f4)
WGS=$(readlink -f "$(dirname "$FA")/../../..")
PROT="$WGS/annotation/$MODE/bakta/$BIN/$BIN.faa"
EOUT="$WGS/annotation/$MODE/eggnog/$BIN"
echo "[$(date -Iseconds)] $SAMPLE/$MODE/$BIN"

[ -s "$PROT" ] || { echo "ERROR: bakta proteins missing ($PROT) -- run repko_1a_bakta.sh first"; exit 2; }

# A COMPLETE emapper run ends with a "## Rate:" trailer. Anything else is a truncated file from a
# timeout/kill and MUST be redone -- a bare -s test would silently accept it (this is exactly how
# the original pipeline left bins with partial/absent annotation and nobody noticed).
if [ -s "$EOUT/$BIN.emapper.annotations" ] && tail -3 "$EOUT/$BIN.emapper.annotations" | grep -q '^## Rate:'; then
    echo "  already done (complete) -- skip"; exit 0
fi
if [ -e "$EOUT/$BIN.emapper.annotations" ]; then
    echo "  WARNING: partial emapper output found -- discarding and re-running"
    rm -f "$EOUT/$BIN.emapper."*
fi

source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate eggnog
mkdir -p "$EOUT"; cd "$EOUT"
emapper.py -i "$PROT" --output "$BIN" --output_dir "$EOUT" --data_dir "$EGGNOG_DB" --dbmem \
           --cpu "${SLURM_CPUS_PER_TASK:-16}" --itype proteins --override

tail -3 "$EOUT/$BIN.emapper.annotations" 2>/dev/null | grep -q '^## Rate:' \
    || { echo "ERROR: eggnog output for $BIN is missing or truncated (no '## Rate:' trailer)"; exit 4; }
echo "[$(date -Iseconds)] DONE $SAMPLE/$BIN  annotated=$(grep -vc '^#' "$EOUT/$BIN.emapper.annotations")  $(tail -1 "$EOUT/$BIN.emapper.annotations")"
