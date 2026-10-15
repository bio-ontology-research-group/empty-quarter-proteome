#!/bin/bash
# repko step 1 -- ANNOTATE the 124 raw bins that were never annotated, writing into the
# per-sample tree exactly where every other annotation lives:
#     <sample>/wgs/annotation/<mode>/bakta/<bin>/   and   .../eggnog/<bin>/
#
# Targets (results/annotate_targets.tsv, built from the unannotated reps):
#   80  promoted_LQ    -- raw bin was LQ pre-refinement, so it never made the HQ+MQ list
#   44  parent_of_sub  -- parent of a _sub representative; annotating it lets step 2 subset it
#
# Tool stack is IDENTICAL to the original 03_bakta_array.sh / 04_eggnog_array.sh so the new
# annotations are directly comparable with the existing 5,270:
#   bakta   --db /ibex/reference/KSL/bakta/6.0/db --skip-plot --force
#   emapper --data_dir /ibex/reference/KSL/eggnog-mapper/2.1.12 --itype proteins --override
#
# Bakta and eggNOG run in ONE task (the originals were two arrays; that split is exactly why
# 04_eggnog_array.sh could `exit 0` when Bakta output was absent, silently leaving bins
# unannotated). Both steps here are resumable and fail loudly.
#
# Resources are set from MEASURED usage of the original jobs (sacct): bakta ~5-10 min @24G/8cpu,
# eggnog ~5-10 min @32G/16cpu per bin; our 124 bins are small (mean 1.1 MB gz, largest 3.7 MB).
# 2h/32G/16cpu therefore carries a 3-6x margin and backfills far better than an oversized request.
#
#   N=$(tail -n +2 results/annotate_targets.tsv | wc -l)
#   sbatch --array=1-${N}%25 scripts/repko_1_annotate.sh
#
#SBATCH --job-name=repko_annot
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/annot_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/logs/annot_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=batch

set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
W=$M/04_functional/repko_clean
TARGETS=$W/results/annotate_targets.tsv
BAKTA_DB=/ibex/reference/KSL/bakta/6.0/db
EGGNOG_DB=/ibex/reference/KSL/eggnog-mapper/2.1.12

# row 1 is the header -> offset by one
LINE=$(tail -n +2 "$TARGETS" | sed -n "${SLURM_ARRAY_TASK_ID}p")
FA=$(echo "$LINE"     | cut -f1)
SAMPLE=$(echo "$LINE" | cut -f2)
MODE=$(echo "$LINE"   | cut -f3)
BIN=$(echo "$LINE"    | cut -f4)
[ -s "$FA" ] || { echo "ERROR: missing input FASTA $FA"; exit 2; }

# <...>/wgs/binning/<mode>/output_bins/<bin>.fa.gz  ->  <...>/wgs
WGS=$(readlink -f "$(dirname "$FA")/../../..")
BOUT="$WGS/annotation/$MODE/bakta/$BIN"
EOUT="$WGS/annotation/$MODE/eggnog/$BIN"
echo "[$(date -Iseconds)] $SAMPLE/$MODE/$BIN"
echo "  bakta  -> $BOUT"
echo "  eggnog -> $EOUT"

# ---------------------------------------------------------------- Bakta
if [ -s "$BOUT/$BIN.faa" ]; then
    echo "  bakta already done -- skip"
else
    source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate bakta
    mkdir -p "$BOUT"
    bakta --db "$BAKTA_DB" --output "$BOUT" --prefix "$BIN" --threads "${SLURM_CPUS_PER_TASK:-16}" \
          --skip-plot --force "$FA"
    conda deactivate 2>/dev/null || true
fi
[ -s "$BOUT/$BIN.faa" ] || { echo "ERROR: bakta produced no proteins for $BIN"; exit 3; }

# ---------------------------------------------------------------- eggNOG-mapper
# a COMPLETE emapper run ends with a "## Rate:" trailer; anything else is a truncated file
# from a timeout/kill and MUST be redone (a bare -s test would silently accept it)
if [ -s "$EOUT/$BIN.emapper.annotations" ] && tail -3 "$EOUT/$BIN.emapper.annotations" | grep -q '^## Rate:'; then
    echo "  eggnog already done (complete) -- skip"
else
    if [ -e "$EOUT/$BIN.emapper.annotations" ]; then
        echo "  WARNING: partial/truncated emapper output found -- discarding and re-running"
        rm -f "$EOUT/$BIN.emapper.annotations" "$EOUT/$BIN.emapper.hits" "$EOUT/$BIN.emapper.seed_orthologs"
    fi
    source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate eggnog
    mkdir -p "$EOUT"; cd "$EOUT"
    emapper.py -i "$BOUT/$BIN.faa" --output "$BIN" --output_dir "$EOUT" \
               --data_dir "$EGGNOG_DB" --cpu "${SLURM_CPUS_PER_TASK:-16}" --itype proteins --override
fi
# fail loudly rather than exit 0 (the defect in the original 04_eggnog_array.sh)
tail -3 "$EOUT/$BIN.emapper.annotations" 2>/dev/null | grep -q '^## Rate:' \
    || { echo "ERROR: eggnog output for $BIN is missing or truncated (no '## Rate:' trailer)"; exit 4; }

echo "[$(date -Iseconds)] DONE $SAMPLE/$BIN  proteins=$(grep -c '^>' "$BOUT/$BIN.faa")  annotated=$(grep -vc '^#' "$EOUT/$BIN.emapper.annotations")"
