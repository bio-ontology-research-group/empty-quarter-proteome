#!/bin/bash
# Phase A.2 — recruit each soil sample's reads to the dark-family gene catalogue (CoverM contig
# mode). Per gene: TPM + covered_fraction. One array task per sample.
#   sbatch --dependency=afterok:<90jobid> --array=0-550%30 scripts/ibex/91_darkfam_coverm_array.sh
# (array upper bound is set from samples_darkfam.tsv; use $(( $(wc -l < samples) - 2 )) )
#SBATCH --job-name=hohndor_darkfam_cov
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_cov_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_cov_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_ecology
CAT=$OUT/catalog/darkfam_catalog.fna
MAN=$OUT/samples_darkfam.tsv
T=${SLURM_CPUS_PER_TASK:-16}
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes

# +2: skip header (line 1) and convert 0-based array id to 1-based data line
LINE=$(awk -v i=$((SLURM_ARRAY_TASK_ID+2)) 'NR==i' "$MAN" | tr -d '\r')
SAMPLE=$(echo "$LINE" | cut -f1); R1=$(echo "$LINE" | cut -f2); R2=$(echo "$LINE" | cut -f3)
[ -n "$SAMPLE" ] || { echo "no sample at task $SLURM_ARRAY_TASK_ID"; exit 0; }
SOUT=$OUT/cov/$SAMPLE
mkdir -p "$SOUT"
if [ -f "$SOUT/DONE" ]; then echo "$SAMPLE already done"; exit 0; fi
echo "[$(date -Iseconds)] coverM $SAMPLE"

coverm contig --reference "$CAT" \
    --coupled "$R1" "$R2" \
    -m tpm covered_fraction \
    -p minimap2-sr -t "$T" \
    --output-file "$SOUT/cov.tsv" 2> "$SOUT/coverm.log"
touch "$SOUT/DONE"
echo "[$(date -Iseconds)] done $SAMPLE"
