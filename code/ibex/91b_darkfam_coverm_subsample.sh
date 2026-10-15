#!/bin/bash
# Recover the handful of samples whose reads make minimap2 pathologically slow (low-complexity
# content explodes chaining; same file size as fast samples but >14 h with no progress). TPM is
# depth-normalised, so a random subsample gives an unbiased relative-abundance estimate. Subsample
# to a fixed read budget, then run CoverM as usual.
#   sbatch --array=0-2 scripts/ibex/91b_darkfam_coverm_subsample.sh
#SBATCH --job-name=hohndor_darkfam_cov_sub
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_cov_sub_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_cov_sub_%A_%a.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_ecology
CAT=$OUT/catalog/darkfam_catalog.fna
MAN=$OUT/samples_darkfam.tsv
T=${SLURM_CPUS_PER_TASK:-32}
NPAIRS=30000000          # subsample budget (read pairs)
SEED=100
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes

SAMPLES=(39Sr2 F4Sr1 F4Sr3)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R1=$(awk -F"\t" -v s="$S" '$1==s{print $2}' "$MAN"); R2=$(awk -F"\t" -v s="$S" '$1==s{print $3}' "$MAN")
SOUT=$OUT/cov/$S; mkdir -p "$SOUT"
TMP=$OUT/_subtmp/$S; mkdir -p "$TMP"
echo "[$(date -Iseconds)] subsample+coverM $S (budget ${NPAIRS} pairs)"

seqtk sample -s$SEED "$R1" $NPAIRS | gzip > "$TMP/R1.fastq.gz" &
seqtk sample -s$SEED "$R2" $NPAIRS | gzip > "$TMP/R2.fastq.gz" &
wait
echo "[$(date -Iseconds)] subsample done: $(ls -lh "$TMP"/R1.fastq.gz | awk '{print $5}') / $(ls -lh "$TMP"/R2.fastq.gz | awk '{print $5}')"

coverm contig --reference "$CAT" \
    --coupled "$TMP/R1.fastq.gz" "$TMP/R2.fastq.gz" \
    -m tpm covered_fraction -p minimap2-sr -t "$T" \
    --output-file "$SOUT/cov.tsv" 2> "$SOUT/coverm.log"
echo "subsampled_pairs=$NPAIRS seed=$SEED" > "$SOUT/SUBSAMPLED"
touch "$SOUT/DONE"
rm -rf "$TMP"
echo "[$(date -Iseconds)] done $S"
