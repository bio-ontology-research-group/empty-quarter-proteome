#!/bin/bash
# 08 - how many blank reads are the E. coli MAG? Step 04 is the paper's invocation verbatim, and the
# paper's invocation reports relative abundance, which is coverage-weighted and so does not say how
# many reads a genome actually took. For the negative control that number is the claim, so this
# reruns the identical mapping with CoverM's count and covered_fraction methods added.
#
#   sbatch scripts/08_coverm_sntc_counts.sbatch
#
# SNTC only: it is 335,190 reads and costs two minutes. Reference, index, genome definition and
# mapper are unchanged from step 04. Two DEVIATIONS, both forced and both in one direction:
#   -m count instead of relative_abundance, and --min-covered-fraction 0 instead of the CoverM 0.7.0
#   default of 10. CoverM refuses `count` while a breadth gate is set, so the two cannot be run
#   together. That gate is the one the exploratory check listed as an unverified caveat: the cohort
#   profiles, and step 04, report zero for any genome covered over less than 10 % of its length.
#   Running counts ungated therefore also measures what that gate hides in the blank.
#SBATCH --job-name=ale_cov_sntc_n
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/coverm_sntc_counts_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/coverm_sntc_counts_%j.err
set -euo pipefail

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean
PREP=/ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/picrust2_current_20260803/coverm_current/prepared
READS=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads
COVERM=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/coverm

source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate eco_probes
OUT=$P/results/coverm/SNTC_counts
mkdir -p "$OUT"

"$COVERM" genome \
    --reference "$PREP/reps.mmi" \
    --minimap2-reference-is-index \
    --genome-definition "$PREP/genome_def.tsv" \
    --coupled "$READS/SNTC_R1.fastq.gz" "$READS/SNTC_R2.fastq.gz" \
    -m count \
    --min-covered-fraction 0 \
    -p minimap2-sr \
    -t "${SLURM_CPUS_PER_TASK:-16}" \
    > "$OUT/counts.tsv" 2> "$OUT/coverm.log"

echo "--- genomes with at least one read ---"
head -1 "$OUT/counts.tsv"
awk -F'\t' 'NR>1 && $2+0>0' "$OUT/counts.tsv" | sort -t$'\t' -k2,2gr
echo "genomes with >=1 read: $(awk -F'\t' 'NR>1 && $1!="unmapped" && $2+0>0' "$OUT/counts.tsv" | wc -l)"
echo "=== DONE ==="
