#!/bin/bash
# 04 - profile the two shotgun control libraries against the 5,438 species representatives,
#      with the paper's own CoverM invocation, verbatim.
#
#   sbatch scripts/04_coverm_controls.sbatch
#
# The cohort profiles in
#   rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/<sample>/abund.tsv
# were produced by scripts/run_current_coverm_profile.py (sbatch 193_current_coverm_full_recompute.sbatch)
# with this command line, recorded verbatim in every sample's run_provenance.json:
#
#   coverm genome --reference <prepared/reps.mmi> --minimap2-reference-is-index \
#       --genome-definition <prepared/genome_def.tsv> \
#       --coupled <S>_R1.fastq.gz <S>_R2.fastq.gz -m relative_abundance -p minimap2-sr -t 16
#
# Reused here unchanged, including the prebuilt minimap2 index (-x sr -I 64G over
# 04_functional/refs/reps.fa = the concatenated 5,438 representatives, 2,485,945 contigs) and the
# contig-to-genome definition, so the control profiles are numerically comparable to the 550 soil
# profiles rather than merely similar. Deliberately NOT changed: no --min-covered-fraction, no
# --min-read-percent-identity, no --exclude-supplementary; CoverM 0.7.0 defaults apply, as in the
# cohort run. The _SE.fastq.gz singletons are ignored, as they are for every cohort sample.
#
# Resources copied from 193_current_coverm_full_recompute.sbatch (16 CPUs, 250G, batch). 250G is the
# floor, not padding: the index alone is 61.9 GB.
#SBATCH --job-name=ale_cov_ctrl
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/coverm_ctrl_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/coverm_ctrl_%A_%a.err
set -euo pipefail

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean
PREP=/ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/picrust2_current_20260803/coverm_current/prepared
READS=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads
COVERM=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/coverm

source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate eco_probes

SAMPLES=(SNTC SPositive)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
OUT=$P/results/coverm/$S
mkdir -p "$OUT"

echo "sample: $S"
"$COVERM" --version
echo "index:  $PREP/reps.mmi"
echo "genome-definition: $PREP/genome_def.tsv  ($(wc -l < "$PREP/genome_def.tsv") contigs)"
date

"$COVERM" genome \
    --reference "$PREP/reps.mmi" \
    --minimap2-reference-is-index \
    --genome-definition "$PREP/genome_def.tsv" \
    --coupled "$READS/${S}_R1.fastq.gz" "$READS/${S}_R2.fastq.gz" \
    -m relative_abundance \
    -p minimap2-sr \
    -t "${SLURM_CPUS_PER_TASK:-16}" \
    > "$OUT/abund.tsv" 2> "$OUT/coverm.log"

echo "rows: $(wc -l < "$OUT/abund.tsv")"
echo "unmapped: $(awk -F'\t' '$1=="unmapped"{print $2}' "$OUT/abund.tsv")"
echo "nonzero representatives: $(awk -F'\t' 'NR>1 && $1!="unmapped" && $2+0>0' "$OUT/abund.tsv" | wc -l)"
command grep -E "reads mapped out of|In sample" "$OUT/coverm.log" || true
echo "=== DONE $S ==="; date
