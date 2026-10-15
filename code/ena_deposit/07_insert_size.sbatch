#!/bin/bash
#SBATCH --job-name=ena_wgs_insert
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/insert_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/insert_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#
# Measured insert size per sample from the binning BAMs (each sample's reads mapped to its own assembly,
# wgs_assemblies/binning/<sample>/<sample>_metaG.reads.sorted.bam). One array task per canonical sample;
# samtools stats "insert size average/sd". Appends to WORK/logs/insert_sizes.tsv (read by 00_build_runs_table.py).
#
# Usage:  sbatch --array=1-551 code/ena_deposit/07_insert_size.sbatch
set -u
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs
BAMDIR=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/binning
COHORT=/ibex/scratch/projects/c2014/rob/metagenomics/metadata/sample_metadata.tsv
OUT=${WORK}/logs/insert_sizes.tsv
module load samtools/1.16.1
[ -f "${OUT}" ] || printf 'sample\tinsert_mean\tinsert_sd\tpairs\n' > "${OUT}"
sample=$(awk -F'\t' -v n=$(( SLURM_ARRAY_TASK_ID + 1 )) 'NR==n{print $1}' "${COHORT}")
[ -n "${sample}" ] || { echo "no sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
grep -q "^${sample}	" "${OUT}" && { echo "[skip] ${sample} done"; exit 0; }
bam=${BAMDIR}/${sample}/${sample}_metaG.reads.sorted.bam
[ -f "${bam}" ] || { echo "no BAM for ${sample}"; printf '%s\t\t\t0\n' "${sample}" >> "${OUT}"; exit 0; }
stats=$(samtools stats -@ 2 "${bam}" | grep '^SN')
mean=$(echo "${stats}" | awk -F'\t' '$2=="insert size average:"{print $3}')
sd=$(echo "${stats}" | awk -F'\t' '$2=="insert size standard deviation:"{print $3}')
pairs=$(echo "${stats}" | awk -F'\t' '$2=="reads properly paired:"{print $3}')
printf '%s\t%s\t%s\t%s\n' "${sample}" "${mean}" "${sd}" "${pairs}" >> "${OUT}"
echo "${sample} insert ${mean} +/- ${sd} (${pairs} properly paired reads)"
