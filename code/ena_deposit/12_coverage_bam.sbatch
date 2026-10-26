#!/bin/bash
#SBATCH --job-name=ena_asm_cov
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/cov_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/cov_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#
# Per-contig read depth for every assembly, one array task per line of WORK/assemblies/sample_list.tsv (sample,
# contig FASTA, BAM; from 10_build_assembly_table.py --prepare): samtools coverage on the SemiBin2 self-supervised
# BAM (BWA-MEM, the sample's own reads on its own assembly). 10_build_assembly_table.py turns it into the
# length-weighted mean depth = the COVERAGE field. Idempotent: an existing output is skipped.
#
#   N=$(wc -l < /ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/sample_list.tsv)
#   sbatch --array=1-${N}%64 code/ena_deposit/12_coverage_bam.sbatch

set -u
module load samtools/1.16.1
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK}/sample_list.tsv")
[ -n "${line}" ] || { echo "ERROR: no sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
s=$(echo "${line}" | cut -f1)
bam=$(echo "${line}" | cut -f3)
out=${WORK}/coverage/${s}.samtools_coverage.tsv
[ -s "${bam}" ] || { echo "ERROR: ${bam} missing"; exit 1; }
if [ -s "${out}" ]; then echo "[skip] ${out} exists"; exit 0; fi
samtools coverage "${bam}" > "${out}.partial" && mv "${out}.partial" "${out}"
echo "[coverage] ${s}: $(wc -l < "${out}") contigs"
