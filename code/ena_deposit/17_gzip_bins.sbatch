#!/bin/bash
#SBATCH --job-name=ena_bin_gzip
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/bins/logs/gzip_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/bins/logs/gzip_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#
# gzip every bin of wgs_bins.tsv (column 1 = bin name; FASTA = BINS/<bin>.fa) into WORK/bins/fasta/<bin>.fa.gz, the
# file Webin-CLI uploads at submit time. One array task per 500 bins (25 tasks cover 12,434); idempotent (complete
# .fa.gz skipped, written via .partial + mv).
#   sbatch --array=1-25 code/ena_deposit/17_gzip_bins.sbatch
set -u
TABLE=/ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-proteome/code/ena_deposit/wgs_bins.tsv
BINS=/ibex/scratch/projects/c2014/rob/metagenomics/05_refine/_refined_bins_all
OUT=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/bins/fasta
mkdir -p "$OUT"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * 500 + 2 ))     # +2: skip the header, 1-based lines
end=$(( start + 499 ))
n=0
while read -r bin; do
    gz=$OUT/$bin.fa.gz
    fasta=$BINS/$bin.fa
    [ -s "$gz" ] && continue
    [ -s "$fasta" ] || { echo "ERROR: $fasta missing"; exit 1; }
    gzip -c "$fasta" > "$gz.partial" && mv "$gz.partial" "$gz"
    n=$((n+1))
done < <(sed -n "${start},${end}p" "$TABLE" | cut -f1)
echo "task ${SLURM_ARRAY_TASK_ID}: compressed $n"
