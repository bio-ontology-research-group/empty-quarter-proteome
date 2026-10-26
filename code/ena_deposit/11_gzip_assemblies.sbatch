#!/bin/bash
#SBATCH --job-name=ena_asm_gzip
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/gzip_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/gzip_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#
# One array task per line of WORK/assemblies/sample_list.tsv (sample, contig FASTA, BAM; from
# 10_build_assembly_table.py --prepare): gzip the contig file into WORK/assemblies/fasta/<s>_assembly.fa.gz
# (Webin-CLI uploads this file itself at submit time, so nothing is pre-uploaded) and write
# stats/<s>.tsv = "<n_contigs> <total_bp>" from the source. Short wall time on purpose: tasks take ~5 min and
# backfill into free slots. Idempotent: a complete .fa.gz (written via .partial + mv) and an existing stats file are skipped.
#
#   N=$(wc -l < /ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/sample_list.tsv)
#   sbatch --array=1-${N}%64 code/ena_deposit/11_gzip_assemblies.sbatch

set -u
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK}/sample_list.tsv")
[ -n "${line}" ] || { echo "ERROR: no sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
s=$(echo "${line}" | cut -f1)
fa=$(echo "${line}" | cut -f2)
gz=${WORK}/fasta/${s}_assembly.fa.gz
st=${WORK}/stats/${s}.tsv
[ -s "${fa}" ] || { echo "ERROR: ${fa} missing"; exit 1; }

if [ -s "${gz}" ]; then
    echo "[skip] ${gz} exists"
else
    if command -v pigz >/dev/null; then pigz -p "${SLURM_CPUS_PER_TASK}" -c "${fa}" > "${gz}.partial"
    else gzip -c "${fa}" > "${gz}.partial"; fi
    mv "${gz}.partial" "${gz}"
    echo "[gzip] ${s} -> $(du -h "${gz}" | cut -f1)"
fi
if [ -s "${st}" ]; then
    echo "[skip] ${st} exists"
else
    awk 'BEGIN{n=0;bp=0} /^>/{n++; next} {bp+=length($0)} END{print n, bp}' "${fa}" > "${st}.partial"
    mv "${st}.partial" "${st}"
    echo "[stats] ${s}: $(cat "${st}")"
fi
