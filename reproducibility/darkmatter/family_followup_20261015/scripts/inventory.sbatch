#!/bin/bash
#SBATCH --job-name=eq-family-inventory
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
mkdir -p inputs results
find /ibex/scratch/projects/c2014/rob/metagenomics -maxdepth 5 \( -iname '*bam' -o -iname '*metadata*' -o -iname '*taxonom*' -o -iname '*assembly*' -o -iname '*checkm*' \) > results/available_paths.txt
curl -fL --retry 2 --max-time 120 https://raw.githubusercontent.com/Caiyulu-818/SMAG/main/Supplementary/Extended%20Data%20Table%202.xlsx -o inputs/smag_extended_table2.xlsx
curl -fL --retry 2 --max-time 120 https://raw.githubusercontent.com/Caiyulu-818/SMAG/main/Supplementary/SMAG_README.md -o inputs/smag_readme.md
curl -fL --retry 2 --max-time 120 https://portal.nersc.gov/GEM/ -o inputs/gem_index.html
sha256sum inputs/* > results/download_sha256.txt
find /ibex/scratch/projects/c2014/rob/miniconda3/envs -maxdepth 3 \( -name '*tmhmm*' -o -name '*phobius*' -o -name '*samtools*' -o -name '*mafft*' \) > results/available_tools.txt
echo "$SLURM_JOB_ID" > results/inventory_complete.txt
