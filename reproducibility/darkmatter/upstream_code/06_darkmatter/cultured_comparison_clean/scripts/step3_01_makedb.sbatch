#!/bin/bash
# step3_01 - DIAMOND database of the 3,024,784 tagged cultured Bakta proteins (all 722 genomes).
#   sbatch scripts/step3_01_makedb.sbatch
#SBATCH --job-name=ale_s3_makedb
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_makedb_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_makedb_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=batch
set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate eggnog
P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean
diamond --version
diamond makedb --in "$P/tmp/cultured_all.faa" -d "$P/tmp/cultured_db" --threads 8
echo "DONE_STEP3_01"
