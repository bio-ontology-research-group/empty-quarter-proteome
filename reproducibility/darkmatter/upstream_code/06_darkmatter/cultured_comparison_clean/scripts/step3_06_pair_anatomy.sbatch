#!/bin/bash
# step3_06 - anatomy of the 23 same-species MAG / cultured pairs.
# One node does all four stages: the FASTA preparation, the parallel byte-range scan of
# the 36 GiB all_cluster90.tsv that gives every MAG protein its dark flag, the 2 x 23 (+
# 2 x 14 pooled) small DIAMOND runs, and the tables. Re-running skips finished stages.
#   sbatch scripts/step3_06_pair_anatomy.sbatch
#SBATCH --job-name=ale_s3_anat
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_06_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_06_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --partition=batch
set -euo pipefail

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean
PY=/home/velazqal/miniconda3/envs/python-utils/bin/python

echo "[$(date -Iseconds)] start on $(hostname), ${SLURM_CPUS_PER_TASK} cpus"
"$PY" "$P/scripts/step3_06_pair_anatomy.py" --stage all --workers "${SLURM_CPUS_PER_TASK:-8}"
echo "[$(date -Iseconds)] done"
