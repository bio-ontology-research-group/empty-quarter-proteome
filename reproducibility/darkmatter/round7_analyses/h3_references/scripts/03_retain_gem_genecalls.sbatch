#!/bin/bash
#SBATCH --job-name=h3_gem_retention
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/gem_retention_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/gem_retention_%j.err
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python -u scripts/03_retain_gem_genecalls.py "${SLURM_CPUS_PER_TASK}"
