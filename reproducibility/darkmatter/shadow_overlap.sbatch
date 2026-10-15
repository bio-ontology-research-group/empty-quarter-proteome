#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=03:00:00
set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python "$R/scripts/shadow_overlap.py" "$R/results/shadow_overlap"
