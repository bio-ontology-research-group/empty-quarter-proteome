#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
BASE=/ibex/scratch/projects/c2014/rob/metagenomics
python3 "$BASE/06_darkmatter/p4_review_20261013/scripts/selection_covariates.py" "$BASE/06_darkmatter/darkfam_basis_clean/results/selection" "$BASE/06_darkmatter/p4_review_20261013/results/selection_covariates"
