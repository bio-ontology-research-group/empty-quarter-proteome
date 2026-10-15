#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
ROOT=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python "$ROOT/scripts/site_level_sensitivity.py" "$ROOT/results/site_level_sensitivity"
