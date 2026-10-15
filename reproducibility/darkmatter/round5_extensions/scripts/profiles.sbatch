#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/logs/round5-profiles-%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=8
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/scripts/round5/profile_extensions.py
