#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/logs/round5-culture-%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/scripts/round5/culture_caller.py
