#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=03:00:00
set -euo pipefail
python3 /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/scripts/hash_large_analysis_inputs.py
