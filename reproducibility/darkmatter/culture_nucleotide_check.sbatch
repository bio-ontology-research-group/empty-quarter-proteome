#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
BASE=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013
python3 "$BASE/scripts/culture_nucleotide_check.py" "$BASE/results/culture_nucleotide_cigar"
