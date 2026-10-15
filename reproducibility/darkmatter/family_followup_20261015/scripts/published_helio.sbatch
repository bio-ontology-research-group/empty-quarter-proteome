#!/bin/bash
#SBATCH --job-name=eq-helio-published-check
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
module load hmmer/3.3.2
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python scripts/published_helio.py
