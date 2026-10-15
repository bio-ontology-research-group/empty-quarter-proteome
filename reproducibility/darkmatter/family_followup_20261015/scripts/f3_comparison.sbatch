#!/bin/bash
#SBATCH --job-name=eq-helio-controls
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python scripts/f3_comparison.py
