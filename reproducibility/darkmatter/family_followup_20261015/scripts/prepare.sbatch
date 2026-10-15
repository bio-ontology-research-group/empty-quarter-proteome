#!/bin/bash
#SBATCH --job-name=eq-family-prepare
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
cp scripts/prepare.py results/prepare_executed_${SLURM_JOB_ID}.py
/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python scripts/prepare.py
