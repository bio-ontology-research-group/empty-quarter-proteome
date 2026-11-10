#!/bin/bash
#SBATCH --job-name=eq-hs3-map
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/heat_stress_20261016
/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/python scripts/hs3_map.py
