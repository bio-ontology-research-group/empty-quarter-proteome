#!/bin/bash
#SBATCH --job-name=eq-hs3-select
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/heat_stress_20261016
/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/python scripts/hs3_select.py
