#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
W=/ibex/scratch/projects/c2014/rob/metagenomics
R=$W/06_darkmatter/p4_review_20261013
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" "$R/scripts/selection_estimands.py" "$W" "$R/results/selection_estimands"
"$PY" "$R/scripts/structural_raw_audit.py" "$W/04_structure/00_structure_clean/05_validate" "$R/results/structural_raw_audit"
