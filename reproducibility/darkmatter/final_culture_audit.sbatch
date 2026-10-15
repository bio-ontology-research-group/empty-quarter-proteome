#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" "$R/scripts/culture_locus_audit.py"
"$PY" "$R/scripts/culture_family_membership.py"
