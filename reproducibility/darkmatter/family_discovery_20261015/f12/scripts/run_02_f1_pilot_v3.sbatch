#!/bin/bash
#SBATCH --job-name=f12_02_f1v3
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/02_f1v3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/02_f1v3_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

./require_slurm_worker.sh "$PYTHON" scripts/02_f1_pilot_v3.py
