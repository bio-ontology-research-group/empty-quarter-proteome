#!/bin/bash
#SBATCH --job-name=f3_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/03_f3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/03_f3_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

echo "=== [Track F3 Rhodopsin Pilot] Started at $(date -u) ==="
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID}"
./require_slurm_worker.sh "$PYTHON" scripts/03_f3_rhodopsin_pilot.py
echo "=== Track F3 Pilot Finished at $(date -u) ==="
