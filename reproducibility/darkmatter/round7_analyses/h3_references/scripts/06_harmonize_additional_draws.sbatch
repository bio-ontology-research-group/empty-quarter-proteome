#!/bin/bash
#SBATCH --job-name=h3_gem_draws
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/gem_draws_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/gem_draws_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
PYRODIGAL="/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

echo "=== H3 Additional Harmonized Draws Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${THREADS}"
echo "Python: $(${PYTHON} --version 2>&1)"
echo "Pyrodigal: $(${PYRODIGAL} --version 2>&1)"

"${PYTHON}" -u scripts/06_harmonize_additional_draws.py "${THREADS}"

echo "=== H3 Additional Harmonized Draws Completed at $(date -Is) ==="
