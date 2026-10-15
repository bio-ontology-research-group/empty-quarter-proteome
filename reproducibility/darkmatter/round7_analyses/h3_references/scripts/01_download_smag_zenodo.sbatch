#!/bin/bash
#SBATCH --job-name=h3_smag_dl
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/download_smag_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/download_smag_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

echo "=== H3 Zenodo SMAG Archive Download Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"
echo "Python: $(${PYTHON} --version 2>&1)"

"${PYTHON}" scripts/01_download_smag_zenodo.py "${SLURM_CPUS_PER_TASK:-4}"

echo "=== H3 Zenodo SMAG Archive Download Completed at $(date -Is) ==="
