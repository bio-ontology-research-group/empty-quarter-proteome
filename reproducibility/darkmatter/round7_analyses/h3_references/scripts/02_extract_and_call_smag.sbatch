#!/bin/bash
#SBATCH --job-name=h3_smag_ext
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/extract_smag_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/extract_smag_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
PYRODIGAL="/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

echo "=== H3 SMAG Extraction & Gene Calling Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${THREADS}"
echo "Python: $(${PYTHON} --version 2>&1)"
echo "Pyrodigal: $(${PYRODIGAL} --version 2>&1)"

"${PYTHON}" scripts/02_extract_and_call_smag.py "${THREADS}"

echo "=== H3 SMAG Extraction & Gene Calling Completed at $(date -Is) ==="
