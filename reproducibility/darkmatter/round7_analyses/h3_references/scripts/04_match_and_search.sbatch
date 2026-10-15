#!/bin/bash
#SBATCH --job-name=h3_search
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/search_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/search_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

echo "=== H3 Matching & Full Census Search Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${THREADS}"
echo "Python: $(${PYTHON} --version 2>&1)"
echo "DIAMOND: $(${DIAMOND} version 2>&1)"

echo -e "\n--- Step 4: Exact Common-Support 1-aa Length Quota Matching ---"
"${PYTHON}" scripts/04_exact_length_matching.py

echo -e "\n--- Step 5: Paired DIAMOND Searches Across ALL 212,084 Queries & Statistical Census ---"
"${PYTHON}" scripts/05_diamond_search_full_census.py "${THREADS}"

echo "=== H3 Matching & Full Census Search Completed at $(date -Is) ==="
