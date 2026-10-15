#!/bin/bash
#SBATCH --job-name=h3_r7_full
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/r7_pipeline_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/r7_pipeline_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
PYRODIGAL="/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
DIAMOND="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

echo "=========================================================================="
echo "=== H3 Round 7 Full Analysis Pipeline Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs allocated: ${THREADS}"
echo "Memory requested: 128G"
echo "Python: $(${PYTHON} --version 2>&1)"
echo "Pyrodigal: $(${PYRODIGAL} --version 2>&1)"
echo "DIAMOND: $(${DIAMOND} version 2>&1)"
echo "=========================================================================="

echo -e "\n[$(date -Is)] >>> STAGE 1: Verify & Download SMAG Zenodo 8223844 Archive Parts (373 parts) <<<"
"${PYTHON}" scripts/01_download_smag_zenodo.py 6

echo -e "\n[$(date -Is)] >>> STAGE 2: Extract & Gene Call 1,500 SMAG Genomes with Pyrodigal v3.7.1 <<<"
"${PYTHON}" scripts/02_extract_and_call_smag.py "${THREADS}"

echo -e "\n[$(date -Is)] >>> STAGE 3: Fetch GEM Nucleotide Assemblies, Call Pyrodigal v3.7.1 & Load Published Proteomes <<<"
"${PYTHON}" scripts/03_fetch_and_call_gem.py 4500000 "${THREADS}"

echo -e "\n[$(date -Is)] >>> STAGE 4: Exact Common-Support 1-aa Length Quota Matching Across All Configurations <<<"
"${PYTHON}" scripts/04_exact_length_matching.py

echo -e "\n[$(date -Is)] >>> STAGE 5: Paired DIAMOND Searches & Statistical Census Across ALL 212,084 Queries <<<"
"${PYTHON}" scripts/05_diamond_search_full_census.py "${THREADS}"

echo -e "\n=========================================================================="
echo "=== H3 Round 7 Full Analysis Pipeline Successfully Completed at $(date -Is) ==="
echo "=========================================================================="
