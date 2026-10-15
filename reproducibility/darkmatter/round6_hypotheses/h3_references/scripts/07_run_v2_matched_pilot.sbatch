#!/bin/bash
#SBATCH --job-name=h3_v2_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_v2_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_v2_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

echo "=== H3 v2 Matched Pilot Started at $(date) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}"

# Step 1: Download public SMAG parts from Zenodo (if not already fetched)
echo "--- Step 1: Fetching public SMAG archive parts ---"
"${PYTHON}" scripts/04_fetch_smag_parts.py 2

# Step 2: Extract complete genomes and call genes with Pyrodigal v3.7.1
echo "--- Step 2: Extracting genomes and calling genes with Pyrodigal ---"
"${PYTHON}" scripts/05_extract_and_genecall_smag.py 8

# Step 3: Match GEM subset, audit redundancy/overlap, and run paired DIAMOND search
echo "--- Step 3: Matching GEM subset, auditing overlap, and running DIAMOND search ---"
"${PYTHON}" scripts/06_match_gem_and_search.py 8

echo "=== H3 v2 Matched Pilot Completed at $(date) ==="
