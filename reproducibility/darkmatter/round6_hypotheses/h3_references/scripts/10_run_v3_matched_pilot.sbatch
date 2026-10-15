#!/bin/bash
#SBATCH --job-name=h3_v3_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_v3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_v3_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
PYRODIGAL="/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
DIAMOND="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"

echo "=== H3 v3 Matched Pilot Started at $(date -Is) on $(hostname) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}"
echo "Python: $(${PYTHON} --version 2>&1)"
echo "Pyrodigal: $(${PYRODIGAL} --version 2>&1)"
echo "DIAMOND: $(${DIAMOND} version 2>&1)"

echo "--- Step 1: Safe streaming tar extraction & deterministic Pyrodigal gene calling ---"
"${PYTHON}" scripts/08_v3_extract_and_call_smag.py "${SLURM_CPUS_PER_TASK:-16}"

echo "--- Step 2: Canonicalization, 3-seed exact length matching, and paired DIAMOND searches ---"
"${PYTHON}" scripts/09_v3_match_gem_and_search.py "${SLURM_CPUS_PER_TASK:-16}"

echo "=== H3 v3 Matched Pilot Completed at $(date -Is) ==="
