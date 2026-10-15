#!/bin/bash
#SBATCH --job-name=h2_context_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/pilot_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/pilot_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context"
cd "$WORKDIR"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

echo "=== Hypothesis 2 Feasibility Pilot ==="
echo "Date: $(date -u)"
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID:-local}"
echo "Python: $("$PYTHON" --version)"

mkdir -p results logs

echo "--- Step 1: Select Pilot Cohort ---"
"$PYTHON" scripts/01_select_cohort.py results

echo "--- Step 2: Extract Neighbourhoods & Join Annotations ---"
"$PYTHON" scripts/02_extract_neighbourhoods.py results

echo "--- Step 3: Analyse Conservation & Permutation Null ---"
"$PYTHON" scripts/03_analyse_conservation.py results

echo "=== Pilot Complete: $(date -u) ==="
