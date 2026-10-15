#!/bin/bash
#SBATCH --job-name=h4_v7_full
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/h4_v7_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/h4_v7_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting Hypothesis 4 V7 Full-Cohort (N=5,851) Analysis on Ibex"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"

W=/ibex/scratch/projects/c2014/rob/metagenomics
H4_DIR=$W/06_darkmatter/p7_analyses_clean/h4_structure
mkdir -p "$H4_DIR/logs" "$H4_DIR/results"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
module load foldseek/10-941cd33

echo "[$(date -Iseconds)] Foldseek version:"
foldseek version

echo "[$(date -Iseconds)] Python environment:"
$PY --version

# Run V7 Full-Cohort Cropped & Whole-Chain Analysis (N=5,851 models, --max-seqs 1000)
echo "[$(date -Iseconds)] Executing v7_crop_and_search_regions.py..."
$PY "$H4_DIR/scripts/v7_crop_and_search_regions.py" --mode full --threads "${SLURM_CPUS_PER_TASK}" --max-seqs 1000

# Generate sha256 checksums of v7 outputs
echo "[$(date -Iseconds)] Generating sha256 checksums of v7 outputs..."
cd "$H4_DIR/results"
sha256sum v7_*.tsv v7_*.json > v7_results.sha256.tsv

# Run independent verification closure script
echo "[$(date -Iseconds)] Running independent closure verification..."
$PY "$H4_DIR/scripts/v7_verify_closure.py"

echo "[$(date -Iseconds)] SUCCESS: Hypothesis 4 V7 full-cohort computation complete."
