#!/bin/bash
#SBATCH --job-name=h4_v2_crop
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_v2_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_v2_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting Hypothesis 4 V2 Cropped Region Analysis on Ibex"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"

W=/ibex/scratch/projects/c2014/rob/metagenomics
H4_DIR=$W/06_darkmatter/p6_hypotheses_clean/h4_structure
mkdir -p "$H4_DIR/logs" "$H4_DIR/results"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
module load foldseek/10-941cd33

echo "[$(date -Iseconds)] Foldseek version:"
foldseek version

echo "[$(date -Iseconds)] Python environment:"
$PY --version

# Run V2 Cropped Region Analysis (Deterministic Pilot, N=300 models)
echo "[$(date -Iseconds)] Executing v2_crop_and_search_regions.py..."
$PY "$H4_DIR/scripts/v2_crop_and_search_regions.py" --mode pilot --threads 16

# Checksums and verification
echo "[$(date -Iseconds)] Generating sha256 checksums of v2 outputs..."
cd "$H4_DIR/results"
sha256sum v2_*.tsv v2_*.json > v2_results.sha256.tsv

echo "[$(date -Iseconds)] SUCCESS: Hypothesis 4 V2 computation complete."
