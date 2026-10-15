#!/bin/bash
#SBATCH --job-name=h4_v3_full
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_v3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_v3_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting Hypothesis 4 V3 Full-Cohort (N=5,851) Analysis on Ibex"
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

# Run V3 Full-Cohort Cropped & Whole-Chain Analysis (N=5,851 models, --max-seqs 1000)
echo "[$(date -Iseconds)] Executing v3_crop_and_search_regions.py..."
$PY "$H4_DIR/scripts/v3_crop_and_search_regions.py" --mode full --threads 16 --max-seqs 1000

# Remove temporary Foldseek scratch directories while preserving v3_crops, v3_wc_queries, and raw hits
rm -rf "$H4_DIR/v3_fs_wc_tmp" "$H4_DIR/v3_fs_crop_tmp"

# Checksums and verification
echo "[$(date -Iseconds)] Generating sha256 checksums of v3 outputs..."
cd "$H4_DIR/results"
sha256sum v3_*.tsv v3_*.json > v3_results.sha256.tsv

echo "[$(date -Iseconds)] SUCCESS: Hypothesis 4 V3 full-cohort computation complete."
