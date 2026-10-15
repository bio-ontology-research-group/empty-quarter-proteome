#!/bin/bash
#SBATCH --job-name=h4_structure
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/logs/h4_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting Hypothesis 4 Structure Analysis on Ibex"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"

W=/ibex/scratch/projects/c2014/rob/metagenomics
H4_DIR=$W/06_darkmatter/p6_hypotheses_clean/h4_structure
LOG_DIR=$H4_DIR/logs
mkdir -p "$LOG_DIR" "$H4_DIR/results"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
module load foldseek/10-941cd33

echo "[$(date -Iseconds)] Foldseek version:"
foldseek version

echo "[$(date -Iseconds)] Python environment:"
$PY --version

# Step 1: Run full audit of exact hits, disconnected maxima, and control leakage
echo "[$(date -Iseconds)] Step 1: Running exact hits and leakage audit..."
$PY "$H4_DIR/scripts/01_audit_exact_hits_and_leakage.py"

# Step 2: Run deterministic pilot of independent domain segmentation & calibration
echo "[$(date -Iseconds)] Step 2: Running deterministic domain segmentation pilot..."
$PY "$H4_DIR/scripts/02_pilot_domain_segmentation.py"

# Step 3: Checksums and integrity verification
echo "[$(date -Iseconds)] Step 3: Generating sha256 checksums of results..."
cd "$H4_DIR/results"
sha256sum *.tsv *.json > results.sha256.tsv

echo "[$(date -Iseconds)] SUCCESS: Hypothesis 4 computation complete."
