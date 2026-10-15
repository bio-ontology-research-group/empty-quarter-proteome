#!/bin/bash
#SBATCH --job-name=h3_raw_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/final_raw_audit/h3_raw_audit_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/final_raw_audit/h3_raw_audit_%j.err

set -euo pipefail

echo "Starting H3 Raw-Evidence Audit on Slurm allocation..."
echo "Date: $(date -u)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"

/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/require_slurm_worker.sh \
  /ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python \
  /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/final_raw_audit/audit_h3_final_raw.py

echo "Audit completed successfully: $(date -u)"
