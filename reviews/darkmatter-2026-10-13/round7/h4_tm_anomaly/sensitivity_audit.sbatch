#!/bin/bash
#SBATCH --job-name=h4_tm_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/tm_anomaly/sensitivity_audit_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/tm_anomaly/sensitivity_audit_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting H4 Foldseek TM-score Anomaly Sensitivity Audit on Ibex"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
$PY -u /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/tm_anomaly/sensitivity_audit.py

echo "[$(date -Iseconds)] Sensitivity audit complete."
