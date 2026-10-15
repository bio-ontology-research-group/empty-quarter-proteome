#!/bin/bash
#SBATCH --job-name=h1_lineage
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/logs/h1_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/logs/h1_%j.err

set -euo pipefail

echo "========================================================================"
echo "Starting H1 Lineage Feasibility Pilot on IBEX"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "========================================================================"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
WORKDIR=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage

mkdir -p "${WORKDIR}/results" "${WORKDIR}/logs"

$PY "${WORKDIR}/scripts/run_h1_pipeline.py"

echo "========================================================================"
echo "H1 Pipeline Finished: $(date)"
echo "========================================================================"
