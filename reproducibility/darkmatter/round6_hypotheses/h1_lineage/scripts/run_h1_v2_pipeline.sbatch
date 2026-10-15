#!/bin/bash
#SBATCH --job-name=h1_v2_ancestry
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/logs/h1_v2_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/logs/h1_v2_%j.err

set -euo pipefail

echo "========================================================================"
echo "Starting H1 Lineage V2 Ancestry-Aware Analysis on IBEX"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "========================================================================"

WORKDIR=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python

mkdir -p "${WORKDIR}/results" "${WORKDIR}/logs"

$PY "${WORKDIR}/scripts/run_h1_v2_pipeline.py"

echo "========================================================================"
echo "H1 V2 Pipeline Finished: $(date)"
echo "========================================================================"
