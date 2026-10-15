#!/bin/bash
#SBATCH --job-name=root_f345_plots
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f345_plots_%j.out
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/f345_embedding_analysis.py" "$D/root_audit/f345_embedding_analysis_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python "$D/root_audit/f345_embedding_analysis_executed_${SLURM_JOB_ID}.py"
