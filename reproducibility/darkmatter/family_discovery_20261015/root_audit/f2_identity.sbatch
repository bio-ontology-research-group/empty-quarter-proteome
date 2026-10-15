#!/bin/bash
#SBATCH --job-name=root_f2_identity
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f2_identity_%j.out
set -euo pipefail
module load hmmer/3.3.2
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/f2_identity.py" "$D/root_audit/f2_identity_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python "$D/root_audit/f2_identity_executed_${SLURM_JOB_ID}.py"
