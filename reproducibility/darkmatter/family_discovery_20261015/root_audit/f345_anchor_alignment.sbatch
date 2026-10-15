#!/bin/bash
#SBATCH --job-name=root_f345_anchor_alignment
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f345_anchor_alignment_%j.out
set -euo pipefail
module load hmmer/3.3.2
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/f345_anchor_alignment.py" "$D/root_audit/f345_anchor_alignment_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python "$D/root_audit/f345_anchor_alignment_executed_${SLURM_JOB_ID}.py"
