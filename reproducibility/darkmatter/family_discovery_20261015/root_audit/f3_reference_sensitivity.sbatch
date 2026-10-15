#!/bin/bash
#SBATCH --job-name=root_f3_ref_sensitivity
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f3_reference_sensitivity_%j.out
set -euo pipefail
module load hmmer/3.3.2
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/f3_reference_sensitivity.py" "$D/root_audit/f3_reference_sensitivity_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python "$D/root_audit/f3_reference_sensitivity_executed_${SLURM_JOB_ID}.py"
