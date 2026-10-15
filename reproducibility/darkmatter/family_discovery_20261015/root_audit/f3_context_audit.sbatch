#!/bin/bash
#SBATCH --job-name=root_f3_context_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f3_context_audit_%j.out
set -euo pipefail
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/f3_context_audit.py" "$D/root_audit/f3_context_audit_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python "$D/root_audit/f3_context_audit_executed_${SLURM_JOB_ID}.py"
