#!/bin/bash
#SBATCH --job-name=root_family_closure
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/family_closure_%j.out
set -euo pipefail
module load hmmer/3.3.2
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
cp "$D/root_audit/finalize_family_bundle.py" "$D/root_audit/finalize_family_bundle_executed_${SLURM_JOB_ID}.py"
"$D/f12/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python "$D/root_audit/finalize_family_bundle_executed_${SLURM_JOB_ID}.py"
