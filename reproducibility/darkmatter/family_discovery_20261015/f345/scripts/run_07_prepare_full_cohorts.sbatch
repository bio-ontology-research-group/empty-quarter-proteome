#!/bin/bash
#SBATCH --job-name=f345_cohorts
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/07_cohorts_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/07_cohorts_%j.err

set -euo pipefail

SCRIPT_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
GUARD="${SCRIPT_DIR}/require_slurm_worker.sh"
PY_ENV="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

echo "=== [07 Prepare Full Cohorts] Started at $(date -u) ==="
echo "Host: $(hostname -s)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/07_prepare_full_cohorts.py"

echo "=== [07 Prepare Full Cohorts] Completed at $(date -u) ==="
