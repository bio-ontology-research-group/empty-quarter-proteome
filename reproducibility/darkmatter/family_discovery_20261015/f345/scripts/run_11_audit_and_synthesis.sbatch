#!/bin/bash
#SBATCH --job-name=f345_audit
#SBATCH --partition=debug,batch
#SBATCH --account=pi-hohndor
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/11_audit_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/11_audit_%j.err

set -euo pipefail

SCRIPT_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
GUARD="${SCRIPT_DIR}/require_slurm_worker.sh"
PY_ENV="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

echo "=== [11 Cross-Track Audit & Synthesis] Started at $(date -u) ==="
echo "Host: $(hostname -s)"
echo "Job ID: ${SLURM_JOB_ID}"

"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/11_full_audit_and_synthesis.py"

echo "=== [11 Cross-Track Audit & Synthesis] Completed at $(date -u) ==="
