#!/bin/bash
#SBATCH --job-name=f1_eco_full
#SBATCH --partition=batch,debug
#SBATCH --account=pi-hohndor
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log
#SBATCH --chdir=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/ecology

set -euo pipefail

ROOT_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/ecology"
cd "${ROOT_DIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python"
CANDIDATES="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/candidates.tsv"
REP_TRAITS="/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism/uv_survival_clean/results/rep_traits.tsv"
METADATA="/ibex/scratch/projects/c2014/rob/metagenomics/metadata/sample_metadata.tsv"
PROFILES_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
OUTPUT_DIR="${ROOT_DIR}/results_root_v3"
SCRIPT="${ROOT_DIR}/f1_ecology_root_v3.py"

mkdir -p "${OUTPUT_DIR}"

echo "Starting F1 Ecology Full Census under Slurm Job ${SLURM_JOB_ID}"
echo "Node: $(hostname -s)"
echo "Submit Dir: ${ROOT_DIR}"
echo "Output Dir: ${OUTPUT_DIR}"
echo "Script SHA256: $(sha256sum "${SCRIPT}" | awk '{print $1}')"

# Execute compute via require_slurm_worker.sh guard
./require_slurm_worker.sh "${PYTHON}" "${SCRIPT}" \
    --candidates "${CANDIDATES}" \
    --rep-traits "${REP_TRAITS}" \
    --metadata "${METADATA}" \
    --profiles-dir "${PROFILES_DIR}" \
    --output-dir "${OUTPUT_DIR}"

echo "Full census run finished successfully."
