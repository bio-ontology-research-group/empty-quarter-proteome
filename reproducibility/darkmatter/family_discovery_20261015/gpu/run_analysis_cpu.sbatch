#!/bin/bash
#SBATCH --job-name=eq_family_analysis_cpu
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/analysis_cpu_%j.out
#SBATCH --error=logs/analysis_cpu_%j.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/home/hohndor/eq-family-discovery-20261015/cpu}"
cd "${ROOT_DIR}"

mkdir -p logs

echo "=== IBEX CPU Family Embedding Analysis Job ==="
echo "Date:          $(date -u +%FT%TZ)"
echo "Host:          $(hostname -s)"
echo "Slurm Job ID:  ${SLURM_JOB_ID:-none}"
echo "Partition:     ${SLURM_JOB_PARTITION:-none}"
echo "Working Dir:   ${ROOT_DIR}"
echo "==============================================="

# 1. Enforce worker node execution guard
./require_slurm_worker.sh echo "Slurm worker execution verified on $(hostname -s)"

# 2. Select verified python environment with numpy, scipy, sklearn, matplotlib
PYTHON_EXEC="${PYTHON_BIN:-/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python}"
if [[ ! -x "${PYTHON_EXEC}" ]]; then
    echo "ERROR: Specified Python executable not found or not executable: ${PYTHON_EXEC}" >&2
    exit 1
fi

echo "Using Python: ${PYTHON_EXEC}"
"${PYTHON_EXEC}" -c "import numpy, scipy, sklearn, matplotlib; print(f'Libraries verified: numpy={numpy.__version__}, scipy={scipy.__version__}, sklearn={sklearn.__version__}, matplotlib={matplotlib.__version__}')"

# 3. Parameters and inputs
EMBEDDINGS="${EMBEDDINGS:-${ROOT_DIR}/embeddings.npz}"
METADATA="${METADATA:-${ROOT_DIR}/smoke_test_metadata.tsv}"
OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/output}"

if [[ ! -f "${EMBEDDINGS}" ]]; then
    echo "ERROR: Embeddings file not found: ${EMBEDDINGS}" >&2
    exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: Metadata file not found: ${METADATA}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Record job environment receipt
cat <<EOF > "${OUTPUT_DIR}/slurm_execution.txt"
slurm_job_id=${SLURM_JOB_ID}
slurm_job_name=${SLURM_JOB_NAME:-eq_family_analysis_cpu}
slurm_partition=${SLURM_JOB_PARTITION:-batch}
slurm_nodelist=${SLURM_JOB_NODELIST}
hostname=$(hostname)
date_utc=$(date -u +%FT%TZ)
python_executable=${PYTHON_EXEC}
embeddings_input=${EMBEDDINGS}
metadata_input=${METADATA}
EOF

# 4. Execute analysis script through Slurm worker guard
./require_slurm_worker.sh "${PYTHON_EXEC}" "${ROOT_DIR}/family_embedding_analysis.py" \
    --embeddings "${EMBEDDINGS}" \
    --metadata "${METADATA}" \
    --output-dir "${OUTPUT_DIR}" \
    --seeds "42,123,20261015" \
    --perplexities "5,30,50" \
    --tsne-max-iter 1000

# 5. Generate checksums for output artifacts
(
    cd "${OUTPUT_DIR}"
    sha256sum * > checksums.sha256 2>/dev/null || true
    touch COMPLETE
)

echo "=== Analysis completed successfully at $(date -u +%FT%TZ) ==="
