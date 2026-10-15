#!/bin/bash
#SBATCH --job-name=h2_context_full
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/logs/full_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/logs/full_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context"
cd "$WORKDIR"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
MODE="${1:-full}"
TARGET_PER_ARM="${2:-1461}"
RESULTS_DIR="${3:-results_full}"
UPDATE_TOP_INPUTS="${4:-false}"

echo "=== Hypothesis 2 Neighbourhood Recurrence Pipeline (v3 Protocol - Full Census) ==="
echo "Date: $(date -u)"
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID:-local}"
echo "Mode: ${MODE}, Target per arm: ${TARGET_PER_ARM}, Results dir: ${RESULTS_DIR}"
echo "Python: $("$PYTHON" --version)"
"$PYTHON" -c "import numpy, pandas, scipy, statsmodels; print('Versions: numpy=' + numpy.__version__ + ', pandas=' + pandas.__version__ + ', scipy=' + scipy.__version__ + ', statsmodels=' + statsmodels.__version__)"

mkdir -p "$RESULTS_DIR" logs

echo "--- Step 1: Select Cohort & Resolve Verified R232 Taxonomy ---"
"$PYTHON" scripts/01_select_cohort.py "$RESULTS_DIR" "$MODE" "$TARGET_PER_ARM"

echo "--- Step 2: Extract Neighbourhoods, Verify Contig SHA256+Log, & Hash Inputs ---"
"$PYTHON" scripts/02_extract_neighbourhoods.py "$RESULTS_DIR" 16 "$UPDATE_TOP_INPUTS"

echo "--- Step 3: Analyse Conservation, Positional Reference, Stratified/Matched Inference & Closure Checks ---"
"$PYTHON" scripts/03_analyse_conservation.py "$RESULTS_DIR"

echo "--- Step 4: Record Output File SHA-256 Checksums (excluding self) ---"
(cd "$RESULTS_DIR" && rm -f outputs_sha256.txt && sha256sum $(ls -1 | grep -v '^outputs_sha256\.txt$') > outputs_sha256.txt)

if [ -f "$RESULTS_DIR/executed_inputs_manifest.tsv" ]; then
    cp -f "$RESULTS_DIR/executed_inputs_manifest.tsv" INPUTS.tsv
fi

echo "=== Pipeline Complete: $(date -u) ==="
