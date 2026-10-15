#!/bin/bash
#SBATCH --job-name=h2_context_v3
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/pilot_v3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/pilot_v3_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context"
cd "$WORKDIR"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
MODE="${1:-pilot}"
TARGET_PER_ARM="${2:-150}"
RESULTS_DIR="${3:-results}"
UPDATE_TOP_INPUTS="${4:-true}"

echo "=== Hypothesis 2 Neighbourhood Recurrence Pipeline (v3 Protocol) ==="
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

echo "--- Step 4: Record Output File SHA-256 Checksums ---"
sha256sum "$RESULTS_DIR"/* > "$RESULTS_DIR/outputs_sha256.txt"

echo "=== Pipeline Complete: $(date -u) ==="
