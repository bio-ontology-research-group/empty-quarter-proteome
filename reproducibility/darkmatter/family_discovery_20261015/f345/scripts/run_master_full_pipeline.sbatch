#!/bin/bash
#SBATCH --job-name=f345_full
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/master_full_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/master_full_%j.err

set -euo pipefail

SCRIPT_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
GUARD="${SCRIPT_DIR}/require_slurm_worker.sh"
PY_ENV="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

# Load modules
module load hmmer/3.3.2

# Add eggnog (DIAMOND) and mmseqs2 env bins to PATH
export PATH="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin:/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin:$PATH"

echo "=== [Master Full Discovery Pipeline F3-F5] Started at $(date -u) ==="
echo "Host: $(hostname -s)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Python: $("${PY_ENV}" --version)"
echo "HMMER: $(hmmsearch -h | head -n 2 | tail -n 1)"
echo "DIAMOND: $(diamond version)"
echo "MMseqs2: $(mmseqs version)"

# Step 1: Prepare Full Cohorts (EQ 5,438 + SMAG + GEM)
echo -e "\n>>> STEP 1: PREPARE IMMUTABLE COHORTS <<<"
"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/07_prepare_full_cohorts.py"

# Step 2: Track F3 Microbial & Heliorhodopsins
echo -e "\n>>> STEP 2: TRACK F3 FULL DISCOVERY (RHODOPSINS) <<<"
"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/08_f3_rhodopsin_full.py"

# Step 3: Track F4 Ectoine Synthase EctC & Pathway Operons
echo -e "\n>>> STEP 3: TRACK F4 FULL DISCOVERY (ECTC & OPERONS) <<<"
"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/09_f4_ectc_full.py"

# Step 4: Track F5 Beta- and Gamma-Carbonic Anhydrases & CsoSCA
echo -e "\n>>> STEP 4: TRACK F5 FULL DISCOVERY (CARBONIC ANHYDRASES) <<<"
"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/10_f5_carbonic_anhydrase_full.py"

# Step 5: Cross-Track Audit & Synthesis
echo -e "\n>>> STEP 5: CROSS-TRACK AUDIT & SYNTHESIS <<<"
"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/11_full_audit_and_synthesis.py"

echo -e "\n=== [Master Full Discovery Pipeline F3-F5] Completed Successfully at $(date -u) ==="
