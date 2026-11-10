#!/bin/bash
#SBATCH --job-name=exp1_st1_qc
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st1_qc_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st1_qc_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs" "${BASE}/provenance"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 1 (538/484 Inventory, Biomass Template Audit, QC, Pilot24 Manifest) on $(hostname)"
python3 -m pytest -v "${BASE}/tests/test_exp1_metabolic.py"

python3 "${BASE}/code/stage1_inventory_and_qc.py" \
  --science-inputs /scratch/project/k10018/eq-science-20261108/inputs \
  --external-manifest-tsv "${BASE}/inputs/external_signalp_manifest.tsv" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --excluded-rxns-tsv /project/k10018/sw/envs/gapseq/share/gapseq/dat/biomass/excluded_reactions.tsv \
  --seed-reactions-tsv /project/k10018/sw/envs/gapseq/share/gapseq/dat/seed_reactions_corrected.tsv \
  --seed-metabolites-tsv /project/k10018/sw/envs/gapseq/share/gapseq/dat/seed_metabolites_edited.tsv \
  --out-dir "${BASE}/results" \
  --workers 96

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 1"
