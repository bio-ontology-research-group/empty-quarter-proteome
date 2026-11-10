#!/bin/bash
#SBATCH --job-name=exp1_st4_def_v5
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_v5_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_v5_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
export PATH="${PATH}:/project/k10018/sw/envs/gapseq/bin"

command -v diamond
diamond --version
python3 -c "import pandas, numpy; print('Python deps verified')"

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1

# Preserve unperformed job 15894447 outputs as explicit namespace-failure baseline
BASELINE_DIR="${BASE}/results/provenance/baseline_15894447_unperformed_namespace_failure"
mkdir -p "${BASELINE_DIR}"
for f in \
  "${BASE}/results/exp1_stage4_deficiency_verification_instances.tsv" \
  "${BASE}/results/exp1_stage4_deficiency_verification_summary.tsv" \
  "${BASE}/results/provenance/stage4_execution_receipts.json" \
  "${BASE}/logs/exp1_st4_def_v4_15894447.out" \
  "${BASE}/logs/exp1_st4_def_v4_15894447.err"; do
  if [ -f "${f}" ]; then
    cp -pn "${f}" "${BASELINE_DIR}/" || true
  fi
done

# Preserve failed job 15894840 logs (inventory path typo before search)
FAILED_15894840_DIR="${BASE}/results/provenance/failed_15894840_inventory_path_typo"
mkdir -p "${FAILED_15894840_DIR}"
for f in \
  "${BASE}/logs/exp1_st4_def_v5_15894840.out" \
  "${BASE}/logs/exp1_st4_def_v5_15894840.err"; do
  if [ -f "${f}" ]; then
    cp -pn "${f}" "${FAILED_15894840_DIR}/" || true
  fi
done

python3 "${BASE}/code/stage4_candidate_deficiency_search_v5.py" \
  --instances-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_instances.tsv" \
  --summary-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_summary.tsv" \
  --inventory-tsv "${BASE}/results/exp1_model_inventory_and_qc_484.tsv" \
  --gapseq-seq-archaea-dir /project/k10018/sw/db/gapseq_seq/Archaea \
  --gapseq-dat-dir /project/k10018/sw/envs/gapseq/share/gapseq/dat \
  --gapseq-user-seq-dir /project/k10018/sw/envs/gapseq/share/gapseq/dat/seq/Archaea/user \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --refined-mags-tar /scratch/project/k10018/eq-science-20261108/inputs/thermo237_refined_mags.tar.gz \
  --assemblies-dir /project/k10018/empty-quarter/data/assemblies \
  --gap-context-tsv-gz /scratch/project/k10018/eq-archaeal-function-20261106/results/gap_context_candidates.tsv.gz \
  --c18-candidates-tsv /project/k10018/work/eq-annot-20261030/results/c18_pathway_holes/candidates.tsv \
  --c18-windows-tsv /project/k10018/work/eq-annot-20261030/results/c18_pathway_holes/windows.tsv \
  --out-dir "${BASE}/results" \
  --tmp-dir "${TMPDIR}" \
  --threads 128 \
  --output-version v5

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed v5 Stage 4 Candidate Deficiency Search"
