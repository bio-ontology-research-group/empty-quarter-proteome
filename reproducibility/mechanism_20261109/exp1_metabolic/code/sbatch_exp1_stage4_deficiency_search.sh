#!/bin/bash
#SBATCH --job-name=exp1_st4_def
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
CODE_DIR="${BASE}/code_v2"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 4 Candidate-Conditioned Deficiency Search on $(hostname)"

python3 "${CODE_DIR}/stage4_candidate_deficiency_search.py" \
  --instances-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_instances.tsv" \
  --summary-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_summary.tsv" \
  --gapseq-rxn-fasta-dir /project/k10018/sw/db/gapseq_seq/Archaea/rxn \
  --refined-mags-tar /scratch/project/k10018/eq-science-20261108/inputs/thermo237_refined_mags.tar.gz \
  --assemblies-dir /project/k10018/empty-quarter/data/assemblies \
  --gap-context-tsv-gz /scratch/project/k10018/eq-archaeal-function-20261106/results/gap_context_candidates.tsv.gz \
  --c18-candidates-tsv /project/k10018/work/eq-annot-20261030/results/c18_pathway_holes/candidates.tsv \
  --out-dir "${BASE}/results" \
  --tmp-dir "${TMPDIR}" \
  --threads 128

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 4 Candidate-Conditioned Deficiency Search"
