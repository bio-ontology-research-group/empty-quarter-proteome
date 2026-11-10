#!/bin/bash
#SBATCH --job-name=exp1_st2_wit
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_wit_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_wit_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
CODE_DIR="${BASE}/code_v2"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic 4-witness rebuild and validation on $(hostname)"

python3 "${CODE_DIR}/validate_four_changed_witnesses.py" \
  --ledger-tsv "${BASE}/results/focused_release_validation/exp1_focused_validation_per_solution_ledger.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-dir "${BASE}/results/focused_release_validation"

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic 4-witness rebuild and validation"
