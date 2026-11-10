#!/bin/bash
#SBATCH --job-name=exp1_st2_foc
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_foc_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_foc_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
CODE_DIR="${BASE}/code_v2"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results/focused_release_validation" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Focused Release Validation on $(hostname)"

python3 "${CODE_DIR}/stage2_focused_release_validation.py" \
  --old-summary-tsv "${BASE}/outputs/pilot_15873653_archive/exp1_pilot24_ensemble_summary.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-dir "${BASE}/results/focused_release_validation" \
  --milp-time-limit-s 10.0 \
  --workers 32

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Focused Release Validation (PASSED)"
