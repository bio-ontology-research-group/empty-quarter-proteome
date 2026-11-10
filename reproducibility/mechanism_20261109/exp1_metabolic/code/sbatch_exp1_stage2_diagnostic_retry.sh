#!/bin/bash
#SBATCH --job-name=exp1_st2_diag_retry
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_diag_retry_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_diag_retry_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs" "${BASE}/outputs/diagnostic_retry"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 2 Focused Diagnostic Retry on $(hostname)"

python3 "${BASE}/code/stage2_focused_diagnostic_retry.py" \
  --old-summary-tsv "${BASE}/provenance/run_15873653/exp1_pilot24_ensemble_summary.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-dir "${BASE}/outputs/diagnostic_retry" \
  --milp-time-limit-s 12.0 \
  --workers 32

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 2 Focused Diagnostic Retry"
