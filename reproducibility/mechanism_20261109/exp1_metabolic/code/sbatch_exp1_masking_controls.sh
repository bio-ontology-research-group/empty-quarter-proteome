#!/bin/bash
#SBATCH --job-name=exp1_st2_mask7hq
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_mask7hq_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_mask7hq_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic 7 HQ Masking Controls on $(hostname)"

python3 "${BASE}/code_v3/run_masking_controls_7hq.py" \
  --manifest-tsv "${BASE}/results/exp1_pilot24_manifest.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-tsv "${BASE}/results/exp1_masking_controls_7hq.tsv" \
  --report-json "${BASE}/results/exp1_masking_controls_7hq_report.json" \
  --workers 32 \
  --milp-time-limit-s 15.0

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic 7 HQ Masking Controls (PASSED)"
