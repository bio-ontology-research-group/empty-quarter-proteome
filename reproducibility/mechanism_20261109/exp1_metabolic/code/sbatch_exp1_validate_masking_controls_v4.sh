#!/bin/bash
#SBATCH --job-name=exp1_mask_v4_validate
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_mask_v4_validate_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_mask_v4_validate_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
RESULTS="${BASE}/results"
CODE="${BASE}/code_v5"
EXECUTED="${BASE}/code_v4"
GDIR="/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${RESULTS}" "${BASE}/logs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python3 "${CODE}/validate_masking_controls_v4_outputs.py" \
  --controls-tsv "${RESULTS}/exp1_masking_controls_7hq_v4.tsv" \
  --baselines-tsv "${RESULTS}/exp1_masking_controls_7hq_v4_baselines.tsv" \
  --ledger-tsv "${RESULTS}/exp1_masking_controls_7hq_v4_per_solution_ledger.tsv" \
  --media-json "${RESULTS}/exp1_masking_controls_v4_media_provenance.json" \
  --checkpoint-json "${RESULTS}/exp1_masking_controls_v4_checkpoint.json" \
  --original-report-json "${RESULTS}/exp1_masking_controls_7hq_v4_report.json" \
  --engine-py "${EXECUTED}/exp1_sbml_highs_engine_v4.py" \
  --runner-py "${EXECUTED}/run_masking_controls_7hq_v4_resumable.py" \
  --helper-py "${EXECUTED}/stage2_pilot24_ensemble.py" \
  --pool-pkl "${RESULTS}/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir "${GDIR}" \
  --milp-time-limit-s 15.0 \
  --out-receipt "${RESULTS}/exp1_masking_controls_7hq_v4_validation_receipt.json" \
  --out-dir "${RESULTS}"

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed independent v4 masking post-hoc validation"
