#!/bin/bash
#SBATCH --job-name=exp1_st2_mask7v4
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_mask7hq_v4_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_mask7hq_v4_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
RESULTS="${BASE}/results"
LOGS="${BASE}/logs"
CODE="${BASE}/code_v4"
GDIR="/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"
SUPP="/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"

export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${RESULTS}" "${LOGS}" "${CODE}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic 7 HQ Masking Controls (v4 Resumable) on ${SLURMD_NODENAME:-$(hostname)}"
echo "Code directory: ${CODE}"

python3 "${CODE}/run_masking_controls_7hq_v4_resumable.py" \
    --manifest-tsv "${RESULTS}/exp1_pilot24_manifest.tsv" \
    --pool-pkl "${RESULTS}/exp1_candidate_reaction_pool.pkl" \
    --gapseq-out-dir "${GDIR}" \
    --supp-models-dir "${SUPP}" \
    --out-tsv "${RESULTS}/exp1_masking_controls_7hq_v4.tsv" \
    --baselines-tsv "${RESULTS}/exp1_masking_controls_7hq_v4_baselines.tsv" \
    --ledger-tsv "${RESULTS}/exp1_masking_controls_7hq_v4_per_solution_ledger.tsv" \
    --media-json "${RESULTS}/exp1_masking_controls_v4_media_provenance.json" \
    --checkpoint-json "${RESULTS}/exp1_masking_controls_v4_checkpoint.json" \
    --report-json "${RESULTS}/exp1_masking_controls_7hq_v4_report.json" \
    --workers 32 \
    --milp-time-limit-s 15.0

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic 7 HQ Masking Controls (v4 Resumable)"
