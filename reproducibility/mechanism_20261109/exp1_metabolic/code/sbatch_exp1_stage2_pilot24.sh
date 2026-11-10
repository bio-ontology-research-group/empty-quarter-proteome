#!/bin/bash
#SBATCH --job-name=exp1_st2_pilot24
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_pilot24_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st2_pilot24_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs" "${BASE}/provenance"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 2 (Pilot24 + HQ8 Masking Controls & Scientific QC Gate) on $(hostname)"

python3 "${BASE}/code/stage2_pilot24_ensemble.py" \
  --manifest-tsv "${BASE}/results/exp1_pilot24_manifest.tsv" \
  --pool-pkl "${BASE}/results/exp1_candidate_reaction_pool.pkl" \
  --gapseq-out-dir /project/k10018/work/eq-annot-20261030/results/c11_gapseq/out \
  --supp-models-dir /scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models \
  --out-dir "${BASE}/results" \
  --n-target-solutions 20 \
  --milp-time-limit-s 10.0 \
  --workers 32

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 2 (Pilot24 Scientific QC Gate PASSED)"
