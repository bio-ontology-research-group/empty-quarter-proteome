#!/bin/bash
#SBATCH --job-name=s4_v4stat
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v4stat_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v4stat_%j.err

set -euo pipefail

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=0
export TMPDIR="${EXP3_ROOT}/scratch_tmp/s4_v4stat_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${EXP3_ROOT}/logs" "${EXP3_ROOT}/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
cd "${EXP3_ROOT}"

echo "EXP3 v4 statistical audit job=${SLURM_JOB_ID} host=$(hostname) start=$(date -u +%FT%TZ)"
python -m pytest -q tests/test_stage4_statistical_audit_v4.py \
  --basetemp="${TMPDIR}/pytest"

python code/stage4_statistical_audit_v4.py \
  --exp3-root "${EXP3_ROOT}" \
  --field-map /scratch/project/k10018/eq-ecology-mag-20261105/results/rund_20261106/matched_profiles.tsv \
  --root-disposition "${EXP3_ROOT}/results/stage4_v4_root_disposition_41Sr3.json" \
  --provisional-disposition "${EXP3_ROOT}/results/stage4_v4_provisional_disposition_second_thermo.json" \
  --stage1-summary /scratch/project/k10018/eq-science-20261108/results/stage1_summary.json \
  --stage1b-summary /scratch/project/k10018/eq-science-20261108/results/stage1b_summary.json \
  --n-perms 9999 \
  --seed 20261109

echo "EXP3 v4 statistical audit complete job=${SLURM_JOB_ID} end=$(date -u +%FT%TZ)"
