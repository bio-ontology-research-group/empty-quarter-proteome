#!/bin/bash
#SBATCH --job-name=exp3_raw_cov_qc
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_raw_vs_mapq20_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_raw_vs_mapq20_%j.err

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 4 Raw-vs-MAPQ20 Focal Viral Coverage QC (14 samples x 7 targets)"
echo "Job $SLURM_JOB_ID on $(hostname) Started at: $(date -Is)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/tmp/stage4_raw_vs_mapq20_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${EXP3_ROOT}/results/stage4_raw_vs_mapq20" "${EXP3_ROOT}/logs"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/stage4_raw_vs_mapq20_focal_coverage.py" --threads 14

echo "Completed at: $(date -Is)"
