#!/bin/bash
#SBATCH --job-name=exp3_st4_posthoc
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_posthoc_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_posthoc_%j.err

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 4 Posthoc inStrain & Host-Array Provenance Followup"
echo "Job $SLURM_JOB_ID on $(hostname) Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/stage4_posthoc_instrain_and_mapping_followup.py"

echo "Completed at: $(date)"
