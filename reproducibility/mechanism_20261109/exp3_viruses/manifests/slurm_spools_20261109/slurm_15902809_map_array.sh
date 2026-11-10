#!/bin/bash
#SBATCH --job-name=exp3_st4_map
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --mem=350G
#SBATCH --array=0-13
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_map_%A_%a.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_map_%A_%a.err

set -euo pipefail

SAMPLES=(
    "41PRr2"
    "V60PRr2"
    "55Dr2"
    "9Sr2"
    "55Sr3"
    "S49Dr1"
    "48Sr2"
    "F48Sr2"
    "V48Dr1"
    "S48PRr1"
    "41Sr3"
    "10Dr2"
    "43PRr3"
    "44Dr2"
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "========================================================================"
echo "EXP3 Stage 4 Competitive Mapping Job $SLURM_ARRAY_JOB_ID Task $SLURM_ARRAY_TASK_ID: Sample $SAMPLE"
echo "Host: $(hostname) Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity

cd "${EXP3_ROOT}"

python "${EXP3_ROOT}/code/stage4_map_sample_competitive.py" --sample "${SAMPLE}" --threads 32

echo "Completed mapping for ${SAMPLE} at: $(date)"
