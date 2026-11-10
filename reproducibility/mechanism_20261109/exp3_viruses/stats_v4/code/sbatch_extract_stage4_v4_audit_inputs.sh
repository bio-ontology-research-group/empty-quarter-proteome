#!/bin/bash
#SBATCH --job-name=s4_v4_extract
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v4_extract_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/s4_v4_extract_%j.err

set -euo pipefail

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=0

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
cd "${EXP3_ROOT}"

echo "EXP3 v4 audit extraction job=${SLURM_JOB_ID} host=$(hostname) start=$(date -u +%FT%TZ)"
python provenance/extract_stage4_v4_audit_inputs.py
echo "EXP3 v4 audit extraction complete job=${SLURM_JOB_ID} end=$(date -u +%FT%TZ)"
