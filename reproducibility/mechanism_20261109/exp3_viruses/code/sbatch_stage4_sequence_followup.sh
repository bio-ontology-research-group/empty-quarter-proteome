#!/bin/bash
#SBATCH --job-name=exp3_seq_domains
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_sequence_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_sequence_%j.err

set -euo pipefail
export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1
export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity
cd "${EXP3_ROOT}"
python "${EXP3_ROOT}/code/stage4_structural_hallmark_screen.py" --threads 16
