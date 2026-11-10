#!/bin/bash
#SBATCH --job-name=exp3_host_prov
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_host_prov_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_host_prov_%j.err

set -euo pipefail

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${EXP3_ROOT}/logs" "${EXP3_ROOT}/results"

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

python3 "${EXP3_ROOT}/code/extract_stage4_host_array_contig_provenance.py"
rm -rf "${TMPDIR}"
