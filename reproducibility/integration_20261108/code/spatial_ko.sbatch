#!/bin/bash
#SBATCH -A k10018
#SBATCH -p workq
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -c 1
#SBATCH --mem=370000M
#SBATCH -t 12:00:00
#SBATCH -J eq-prot-spatial-ko
#SBATCH -o /scratch/project/k10018/eq-proteome-integration-20261106/logs/eq-prot-spatial-ko-%j.log
set -euo pipefail
cd /scratch/project/k10018/eq-proteome-integration-20261106
export TMPDIR="$PWD/tmp" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/scratch/project/k10018/eq-ecology-mag-20261105/code/deps
PY=/project/k10018/sw/envs/gspa-v3/bin/python
"$PY" code/spatial_sensitivity.py ko
