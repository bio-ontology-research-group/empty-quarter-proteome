#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_loc_prep
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/logs/%x-%j.out
# Stage 1: protein selection, SignalP join, TMbed reuse by exact sequence, FASTA files, TMbed to-do parts.
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis3_localisation
PY=/project/k10018/sw/envs/mapping/bin/python
$PY $OUT/01_prep.py --root $ROOT --out $OUT/work --parts 32 --qc-n 1500 \
    --tmbed-dir /project/k10018/work/eq-annot-20261030/results/sp_tm_analysis/tmbed
echo PREP_JOB_DONE
