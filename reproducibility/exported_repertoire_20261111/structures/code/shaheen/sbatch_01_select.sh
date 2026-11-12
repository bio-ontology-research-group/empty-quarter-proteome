#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_a4_select
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# Step 1: representatives, mature sequences, conserved neighbours, complex pairs and controls.
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
A4=$ROOT/analysis4_structures
PY=/project/k10018/sw/envs/mapping/bin/python
date -Is; hostname
$PY $A4/scripts/01_select.py --root $ROOT --out $A4/work/step1 --seed 20261111 --min-controls 120
date -Is
