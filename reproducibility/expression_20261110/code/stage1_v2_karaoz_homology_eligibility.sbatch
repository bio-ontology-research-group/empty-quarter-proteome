#!/bin/bash
#SBATCH --job-name=eq_exp4_stage1_v2
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/expression/logs/stage1_v2_karaoz_homology_%j.out
#SBATCH --error=/scratch/project/k10018/eq-protein-external-20261110/expression/logs/stage1_v2_karaoz_homology_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/gspa-tools

echo "Running negative-path unit tests before Stage 1 v2 on $(hostname) at $(date -Is)"
python3 -m unittest discover -v -s /scratch/project/k10018/eq-protein-external-20261110/expression/code -p "test_*.py"

echo "Starting Stage 1 v2 Karaoz 2026 corrected homology & eligibility audit at $(date -Is)"
python3 /scratch/project/k10018/eq-protein-external-20261110/expression/code/stage1_v2_karaoz_homology_eligibility.py
echo "Finished Stage 1 v2 at $(date -Is)"
