#!/bin/bash
#SBATCH --job-name=exp3_followups_v2
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_followups_v2_%j.log
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/logs/stage4_followups_v2_%j.err

set -euo pipefail

echo "========================================================================"
echo "EXP3 Stage 3/4 v2 Followups Job $SLURM_JOB_ID on $(hostname)"
echo "Started at: $(date)"
echo "========================================================================"

export EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
export PYTHONNOUSERSITE=1

source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity

cd "${EXP3_ROOT}"

echo "[Step 1] Running audit_thermo_pair_spacers.py..."
python "${EXP3_ROOT}/code/audit_thermo_pair_spacers.py"

echo "[Step 2] Running stage3_host_link_ledger_v2.py..."
python "${EXP3_ROOT}/code/stage3_host_link_ledger_v2.py"

echo "[Step 3] Running audit_tier1_pairs_evidence.py..."
python "${EXP3_ROOT}/code/audit_tier1_pairs_evidence.py"

echo "[Step 4] Running stage4_downstream_mechanism_v2.py..."
python "${EXP3_ROOT}/code/stage4_downstream_mechanism_v2.py"

echo "[Step 5] Running stage4_targeted_joint_competitive_mapping.py..."
python "${EXP3_ROOT}/code/stage4_targeted_joint_competitive_mapping.py" --threads 32

echo "[Step 6] Running stage4_structural_hallmark_screen.py..."
python "${EXP3_ROOT}/code/stage4_structural_hallmark_screen.py" --threads 16

echo "========================================================================"
echo "All Stage 3/4 v2 Followups Completed at: $(date)"
echo "========================================================================"
