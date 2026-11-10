#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_ext_export_cmp
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/logs/%x-%j.out
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis1_export/harmonised_20261110
/project/k10018/sw/envs/mapping/bin/python "$ROOT/analysis1_export/export_comparison.py" \
    --root "$ROOT" --out "$OUT.tmp.$SLURM_JOB_ID"
[ -f "$OUT.tmp.$SLURM_JOB_ID/receipt.json" ]
[ ! -e "$OUT" ] || mv "$OUT" "$OUT.superseded.$SLURM_JOB_ID"
mv "$OUT.tmp.$SLURM_JOB_ID" "$OUT"
echo EXPORT_COMPARISON_DONE
