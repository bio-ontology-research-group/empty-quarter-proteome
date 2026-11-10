#!/bin/bash
set -euo pipefail

EXP3_ROOT="${EXP3_ROOT:-/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses}"
cd "$EXP3_ROOT"
mkdir -p "$EXP3_ROOT/provenance" "$EXP3_ROOT/logs" "$EXP3_ROOT/results"

S1_COLLECT_JOB_ID="${1:-}"
if [[ -z "$S1_COLLECT_JOB_ID" ]]; then
  echo "Usage: $0 <replacement_stage1_collect_job_id>" >&2
  exit 1
fi

J_S2=$(sbatch --parsable --dependency=afterok:${S1_COLLECT_JOB_ID} "$EXP3_ROOT/code/sbatch_stage2_search_and_verify.sh")
J_S3=$(sbatch --parsable --dependency=afterok:${J_S2} "$EXP3_ROOT/code/sbatch_stage3_ledger.sh")
J_S4=$(sbatch --parsable --dependency=afterok:${J_S3} "$EXP3_ROOT/code/sbatch_stage4_downstream.sh")

cat <<EOF > "$EXP3_ROOT/provenance/stage2_to_4_submission_receipts.txt"
timestamp=$(date -Is)
stage1_collect_dependency_job_id=${S1_COLLECT_JOB_ID}
stage2_search_and_verify_job_id=${J_S2} (afterok:${S1_COLLECT_JOB_ID})
stage3_ledger_job_id=${J_S3} (afterok:${J_S2})
stage4_downstream_job_id=${J_S4} (afterok:${J_S3})
EOF

cat "$EXP3_ROOT/provenance/stage2_to_4_submission_receipts.txt"
