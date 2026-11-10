#!/bin/bash
# Submit Stage 4 competitive mapping pipeline with exact Slurm job dependencies.
# Preserves a live reference build job by default (e.g. 15893933) and snapshots
# repaired viral mapping/validation/aggregation/sequence-followup code.
set -euo pipefail

EXP3_ROOT="${EXP3_ROOT:-/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses}"
cd "${EXP3_ROOT}"

JOB_REF="${1:-15893933}"

echo "========================================================================"
echo "Submitting EXP3 Stage 4 Repaired Competitive Mapping & Follow-up Chain"
echo "Timestamp: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Live/Target Reference Build Job: ${JOB_REF}"
echo "========================================================================"

if [[ "${JOB_REF}" == "--new-builder" ]]; then
    OUT_REF=$(sbatch "${EXP3_ROOT}/code/sbatch_stage4_build_ref.sh")
    JOB_REF=$(echo "${OUT_REF}" | awk '{print $NF}')
    echo "Submitted new Reference Build & Indexing Job: ${JOB_REF}"
else
    scontrol show job "${JOB_REF}" >/dev/null
    echo "Verified existing Reference Build Job in Slurm: ${JOB_REF}"
fi

# 1. Snapshot repaired viral code and record cryptographic hashes
SNAP_DIR="${EXP3_ROOT}/provenance/version_snapshots/pi_viral_repair_20261109"
mkdir -p "${SNAP_DIR}/code" "${SNAP_DIR}/inputs" "${SNAP_DIR}/tests"
OWNED_FILES=(
    "code/stage4_competitive_common.py"
    "code/stage4_map_sample_competitive.py"
    "code/stage4_validate_competitive_reference.py"
    "code/stage4_aggregate_competitive_mapping.py"
    "code/stage4_structural_hallmark_screen.py"
    "code/audit_target_sequence_collisions_v2.py"
    "code/stage4_targeted_joint_competitive_mapping.py"
    "code/sbatch_stage4_validate_ref.sh"
    "code/sbatch_stage4_map_array.sh"
    "code/sbatch_stage4_aggregate.sh"
    "code/sbatch_stage4_sequence_followup.sh"
    "code/submit_stage4_competitive_mapping_pipeline.sh"
    "inputs/stage4_boundary_sensitive_intervals.json"
    "tests/test_viral_mapping_repair.py"
)
for rel in "${OWNED_FILES[@]}"; do
    cp -pf "${EXP3_ROOT}/${rel}" "${SNAP_DIR}/${rel}"
done
sha256sum "${OWNED_FILES[@]}" > "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv"
cp -pf "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv" "${SNAP_DIR}/pi_viral_repair_code_hashes.tsv"

# 2. Submit Post-Build Reference Validator (dependent on Reference Build)
OUT_VAL=$(sbatch --dependency=afterok:${JOB_REF} "${EXP3_ROOT}/code/sbatch_stage4_validate_ref.sh")
JOB_VAL=$(echo "${OUT_VAL}" | awk '{print $NF}')
echo "Job 1 (Post-Build Reference Validation): ${JOB_VAL} (afterok:${JOB_REF})"

# 3. Submit Competitive Mapping Array across 14 samples (dependent on Reference Validator)
OUT_MAP=$(sbatch --dependency=afterok:${JOB_VAL} "${EXP3_ROOT}/code/sbatch_stage4_map_array.sh")
JOB_MAP=$(echo "${OUT_MAP}" | awk '{print $NF}')
echo "Job 2 (Mapping Array across 14 samples, 0-13): ${JOB_MAP} (afterok:${JOB_VAL})"

# 4. Submit Competitive Mapping Aggregation (dependent on Mapping Array)
OUT_AGGR=$(sbatch --dependency=afterok:${JOB_MAP} "${EXP3_ROOT}/code/sbatch_stage4_aggregate.sh")
JOB_AGGR=$(echo "${OUT_AGGR}" | awk '{print $NF}')
echo "Job 3 (Competitive Mapping Aggregation): ${JOB_AGGR} (afterok:${JOB_MAP})"

# 5. Submit Sequence Domain & Existing-Model Follow-up (dependent on Aggregation)
OUT_SEQ=$(sbatch --dependency=afterok:${JOB_AGGR} "${EXP3_ROOT}/code/sbatch_stage4_sequence_followup.sh")
JOB_SEQ=$(echo "${OUT_SEQ}" | awk '{print $NF}')
echo "Job 4 (Sequence Domain & Structural Follow-up): ${JOB_SEQ} (afterok:${JOB_AGGR})"

RECEIPT_FILE="${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission.json"
if [[ -f "${RECEIPT_FILE}" && ! -f "${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission_v1_superseded.json" ]]; then
    cp -pf "${RECEIPT_FILE}" "${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission_v1_superseded.json"
fi

cat > "${RECEIPT_FILE}" <<EOF
{
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "job_ref_build": "${JOB_REF}",
  "job_validate_ref": "${JOB_VAL}",
  "job_map_array": "${JOB_MAP}",
  "job_aggregate": "${JOB_AGGR}",
  "job_sequence_followup": "${JOB_SEQ}",
  "dependencies": {
    "validate_ref": "afterok:${JOB_REF}",
    "map_array": "afterok:${JOB_VAL}",
    "aggregate": "afterok:${JOB_MAP}",
    "sequence_followup": "afterok:${JOB_AGGR}"
  },
  "resource_requests": {
    "validate_ref": {"cpus": 16, "mem": "350G", "time": "04:00:00", "partition": "workq"},
    "map_array": {"cpus": 32, "mem": "350G", "time": "04:00:00", "array": "0-13", "partition": "workq"},
    "aggregate": {"cpus": 4, "mem": "32G", "time": "02:00:00", "partition": "workq"},
    "sequence_followup": {"cpus": 16, "mem": "64G", "time": "04:00:00", "partition": "workq"}
  },
  "samples": [
    "41PRr2", "V60PRr2", "55Dr2", "9Sr2", "55Sr3",
    "S49Dr1", "48Sr2", "F48Sr2", "V48Dr1", "S48PRr1",
    "41Sr3", "10Dr2", "43PRr3", "44Dr2"
  ],
  "snapshot_dir": "${SNAP_DIR}",
  "code_hashes_tsv": "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv",
  "status": "repaired_pipeline_submitted"
}
EOF

echo "Pipeline submission receipt written to: ${RECEIPT_FILE}"
echo "========================================================================"
