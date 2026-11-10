#!/bin/bash
# Submit Stage 4 competitive mapping pipeline with exact Slurm job dependencies.
# Preserves completed Reference Build (15893933), completed Sequence-Domain Follow-up (15894662),
# and historical chains (15894659 -> 15894660 -> 15894661 and 15902808 -> 15902809 -> 15902810),
# cancels only stranded 15902809/15902810, snapshots the boundary-length-repaired code,
# submits new Validator -> 14-sample Mapping Array -> Aggregator, and captures Slurm batch spools.
set -euo pipefail

EXP3_ROOT="${EXP3_ROOT:-/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses}"
cd "${EXP3_ROOT}"

JOB_REF="${1:-15893933}"
JOB_SEQ_COMPLETED="${2:-15894662}"

echo "========================================================================"
echo "Submitting EXP3 Stage 4 Boundary-Length-Repaired Validator & Mapping Chain"
echo "Timestamp: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Completed Reference Build Job: ${JOB_REF}"
echo "Completed Sequence-Domain Job: ${JOB_SEQ_COMPLETED}"
echo "========================================================================"

# 1. Preserve historical snapshots and receipts before creating boundary-length-repair snapshot
mkdir -p "${EXP3_ROOT}/provenance/slurm_spools_20261109"
if [[ -f "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv" && ! -f "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes_v1_15894659.tsv" ]]; then
    cp -pf "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv" "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes_v1_15894659.tsv"
fi
RECEIPT_FILE="${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission.json"
if [[ -f "${RECEIPT_FILE}" && ! -f "${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission_v3_15902808_superseded.json" ]]; then
    cp -pf "${RECEIPT_FILE}" "${EXP3_ROOT}/provenance/stage4_competitive_pipeline_submission_v3_15902808_superseded.json"
fi

# 2. Cancel ONLY stranded failed-dependent jobs 15902809 and 15902810
scancel 15902809 15902810 || true
echo "Cancelled stranded dependent jobs: 15902809 15902810"

# 3. Snapshot boundary-length-repaired viral code and record cryptographic hashes
SNAP_DIR="${EXP3_ROOT}/provenance/version_snapshots/pi_viral_boundary_length_repair_20261109"
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
sha256sum "${OWNED_FILES[@]}" > "${EXP3_ROOT}/provenance/pi_viral_boundary_length_repair_code_hashes.tsv"
cp -pf "${EXP3_ROOT}/provenance/pi_viral_boundary_length_repair_code_hashes.tsv" "${SNAP_DIR}/pi_viral_boundary_length_repair_code_hashes.tsv"

# 4. Submit Post-Build Reference Validator (immediately if 15893933 already COMPLETED, else afterok)
if scontrol show job "${JOB_REF}" 2>/dev/null | grep -qE "JobState=(RUNNING|PENDING)"; then
    VAL_DEP="afterok:${JOB_REF}"
    OUT_VAL=$(sbatch --dependency="${VAL_DEP}" "${EXP3_ROOT}/code/sbatch_stage4_validate_ref.sh")
else
    test -f "${EXP3_ROOT}/results/stage4_competitive_reference/REFERENCE_BUILD_DONE"
    test -s "${EXP3_ROOT}/results/stage4_competitive_reference/unified_competitive_reference.receipt.json"
    VAL_DEP="satisfied_completed_${JOB_REF}"
    OUT_VAL=$(sbatch "${EXP3_ROOT}/code/sbatch_stage4_validate_ref.sh")
fi
JOB_VAL=$(echo "${OUT_VAL}" | awk '{print $NF}')
echo "Job 1 (Post-Build Reference Validation): ${JOB_VAL} (${VAL_DEP})"

# 5. Submit Competitive Mapping Array across 14 samples (dependent afterok on Reference Validator)
OUT_MAP=$(sbatch --dependency=afterok:${JOB_VAL} "${EXP3_ROOT}/code/sbatch_stage4_map_array.sh")
JOB_MAP=$(echo "${OUT_MAP}" | awk '{print $NF}')
echo "Job 2 (Mapping Array across 14 samples, 0-13): ${JOB_MAP} (afterok:${JOB_VAL})"

# 6. Submit Competitive Mapping Aggregation (dependent afterok on Mapping Array)
OUT_AGGR=$(sbatch --dependency=afterok:${JOB_MAP} "${EXP3_ROOT}/code/sbatch_stage4_aggregate.sh")
JOB_AGGR=$(echo "${OUT_AGGR}" | awk '{print $NF}')
echo "Job 3 (Competitive Mapping Aggregation): ${JOB_AGGR} (afterok:${JOB_MAP})"

# Re-snapshot submit_stage4_competitive_mapping_pipeline.sh and hashes so snapshot is exact
cp -pf "${EXP3_ROOT}/code/submit_stage4_competitive_mapping_pipeline.sh" "${SNAP_DIR}/code/submit_stage4_competitive_mapping_pipeline.sh"
sha256sum "${OWNED_FILES[@]}" > "${EXP3_ROOT}/provenance/pi_viral_boundary_length_repair_code_hashes.tsv"
cp -pf "${EXP3_ROOT}/provenance/pi_viral_boundary_length_repair_code_hashes.tsv" "${SNAP_DIR}/pi_viral_boundary_length_repair_code_hashes.tsv"

# 7. Capture Slurm batch spools and job metadata (non-fatal after job submission to prevent duplicate sbatch)
scontrol write batch_script "${JOB_VAL}" "${EXP3_ROOT}/provenance/slurm_spools_20261109/slurm_${JOB_VAL}_validate_ref.sh" || true
scontrol write batch_script "${JOB_MAP}" "${EXP3_ROOT}/provenance/slurm_spools_20261109/slurm_${JOB_MAP}_map_array.sh" || true
scontrol write batch_script "${JOB_AGGR}" "${EXP3_ROOT}/provenance/slurm_spools_20261109/slurm_${JOB_AGGR}_aggregate.sh" || true
{
    scontrol show job "${JOB_VAL}" || true
    scontrol show job "${JOB_MAP}" || true
    scontrol show job "${JOB_AGGR}" || true
} > "${EXP3_ROOT}/provenance/slurm_spools_20261109/scontrol_show_submitted_${JOB_VAL}_${JOB_MAP}_${JOB_AGGR}.txt" || true

cat > "${RECEIPT_FILE}" <<EOF
{
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "job_ref_build": "${JOB_REF}",
  "job_validate_ref": "${JOB_VAL}",
  "job_map_array": "${JOB_MAP}",
  "job_aggregate": "${JOB_AGGR}",
  "job_sequence_followup": "${JOB_SEQ_COMPLETED}",
  "dependencies": {
    "validate_ref": "${VAL_DEP}",
    "map_array": "afterok:${JOB_VAL}",
    "aggregate": "afterok:${JOB_MAP}",
    "sequence_followup": "completed_independent_of_mapping_15894662"
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
  "code_hashes_tsv": "${EXP3_ROOT}/provenance/pi_viral_boundary_length_repair_code_hashes.tsv",
  "historical_chain_15894659": {
    "job_ref_build": "15893933",
    "job_validate_ref": "15894659",
    "job_map_array_cancelled": "15894660",
    "job_aggregate_cancelled": "15894661",
    "job_sequence_followup_completed": "15894662",
    "historical_snapshot_dir": "${EXP3_ROOT}/provenance/version_snapshots/pi_viral_repair_20261109",
    "historical_code_hashes_tsv": "${EXP3_ROOT}/provenance/pi_viral_repair_code_hashes.tsv",
    "failure_reason": "Validator 15894659 failed at line 229 because touch-created 0-byte REFERENCE_BUILD_DONE sentinel was checked with st_size == 0."
  },
  "historical_chain_15902808": {
    "job_ref_build": "15893933",
    "job_validate_ref": "15902808",
    "job_map_array_cancelled": "15902809",
    "job_aggregate_cancelled": "15902810",
    "job_sequence_followup_completed": "15894662",
    "historical_snapshot_dir": "${EXP3_ROOT}/provenance/version_snapshots/pi_viral_sentinel_repair_20261109",
    "historical_code_hashes_tsv": "${EXP3_ROOT}/provenance/pi_viral_sentinel_repair_code_hashes.tsv",
    "failure_reason": "Validator 15902808 passed sentinel, crosswalk, FAI, and 25.5 GB reference scan, then failed at line 285 with KeyError: 'target_length' because load_boundary_candidates validated target_length=67924 in inputs/stage4_boundary_sensitive_intervals.json but omitted it from records.append; fixed in pi_viral_boundary_length_repair_20261109."
  },
  "status": "repaired_boundary_length_pipeline_submitted"
}
EOF

echo "Pipeline submission receipt written to: ${RECEIPT_FILE}"
echo "========================================================================"
