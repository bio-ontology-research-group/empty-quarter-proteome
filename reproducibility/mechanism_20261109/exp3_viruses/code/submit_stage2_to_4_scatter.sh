#!/bin/bash
# submit_stage2_to_4_scatter.sh
# Authoritative launcher for Experiment 3 scalable scattered Stage 2 to Stage 4 pipeline
set -euo pipefail

EXP3_ROOT="/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"
cd "$EXP3_ROOT"

echo "[$(date -Is)] Verifying upstream prerequisites..."
if [[ ! -f "$EXP3_ROOT/results/STAGE0_DONE" ]]; then
    echo "ERROR: Missing $EXP3_ROOT/results/STAGE0_DONE" >&2
    exit 1
fi
if [[ ! -f "$EXP3_ROOT/results/STAGE1_COLLECT_DONE" ]]; then
    echo "ERROR: Missing $EXP3_ROOT/results/STAGE1_COLLECT_DONE" >&2
    exit 1
fi
echo "[$(date -Is)] Upstream Stage 0 and Stage 1 verified complete."

# 1. Submit Stage 2 catalog search (6 pools, 1 node, 45m)
echo "[$(date -Is)] Submitting Stage 2 catalog pool search..."
J_CAT=$(sbatch --parsable "$EXP3_ROOT/code/sbatch_stage2_catalog.sh")
echo "Submitted Stage 2 Catalog Job: ${J_CAT}"

# 2. Submit Stage 2 assembly array (8 shards, %4 concurrency, indexed makeblastdb + blastn -db)
echo "[$(date -Is)] Submitting Stage 2 raw assembly search array (8 shards)..."
J_ARR=$(sbatch --parsable "$EXP3_ROOT/code/sbatch_stage2_asm_array.sh")
echo "Submitted Stage 2 Assembly Array Job: ${J_ARR}"

# 3. Submit Stage 2 collect & annotation (dependent on both catalog and all array shards)
echo "[$(date -Is)] Submitting Stage 2 collect & annotation..."
J_COL=$(sbatch --parsable --dependency=afterok:${J_CAT}:${J_ARR} "$EXP3_ROOT/code/sbatch_stage2_collect.sh")
echo "Submitted Stage 2 Collect Job: ${J_COL} (dependency=afterok:${J_CAT}:${J_ARR})"

# 4. Submit Stage 3 host-link ledger (dependent on Stage 2 collect)
echo "[$(date -Is)] Submitting Stage 3 host-link ledger..."
J_S3=$(sbatch --parsable --dependency=afterok:${J_COL} "$EXP3_ROOT/code/sbatch_stage3_ledger.sh")
echo "Submitted Stage 3 Ledger Job: ${J_S3} (dependency=afterok:${J_COL})"

# 5. Submit Stage 4 downstream mechanism analyses (dependent on Stage 3 ledger)
echo "[$(date -Is)] Submitting Stage 4 downstream mechanism analyses..."
J_S4=$(sbatch --parsable --dependency=afterok:${J_S3} "$EXP3_ROOT/code/sbatch_stage4_downstream.sh")
echo "Submitted Stage 4 Mechanism Job: ${J_S4} (dependency=afterok:${J_S3})"

# Record submission receipt
RECEIPT_FILE="$EXP3_ROOT/provenance/stage2_to_4_submission_receipts.txt"
mkdir -p "$EXP3_ROOT/provenance"
cat <<EOF >> "$RECEIPT_FILE"
--------------------------------------------------------------------------------
Timestamp UTC: $(date -u -Is)
Timestamp Local (Riyadh): $(date -Is)
Launcher: submit_stage2_to_4_scatter.sh
Architecture: Scattered Slurm Execution with Indexed DB Parallelism
Stage 2 Catalog Search Job: ${J_CAT} (1 node, 192 CPUs, 45m)
Stage 2 Assembly Array Job:  ${J_ARR} (8 shards [0-7%4], 192 CPUs/node, indexed makeblastdb + 8-thread blastn)
Stage 2 Collector Job:      ${J_COL} (dependency=afterok:${J_CAT}:${J_ARR})
Stage 3 Host-link Ledger:   ${J_S3} (dependency=afterok:${J_COL})
Stage 4 Downstream Analysis:${J_S4} (dependency=afterok:${J_S3})
Raw Output Root:            results/stage2_raw_outputs_v2
Manifest Author:            Single Authoritative Collector (no shard collisions)
Exact Hashing:              Streamed SHA-256 for all queries and subjects without size cutoff
Container Verified Digest:  6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f
Tool Probed Version:        blastn: 2.17.0+
--------------------------------------------------------------------------------
EOF

echo "[$(date -Is)] Successfully submitted scattered pipeline."
echo "Receipt appended to: ${RECEIPT_FILE}"
