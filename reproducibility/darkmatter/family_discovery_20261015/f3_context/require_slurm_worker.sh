#!/bin/bash
# Execute compute only on a worker in the current running Slurm allocation.
set -euo pipefail
compute_host=$(hostname -s)
case "$compute_host" in
  *login*|ilogin*) echo "Refusing compute on login host: $compute_host" >&2; exit 64 ;;
esac
if [[ -z ${SLURM_JOB_ID:-} || ! ${SLURM_JOB_ID} =~ ^[0-9]+$ ]]; then
  echo 'Refusing compute without a Slurm job ID.' >&2; exit 64
fi
job_record=$(scontrol show job "$SLURM_JOB_ID" -o)
if [[ " $job_record " != *' JobState=RUNNING '* ]]; then
  echo 'Refusing compute: allocation is not RUNNING.' >&2; exit 64
fi
allocated_nodes=$(sed -n 's/.* NodeList=\([^ ]*\).*/\1/p' <<< "$job_record")
if [[ -z $allocated_nodes ]] || ! scontrol show hostnames "$allocated_nodes" | grep -Fxq "$compute_host"; then
  echo "Refusing compute: $compute_host is outside allocation $SLURM_JOB_ID." >&2; exit 64
fi
if [[ $# -eq 0 ]]; then
  echo 'Usage: require_slurm_worker.sh COMMAND [ARG ...]' >&2; exit 64
fi
echo "Compute allocation verified: job=$SLURM_JOB_ID host=$compute_host" >&2
exec "$@"
