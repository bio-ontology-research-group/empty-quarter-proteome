#!/bin/bash
#SBATCH --job-name=h4-gpu-comparison
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --gpus=nvidia_h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-1%2
#SBATCH --requeue
#SBATCH --output=/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014/logs/search_%A_%a_%r.out
#SBATCH --error=/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014/logs/search_%A_%a_%r.err
set -euo pipefail
ROOT=/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014
if [[ ${H4_WORKER_VERIFIED:-0} != 1 ]]; then
 export H4_WORKER_VERIFIED=1
 exec "$ROOT/require_slurm_worker.sh" /bin/bash "$ROOT/search.sbatch"
fi
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
cd "$ROOT/payload"
test -f "$ROOT/TRANSFER_COMPLETE"
sha256sum -c payload.sha256 --quiet
case "$SLURM_ARRAY_TASK_ID" in 0) ARM=wc;; 1) ARM=crop;; *) exit 64;; esac
FOLDSEEK="$ROOT/payload/software/foldseek/bin/foldseek"
test "$("$FOLDSEEK" version)" = 941cd33ff0771cd2e3f144e3293e22a2b87e9fda
FINAL="$ROOT/results/$ARM"
if [[ -f "$FINAL/COMPLETE" ]]; then
 (cd "$FINAL" && sha256sum -c result.sha256 --quiet)
 exit 0
fi
WORK="$ROOT/results/.${ARM}_${SLURM_ARRAY_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$WORK"
test ! -e "$WORK/STARTED"
touch "$WORK/STARTED"
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv > "$WORK/gpu.csv"
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv -l 5 > "$WORK/gpu_utilization.csv" &
GPU_MONITOR=$!
trap 'kill "$GPU_MONITOR" 2>/dev/null || true' EXIT
FMT='query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,pident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,cigar'
set -x
"$FOLDSEEK" search "$ROOT/payload/$ARM/query" "$ROOT/payload/pdb/target" "$WORK/alignment" "$WORK/tmp" --alignment-type 1 --gpu 1 --prefilter-mode 1 -e 10 --max-seqs 1000 --threads "$SLURM_CPUS_PER_TASK" -a
"$FOLDSEEK" convertalis "$ROOT/payload/$ARM/query" "$ROOT/payload/pdb/target" "$WORK/alignment" "$WORK/hits.raw.tsv" --format-output "$FMT" --threads "$SLURM_CPUS_PER_TASK"
set +x
kill "$GPU_MONITOR" 2>/dev/null || true
wait "$GPU_MONITOR" 2>/dev/null || true
printf 'job=%s task=%s restart=%s\n' "$SLURM_JOB_ID" "$SLURM_ARRAY_TASK_ID" "${SLURM_RESTART_COUNT:-0}" > "$WORK/execution.txt"
cp "$ROOT/payload/provenance.json" "$WORK/input_provenance.json"
cp "$ROOT/search.sbatch" "$WORK/executed_search.sbatch"
cp "$ROOT/require_slurm_worker.sh" "$WORK/executed_guard.sh"
(cd "$WORK" && sha256sum hits.raw.tsv execution.txt input_provenance.json executed_search.sbatch executed_guard.sh gpu.csv gpu_utilization.csv > result.sha256 && sha256sum -c result.sha256 --quiet)
touch "$WORK/COMPLETE"
test ! -e "$FINAL"
mv "$WORK" "$FINAL"
