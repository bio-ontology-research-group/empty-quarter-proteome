#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_loc_tmbed
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/logs/%x-%A_%a.out
# Stage 2b: TMbed 1.0.2 (CPU) for signal-peptide-positive EQ member sequences without a reusable prediction.
# 32 parts, 8 per array task, 8 concurrent processes x 24 threads. One packed output per array task.
set -euo pipefail
OUT=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation
WORK=$OUT/work
E=/project/k10018/sw/envs/tmbed
T=/dev/shm/eqloc_tmbed_${SLURM_JOB_ID}; rm -rf "$T"; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
cp -r /project/k10018/sw/db/eq-localisation-20261111/hf_home $T/hf_home
export HF_HOME=$T/hf_home HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=24 MKL_NUM_THREADS=24 TOKENIZERS_PARALLELISM=false
A=$SLURM_ARRAY_TASK_ID
pids=()
for j in 0 1 2 3 4 5 6 7; do
  p=$(printf %03d $((A*8+j)))
  if [ -s $WORK/tmbed_parts/part_$p.faa ]; then
    ( $E/bin/tmbed predict -f $WORK/tmbed_parts/part_$p.faa -p $T/part_$p.pred --out-format 0 --no-use-gpu \
        --cpu-fallback -t 24 --batch-size 2000 > $T/part_$p.log 2>&1 || { echo "FAILED part $p"; tail -5 $T/part_$p.log; exit 1; } ) &
    pids+=($!)
  else
    : > $T/part_$p.pred
  fi
done
rc=0
for p in "${pids[@]}"; do wait $p || rc=1; done
[ $rc -eq 0 ]
nin=0; nout=0
for j in 0 1 2 3 4 5 6 7; do
  p=$(printf %03d $((A*8+j)))
  a=$(grep -c '>' $WORK/tmbed_parts/part_$p.faa || true); b=$(grep -c '>' $T/part_$p.pred || true)
  nin=$((nin+a)); nout=$((nout+b))
done
echo "task $A sequences in=$nin predicted=$nout"
[ $nin -eq $nout ]
cat $T/part_*.pred | gzip -c > $WORK/tmbed_new_task$A.pred.gz.tmp && mv $WORK/tmbed_new_task$A.pred.gz.tmp $WORK/tmbed_new_task$A.pred.gz
$E/bin/python -c "import tmbed, torch, transformers; print('tmbed', tmbed.__version__ if hasattr(tmbed,'__version__') else '1.0.2', 'torch', torch.__version__, 'transformers', transformers.__version__)"
echo TMBED_TASK_DONE
