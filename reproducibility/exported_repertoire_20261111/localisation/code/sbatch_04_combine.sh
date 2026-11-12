#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_loc_combine
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/logs/%x-%j.out
# Stage 3: per-protein features (02_features.py), family calls, summaries, checks, receipt (03_combine.py).
# Job IDs of the earlier stages are read from jobs.json.
set -euo pipefail
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 PYTHONWARNINGS=ignore
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis3_localisation
WORK=$OUT/work
PY=/project/k10018/sw/envs/mapping/bin/python
for f in pfam_vs_members.domtbl.gz dbcan_vs_members.domtbl.gz sorting_models_vs_eq_all.domtbl.gz merops_pepunit.m8.gz \
         tmbed_new_task0.pred.gz tmbed_new_task1.pred.gz tmbed_new_task2.pred.gz tmbed_new_task3.pred.gz; do [ -e $WORK/$f ]; done
cat $WORK/tmbed_new_task[0-3].pred.gz > $WORK/tmbed_new.pred.gz
$PY $OUT/02_features.py --work $WORK --threads 96
EXTRA=$(cat $OUT/jobs.json)
$PY $OUT/03_combine.py --root $ROOT --out $OUT --clans /project/k10018/sw/db/eq-localisation-20261111/Pfam-A.clans.tsv.gz \
    --draws 2000 --seed 20261111 --receipt-extra "$EXTRA"
rm -rf $WORK/tmbed_parts
echo COMBINE_JOB_DONE
