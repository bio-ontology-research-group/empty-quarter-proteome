#!/bin/bash
# Submit prep -> shard array -> merge with dependencies (run on the IBEX login node; it only calls sbatch).
#   test: RECON_TEST=1 RECON_ROOT=$D/test bash genecall_reconcile_submit.sh
#   full: bash genecall_reconcile_submit.sh
set -euo pipefail
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/genecall_reconcile_20261102
S=$D/scripts/genecall_reconcile.sbatch
ROOT=${RECON_ROOT:-$D}; mkdir -p "$ROOT/logs"
export RECON_ROOT=$ROOT RECON_TEST=${RECON_TEST:-0}
if [ "$RECON_TEST" = 1 ]; then NARR=0; T1=00:30:00; T2=00:30:00; T3=02:00:00
else NARR=$(( ${RECON_NSHARD:-48} - 1 )); T1=03:00:00; T2=06:00:00; T3=12:00:00; fi
j1=$(sbatch --parsable --export=ALL -c 16 --mem=64G -t $T1 -o "$ROOT/logs/prep_%j.out" $S prep)
j2=$(sbatch --parsable --export=ALL -d afterok:$j1 --array=0-$NARR%24 -c 4 --mem=32G -t $T2 \
     -o "$ROOT/logs/shard_%A_%a.out" $S shard)
j3=$(sbatch --parsable --export=ALL -d afterok:$j2 -c 8 --mem=192G -t $T3 -o "$ROOT/logs/merge_%j.out" $S merge)
echo "prep $j1  shard $j2  merge $j3  root $ROOT" | tee -a "$ROOT/logs/submitted.txt"
