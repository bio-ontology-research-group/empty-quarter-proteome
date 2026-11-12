#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_a4_realign
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# Rebuild alignments from the existing search results. TAG=own (own proteins) or TAG=ownenv (plus the
# local colabfold_search rows from UniRef30/ColabFold envdb as unpaired sequences).
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
A4=$ROOT/analysis4_structures
PY=/project/k10018/sw/envs/mapping/bin/python
W=$A4/work/search
TAG=${TAG:-own}
EXTRA=""
[ "$TAG" = ownenv ] && EXTRA="--extra-dir $A4/work/cfsearch"
date -Is
$PY $A4/scripts/03_build_alignments.py --step1 $A4/work/step1 --msa-db $W/msa --qdb-lookup $W/qDB.lookup \
    --out $A4/results/alignments --tag $TAG --threads 96 $EXTRA
date -Is; echo REALIGN_DONE
