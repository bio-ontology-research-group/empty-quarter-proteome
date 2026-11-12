#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_a4_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# Step 2: iterative profile search of every query (mature sequence) against all 3,617,011 proteins,
# a3m per query, then alignments/depths (own proteins only) and Pfam profiles of the queries.
set -euo pipefail
module load singularity >/dev/null 2>&1 || true
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C LC_ALL=C
ROOT=/scratch/project/k10018/eq-protein-external-20261110
A4=$ROOT/analysis4_structures
PY=/project/k10018/sw/envs/mapping/bin/python
SIF=/project/k10018/sw/images/eq-annot/mmseqs2.sif
MM="singularity exec -B /project,/scratch $SIF mmseqs"
T=${SLURM_CPUS_PER_TASK:-192}
ALL=$ROOT/analysis2_families/mmseqs_db/all
W=$A4/work/search
TMP=$A4/work/tmp_search_$SLURM_JOB_ID
mkdir -p "$W" "$TMP"
trap 'rm -rf "$TMP"' EXIT
date -Is; hostname; echo "mmseqs $($MM version)"
if [ ! -s $W/msa.index ]; then
  rm -f $W/qDB* $W/res* $W/msa*
  $MM createdb $A4/work/step1/queries.faa $W/qDB -v 2
  $MM search $W/qDB $ALL $W/res $TMP/s --num-iterations 3 -s 7.5 -e 1e-3 --e-profile 1e-3 \
      --max-seqs 20000 -a 1 --threads $T -v 2
  date -Is
  $MM result2msa $W/qDB $ALL $W/res $W/msa --msa-format-mode 6 --filter-msa 0 --threads $T -v 2
  $MM convertalis $W/qDB $ALL $W/res $W/hits_last_iteration.m8 \
      --format-output query,target,pident,alnlen,qcov,tcov,evalue,bits --threads $T -v 2
  gzip -f $W/hits_last_iteration.m8
fi
date -Is
$PY $A4/scripts/03_build_alignments.py --step1 $A4/work/step1 --msa-db $W/msa --qdb-lookup $W/qDB.lookup \
    --out $A4/results/alignments --tag own --threads 96
date -Is
# Pfam-A profiles (gathering thresholds) on all query sequences; failure here does not fail the job
HM=/project/k10018/sw/images/eq-annot/hmmer.sif
( singularity exec -B /project,/scratch $HM hmmsearch --cut_ga --cpu 64 --noali \
    --domtblout $A4/results/queries_pfam.domtbl -o /dev/null \
    /project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm $A4/work/step1/queries.faa \
  && echo PFAM_DONE ) || echo PFAM_FAILED
date -Is
echo STEP2_DONE
