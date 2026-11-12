#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_fam_search
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis2_families/logs/%x-%j.out
# Stage 2: sensitive search of every desert (EQ) protein against every external protein.
# This is a superset of "dark exported proteins without external family members plus a matched
# sample of the other classes": all classes are searched under identical settings.
set -euo pipefail
module load singularity >/dev/null 2>&1 || true
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C LC_ALL=C
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis2_families
SIF=/project/k10018/sw/images/eq-annot/mmseqs2.sif
PY=/project/k10018/sw/envs/mapping/bin/python
THREADS=${SLURM_CPUS_PER_TASK:-192}
MM="singularity exec -B /project,/scratch $SIF mmseqs"
S=$OUT/search
TMP=$OUT/tmp_search_$SLURM_JOB_ID
mkdir -p "$S" "$TMP"
trap 'rm -rf "$TMP"' EXIT
echo "mmseqs version: $($MM version)"; date -Is
$PY $OUT/prep_search_inputs.py --root $ROOT --out $S
set -x
$MM createdb $S/eq_proteins.faa $S/eqDB -v 2
$MM createdb $S/external_proteins.faa $S/extDB -v 2
rm -f $S/eq_vs_ext*
$MM search $S/eqDB $S/extDB $S/eq_vs_ext $TMP/s -s 7.5 -e 1e-5 --cov-mode 2 -c 0.5 \
    --max-seqs 4000 --threads $THREADS -v 2
date -Is
$MM convertalis $S/eqDB $S/extDB $S/eq_vs_ext $TMP/hits.m8 \
    --format-output query,target,pident,alnlen,qcov,tcov,evalue,bits --threads $THREADS -v 2
set +x
date -Is
wc -l $TMP/hits.m8
# convertalis output is grouped by query; verify before the streaming reduction.
cut -f1 $TMP/hits.m8 | uniq | sort | uniq -d | head -3 > $TMP/dupq.txt
[ ! -s $TMP/dupq.txt ] || { echo "queries not contiguous"; exit 4; }
awk -f $OUT/reduce_hits.awk $S/genome_groups.tsv $TMP/hits.m8 > $S/eq_vs_ext_by_group.tsv.tmp
mv $S/eq_vs_ext_by_group.tsv.tmp $S/eq_vs_ext_by_group.tsv
wc -l $S/eq_vs_ext_by_group.tsv
rm -f $S/eq_proteins.faa $S/external_proteins.faa
date -Is
echo SEARCH_DONE
