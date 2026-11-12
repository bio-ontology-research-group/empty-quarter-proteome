#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_loc_merops
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/logs/%x-%j.out
# Stage 2a (rerun of the MEROPS part only; it failed silently inside sbatch_02_search.sh, job 15955356):
# DIAMOND blastp of EQ member proteins against MEROPS pepunit.
set -euo pipefail
WORK=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/work
IMG=/project/k10018/sw/images/eq-annot
LDB=/project/k10018/sw/db/eq-localisation-20261111
module load singularity >/dev/null 2>&1
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C
sif(){ local n=$1; shift; singularity exec -B /project,/scratch,/dev/shm "$IMG/$n.sif" "$@"; }
T=/dev/shm/eqloc_${SLURM_JOB_ID}; rm -rf "$T"; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
# pepunit.lib wraps sequences and contains a few non-residue characters; normalise to plain FASTA first
awk '/^>/{print; next} {gsub(/[^A-Za-z]/,""); if(length($0)) print toupper($0)}' $LDB/pepunit.lib > $T/pepunit.faa
grep -c '>' $T/pepunit.faa
sif diamond diamond makedb --in $T/pepunit.faa -d $T/pepunit -p 64 2>&1 | tail -5
sif diamond diamond blastp -q $WORK/eq_members_full.faa -d $T/pepunit --very-sensitive -e 1e-3 -k 5 -p 128 \
     --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \
     --tmpdir $T -o $T/merops.m8 2>&1 | tail -8
wc -l $T/merops.m8
gzip -c $T/merops.m8 > $WORK/merops_pepunit.m8.gz.tmp && mv $WORK/merops_pepunit.m8.gz.tmp $WORK/merops_pepunit.m8.gz
echo MEROPS_DONE
