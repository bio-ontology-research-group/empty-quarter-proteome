#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_fam_cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis2_families/logs/%x-%j.out
# Stage 1: MMseqs2 protein families over all 3,617,011 proteins of the 1,957 genomes.
# The FASTA holds exactly the genes of the 1,957 retained genomes (checked in stage 2).
set -euo pipefail
module load singularity >/dev/null 2>&1 || true
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis2_families
SIF=/project/k10018/sw/images/eq-annot/mmseqs2.sif
FAA=$ROOT/genecall/harmonized_proteins.faa
THREADS=${SLURM_CPUS_PER_TASK:-192}
MM="singularity exec -B /project,/scratch $SIF mmseqs"
DB=$OUT/mmseqs_db
TMP=$OUT/tmp_cluster_$SLURM_JOB_ID
mkdir -p "$DB" "$TMP" "$OUT/clusters"
trap 'rm -rf "$TMP"' EXIT
echo "mmseqs version: $($MM version)"; date -Is
set -x
[ -s $DB/all.dbtype ] || $MM createdb "$FAA" $DB/all -v 2
for ID in 30 50; do
  if [ ! -s $OUT/clusters/families_id${ID}_c80.tsv ]; then
    rm -f $DB/clu${ID}*
    $MM cluster $DB/all $DB/clu${ID} $TMP/c${ID} \
        --min-seq-id 0.${ID} -c 0.8 --cov-mode 0 -s 7.5 -e 1e-3 \
        --cluster-mode 0 --cluster-reassign 1 --threads $THREADS -v 2
    $MM createtsv $DB/all $DB/all $DB/clu${ID} $OUT/clusters/families_id${ID}_c80.tsv.tmp --threads $THREADS -v 2
    mv $OUT/clusters/families_id${ID}_c80.tsv.tmp $OUT/clusters/families_id${ID}_c80.tsv
    rm -rf $TMP/c${ID}
  fi
  date -Is
done
set +x
wc -l $OUT/clusters/families_id*_c80.tsv
cut -f1 $OUT/clusters/families_id30_c80.tsv | uniq | sort -u | wc -l
echo CLUSTER_DONE
