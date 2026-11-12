#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_fam_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis2_families/logs/%x-%j.out
# Stage 3: family classes, sharing fractions, controls, block bootstrap.
# Usage: sbatch --export=ALL,CLUSTER_JOB=<id>,SEARCH_JOB=<id> sbatch_03_analysis.sh
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONWARNINGS=ignore
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis2_families
RES=$OUT/results
PY=/project/k10018/sw/envs/mapping/bin/python
[ -s $OUT/clusters/families_id30_c80.tsv ] && [ -s $OUT/clusters/families_id50_c80.tsv ]
[ -s $OUT/search/eq_vs_ext_by_group.tsv ]
EXTRA=$(printf '{"slurm_job_id_cluster":"%s","slurm_job_id_search":"%s","mmseqs_version":"18.8cc5c","mmseqs_image":"/project/k10018/sw/images/eq-annot/mmseqs2.sif","mmseqs_image_sha256":"5b79474ef8045ec82715c357691386d9301b98448285293f94348e8b69dc1165","cluster_params":"mmseqs cluster --min-seq-id {0.30|0.50} -c 0.8 --cov-mode 0 -s 7.5 -e 1e-3 --cluster-mode 0 --cluster-reassign 1","search_params":"mmseqs search (EQ proteins vs external proteins) -s 7.5 -e 1e-5 --cov-mode 2 -c 0.5 --max-seqs 4000"}' "${CLUSTER_JOB:-NA}" "${SEARCH_JOB:-NA}")
$PY $OUT/family_sharing.py --root $ROOT --out $RES.tmp.$SLURM_JOB_ID \
    --families id30_c80=$OUT/clusters/families_id30_c80.tsv id50_c80=$OUT/clusters/families_id50_c80.tsv \
    --search $OUT/search/eq_vs_ext_by_group.tsv --max-seqs 4000 \
    --draws 2000 --split-reps 500 --seed 20261110 --receipt-extra "$EXTRA"
[ -f $RES.tmp.$SLURM_JOB_ID/receipt.json ]
[ ! -e $RES ] || mv $RES $RES.superseded.$SLURM_JOB_ID
mv $RES.tmp.$SLURM_JOB_ID $RES
echo FAMILY_ANALYSIS_DONE
