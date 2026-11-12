#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_a4_cfsearch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# Step 3 (optional depth): local colabfold_search of all queries against the local copies of
# UniRef30 2302 and the ColabFold environmental database 202108. No web service is contacted.
set -euo pipefail
module load singularity >/dev/null 2>&1 || true
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C LC_ALL=C
ROOT=/scratch/project/k10018/eq-protein-external-20261110
A4=$ROOT/analysis4_structures
CF=/project/k10018/sw/images/eq-annot/colabfold.sif
DB=/project/k10018/sw/db/colabfold
OUT=$A4/work/cfsearch
rm -rf "$OUT"; mkdir -p "$OUT"
date -Is; hostname; free -g | head -2
singularity exec -B /project,/scratch $CF mmseqs version
singularity exec -B /project,/scratch $CF colabfold_search $A4/work/step1/queries.faa $DB $OUT \
    --db1 uniref30_2302_db --db3 colabfold_envdb_202108_db --use-env 1 --use-templates 0 --filter 1 \
    --db-load-mode 0 --threads ${SLURM_CPUS_PER_TASK:-192}
date -Is
ls $OUT | head; ls $OUT/*.a3m | wc -l
echo CFSEARCH_DONE
