#!/bin/bash
#SBATCH --job-name=exp1_st4_def_v3
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_v3_%j.out
#SBATCH --error=/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/logs/exp1_st4_def_v3_%j.err

set -euo pipefail
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

BASE="/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic"
export TMPDIR="${BASE}/tmp/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}" "${BASE}/results" "${BASE}/logs"
export OMP_NUM_THREADS=1

echo "[$(date -Is)] Starting exp1_metabolic Stage 4 Candidate Deficiency Search (v3) on $(hostname)"

python3 "${BASE}/code/stage4_candidate_deficiency_search.py" \
  --instances-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_instances.tsv" \
  --summary-tsv "${BASE}/results/exp1_stage4_supported_candidate_deficiencies_summary.tsv" \
  --gapseq-rxn-fasta-dir /project/k10018/sw/db/gapseq_seq/Archaea/rxn \
  --refined-mags-tar /project/k10018/work/eq-mag-refine-20261028/results/c07_mag_refine/refined_mags.tar.gz \
  --assemblies-dir /project/k10018/work/eq-annot-20261030/results/c05_megahit/out \
  --gap-context-tsv-gz /scratch/project/k10018/eq-archaeal-function-20261106/results/gapseq_pathway_hole_synteny_context.tsv.gz \
  --c18-candidates-tsv /scratch/project/k10018/eq-archaeal-function-20261106/results/c18_uncalled_candidate_genes.tsv \
  --out-dir "${BASE}/results" \
  --tmp-dir "${TMPDIR}" \
  --threads 64

rm -rf "${TMPDIR}"
echo "[$(date -Is)] Completed exp1_metabolic Stage 4 Candidate Deficiency Search (v3 PASSED)"
