#!/bin/bash
set -euo pipefail
cd /scratch/project/k10018/eq-proteome-integration-20261106
prep=$(sbatch --parsable code/prepare.sbatch)
known=$(sbatch --parsable code/known.sbatch)
candidates=$(sbatch --parsable code/candidates.sbatch)
pilot=$(sbatch --parsable --dependency=afterok:$prep --array=0 code/hmm.sbatch)
hmm=$(sbatch --parsable --dependency=afterok:$pilot --array=1-49%500 code/hmm.sbatch)
summary=$(sbatch --parsable --dependency=afterok:$pilot:$hmm code/hmm_summary.sbatch)
printf 'prepare\t%s\nknown_panels_and_ecology\t%s\ncandidate_review\t%s\nhmm_pilot\t%s\nhmm_array\t%s\nhmm_summary_and_ecology\t%s\n' "$prep" "$known" "$candidates" "$pilot" "$hmm" "$summary" | tee results/jobs.tsv
