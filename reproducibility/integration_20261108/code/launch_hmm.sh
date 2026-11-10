#!/bin/bash
set -euo pipefail
cd /scratch/project/k10018/eq-proteome-integration-20261106
# Preparation and KO analysis completed with exit 0; Slurm has purged their controller records.
test -s results/protein_source.json
test -s results/paired_guild_tests_ko.tsv
pilot=$(sbatch --parsable --array=0 code/hmm.sbatch)
array=$(sbatch --parsable --dependency=afterok:$pilot --array=1-49%500 code/hmm.sbatch)
summary=$(sbatch --parsable --dependency=afterok:$pilot:$array code/hmm_summary.sbatch)
spatial=$(sbatch --parsable code/spatial_ko.sbatch)
spatialh=$(sbatch --parsable --dependency=afterok:$summary code/spatial_hmm.sbatch)
printf 'prepare\t15838457\nknown_panels_and_ecology\t15838458\ncandidate_review\t15838459\nhmm_pilot\t%s\nhmm_array\t%s\nhmm_summary_and_ecology\t%s\nspatial_sensitivity_ko\t%s\nspatial_sensitivity_hmm\t%s\n' "$pilot" "$array" "$summary" "$spatial" "$spatialh" | tee results/jobs.tsv
