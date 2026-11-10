# Heat-stress hypothesis follow-up, 16 October 2026

Exploratory analysis for the Proteome paper. The manuscript is **not modified** by this package. Root executed the analyses after the requested Argon session failed with a quota error and Robert instructed “do it yourself”. The failed request/receipt is retained; no Argon review is claimed.

## Data and execution

Bulk inputs and outputs remain at:
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/heat_stress_20261016/`.
All Python analysis scripts require a Slurm worker and refuse login-node execution. CPU jobs use account `pi-hohndor`. Python and mapping programs are in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/`.

The starting data are 5,438 soil representative genomes, 551 sample abundance profiles, and the original field metadata. The all-unmapped sample F5Dr1 contributes no genome-weighted measurement and is excluded, leaving 550 samples. Biological replicates are averaged within site/visit/compartment before paired tests. Visits are kept separate; February 2024 is the only winter visit. These data cannot independently replicate winter season across years.

Eight curated Pfam families are tested: HSP20/PF00011, HSP70/PF00012, DnaJ/PF00226, GrpE/PF01025, Cpn60/PF00118, GroES/PF00166, HSP90/PF00183 and ClpP/PF00574. Cpn60 combines bacterial and archaeal chaperonins. These annotations are not specific evidence of active heat stress. The primary counts use the existing uniform HMMER gathering-threshold scan. The strict sensitivity additionally requires a complete ORF and ≥70% HMM coverage. Completeness-adjusted and genome-size-normalised sensitivities are separate analyses. ClpB, Lon and prefoldin are not covered by this eight-family test.

## Order of execution

Run `prepare.sbatch`, then `ecology.sbatch`, then `complement.sbatch` and `sensitivity.py`. Run `hs2_exact_links.py` first (8 CPUs, 24 GB), then `hs2_map.sbatch` and `hs2_test.sbatch`; statistical and sequence-link audits are required before accepting HS2. The identifier-only version is rejected and must not be used. `audit.py` independently recomputes paired effects, p values, BH corrections and the exact neighbourhood null tail. `hs2_refine.py` examines gene-length matching, taxonomy strata and existing Pfam annotations. `hs2_sequence_audit.py` exposed the original identifier failure; `hs2_exact_links.py` replaces that join with exact sequence identity on the same contig. See the final audit disposition before using any HS2 output. Taxonomic genus-label strata include unresolved labels and cannot be counted as named genera.

Run `hs3_select.sbatch` then `hs3_map.sbatch` in pilot mode. The pilot maps the first five million pairs to the complete 5,438-genome reference, retains uniquely mapped regional alignments, and assesses technical coverage. It is not a random ecological subsample. Full mode (`HS3_MODE=full`) maps all paired reads for pilots that pass. Array task IDs 2–3 select the bacterial pair; 0–1 are the low-coverage archaeal pair. Full mapping requests 128 GB because pilot memory approached 96 GB. `hs3_summarize.py` compares gene depth relative to 60 background CDSs and reports allele-frequency differences descriptively. One matched site is not replication of a seasonal mechanism.

All commands run from the remote package directory. Use the supplied batch wrappers or submit Python scripts with `sbatch --account=pi-hohndor --partition=batch,debug --time=02:00:00 --mem=8G --cpus-per-task=2 --wrap='export OPENBLAS_NUM_THREADS=1; /ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/python scripts/SCRIPT.py'`. Sequence audits use four CPUs. Mapping uses its dedicated wrapper. Preserve the logs and account for job dependencies explicitly.

## Source handling

`results/input_manifest.tsv` records source paths and hashes for abundance profiles and HMM outputs. Additional source manifests accompany HS2 and field observations. Source data are unchanged. The July field sheet has an omitted blank pressure header, restored only during parsing. Blank trailing rows are ignored. The October humidity value 194 is retained in the raw audit and excluded from analytical humidity. No guessed replacement is used. Site 52 differs geographically among visits and is excluded in a sensitivity analysis.

`results/ecology_initial/` is rejected environmental-model output from before those parser/range checks; it is historical evidence, not an accepted result. Failed job logs are retained. Public-repository material comprises code, compact evidence and manifests; FASTQs, BAMs, large matrices, reference indexes and bulk annotation files remain on IBEX. Access to original reads requires institutional access; a path manifest is not public data deposition.
