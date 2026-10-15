# H1–H4 extension evidence and execution map

The manuscript source is `manuscript/main.tex`. The local final compact evidence is
`reproducibility/darkmatter/round7_final/`; its claims are FROZEN table extractions,
not local reruns of genome searches. Reviews, independent audit code, scheduler
receipts and author responses are under `reviews/darkmatter-2026-10-13/round7/`.

IBEX workbench root: `/ibex/scratch/projects/c2014/rob/metagenomics`.

| Analysis | Authoritative IBEX directory below the workbench | Reproduction evidence |
| --- | --- | --- |
| H1, phylogeny-conditioned Thermoplasmatota associations | `06_darkmatter/p6_hypotheses_clean/h1_lineage` | Final v3 scripts and results; `v3_input_file_manifest.tsv`; tree and comparator manifests; local `round6/snapshots/h1_final` |
| H2, flanking-orthologue recurrence | `06_darkmatter/p7_analyses_clean/h2_context` | `scripts/` and **`results_full/`**, executed input manifest, output checksums; local `round7/snapshots/h2_final`; exact analysis replay job53740945 |
| H3, exact-length matched references | `06_darkmatter/p7_analyses_clean/h3_references` | Genome/part manifests, sampled IDs, raw calls and sequence maps, matched FASTAs, DIAMOND indexes and hits, completion signatures; `final_raw_audit/` independently rereads raw inputs/outputs |
| H4, paired whole chains and confidence regions | `06_darkmatter/p7_analyses_clean/h4_structure` | `results/`, original raw hits, `scripts_recovery_v2/` executed postprocessing and `scripts_recovery_v3/` repaired census verifier; frozen database/query hash inventories |
| H4, GPU retrieval sensitivity | `06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014` | `gpu_results/{wc,crop}` preserves ORIX raw tables, executed scripts, GPU monitor records and worker hashes; `audit/` independently validates rows and compares CPU/GPU qTM outcomes |

ORIX source directory:
`/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014`.
GPU array315233 used freecycle H200 allocations; whole-chain execution completed
in restart2 and a later restart verified completed hashes. IBEX preparation53760480
staged the frozen target and query databases. The GPU run changes the prefilter
algorithm, and its stricter control exclusions retain finite retrieval-depth limits.

All 373 SMAG archive parts and the sampled GEM genome/call artifacts remain on IBEX.
H3 matching/search completion signatures validate cached outputs against input
hashes and command parameters. Original failed-attempt code and logs are preserved;
use the stated final versions rather than inferring success from filenames.

H4's original strict verifier rejects out-of-range alignment-normalised scores.
The source audit identified an inclusive-span normalisation defect in Foldseek10.
Those original raw values are preserved. The accepted endpoint is qTM, and the
initial 1,031-model dark group pools all historical alnTM subdivisions. The independent
raw validation explicitly records alnTM anomalies rather than loosening their bound.
See `round7/root_checks/H4_SCIENTIFIC_DISPOSITION.md` in the review tree.

CPU analysis, parsing, compilation, extraction and bulk checksums must run inside
IBEX Slurm allocations. Use account `pi-hohndor` and invoke commands through
`require_slurm_worker.sh`; a script's `.sbatch` suffix is not an allocation.
ORIX allocations are for actual GPU searches. Old executed scripts are historical
artifacts; wrap future invocations with the worker guard rather than running them
on a login node. The guard was tested to reject login execution and admit a genuine
allocated worker.

The whole-paper raw-read access/deposition gaps remain separately documented in
`reproducibility/darkmatter/README.md`. Complete extension-specific evidence does
not close those earlier raw-data gaps or make the paper submission-ready.
