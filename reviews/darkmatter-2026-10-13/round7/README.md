# Full extension analyses and results review

Round 7 executes the full H2 eligible-family analysis, broadly sampled matched-reference H3 full-query searches, and H4 paired structural census. H1 already has a complete 237-genome PGLS and was audited for completeness. Gemini Flash 3.8 High implements; Gemini 4 Argon High independently reviews actual results, code, reproducibility and manuscript fit. Round 6 approved feasibility packages remain unchanged.

## Completed and reviewed

H1 and H2 have final Argon verdicts **SUPPORTED_WITH_BOUNDS**, with all specified findings closed. See [final disposition](argon_h1h2_review1/verdicts_final.json), [corrected review](argon_h1h2_review1/review.txt), [verification](argon_h1h2_review1/verify_h1_h2_facts.json), and [closure](argon_h1h2_review1/closure.txt).

- H1: Reuse the complete 237-genome analysis and six sensitivity specifications. The within-Thermoplasmatota conditional associations survive the tested specifications. High fitted lambda describes residual phylogenetic covariance; it does not quantify variance explained or establish gain/loss mechanisms. Summary tables are generated directly from authoritative source files after a transcription defect was identified and repaired.
- H2: Full 1,461-family census (570 K/891 QnotK), 7,604 candidate loci; 353 families satisfy the strict primary context requirement. Within-arm KO recurrence exceeds the specified positional expectation, but the adjusted between-arm contrast is not established (difference 0.035358, permutation p = 0.406859). Strong eligibility attrition, coarse taxonomic adjustment and cross-family dependence limit interpretation. Argon recommends brief main-text mentions and detailed supplementary analyses for H1/H2.

H2 final snapshot: `snapshots/h2_final/`; H1 authoritative results: `../round6/snapshots/h1_final/`. H2 primary Slurm 53738443 completed; replay 53740945 after minor code cleanup produced 14 byte-identical result files. [Replay evidence](root_checks/h2_review_replay_verification.json).

The first Argon draft itself contained factual errors. The original, [root fact audit](REVIEW1_FACT_AUDIT.md), corrected report, saved verification code and [correction ledger](argon_h1h2_review1/FACT_CORRECTIONS.txt) are retained. A favorable verdict alone is not treated as verification.

## Still running

H3: 1,500 SMAG genomes randomly selected from the complete 40,039-genome release; 2,055 GEM genomes selected per comparator draw. Primary GEM DNA downloaded and called with the same Pyrodigal settings planned for SMAG; all three GEM draws are now gene-called with the same Pyrodigal settings, alongside published-proteome caller sensitivities. Additional draw job53742584 completed with 4,220,573 and 4,337,716 unique canonical proteins for seeds20261015/16. Caller/overlap/identity and draw comparisons are queued. GEM raw FAA files, GFF coordinates and mappings retained and all 4,824,709 canonical protein IDs/sequences verified against the original pool. Full SMAG archive transfer is in progress. Full 212,084-query searches have not yet run.

H4: all 5,851 models cropped, 7,773 regions, 161 models without eligible regions; both full paired whole-chain/crop PDB searches completed, followed by signed adaptive searches of 39 whole chains and 55 crops at a larger hit cap. Full census tables are generated and table invariants pass, but independent raw audit detected an out-of-range TM score. Exact-version diagnosis and an exclusion sensitivity are in progress; no final census result is yet accepted.

Live job/session records: [STATE.json](STATE.json) and each analysis area's STATUS.json. Root operational changes are documented in [audit notes](ROOT_AUDIT_NOTES.md). Manuscript remains unchanged pending completion of results review.

## Reproducibility locations

Local: `reproducibility/darkmatter/round7_analyses/`.
IBEX: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/`.

Exact prompts, model selections and CLI receipts are retained under `prompts/` and `receipts/`. Large assemblies, translated proteomes, databases, raw searches and gene-coordinate files remain on IBEX with compact evidence locally. Previously documented public-release/deposition and historical raw-data access gaps remain unresolved by these extension analyses.
