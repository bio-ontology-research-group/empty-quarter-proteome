# Four viable dark-matter paper extensions

Gemini 4 Argon High has independently judged **H1, H2, H3 and H4 VIABLE**. Each first review requested revision; all four revised designs received an independent second review with no critical issues. Final documentation corrections are complete and independently confirmed by Argon; all 265 frozen package files passed local and IBEX verification. Nature Communications remains the manuscript target; this judgment establishes bounded scientific viability, not positive findings or publication readiness.

The work used agy Gemini Flash High for initial implementation, agy Gemini 4 Argon High for critical reviews and substantive repairs, and fresh Argon conversations for independent second reviews. Root independently checked numerical results, eligibility and mapping invariants, source provenance, and file hashes. Original results, failed runs, service errors, reviews, responses and replacement versions are retained.

| Hypothesis | Approved question | Executed scope | Main finding and boundary |
|---|---|---|---|
| H1: lineage | How does Thermoplasmatota dark fraction vary with external-reference AAI and shared-gene coverage after quality and phylogenetic adjustment? | All 237 genomes; branch-length PGLS, lambda profile and seven sensitivities; independent statsmodels verification | Shared-gene coverage has a negative association with dark fraction; AAI has a small positive association. These are conditional associations, not ancestral gain/loss or accessory expansion. |
| H2: neighbourhoods | Do low-query-TM and fold-matched dark families differ in flanking KO recurrence after accounting for eligible positional backgrounds and measured covariates? | 300 families, 1,518 candidate loci, 685 sequence-verified bins; 72 primary and 133 sensitivity-evaluable families | Both groups show recurrence relative to the chosen positional reference. The primary adjusted between-group contrast is not significant (difference −0.0648; permutation p=0.4907). Dependence/exchangeability remain assumptions with sensitivity checks. |
| H3: references | Does reference recovery differ between verified SMAG and GEM pools matched exactly on canonical sequence length, count and redundancy? | 212 complete SMAG genomes; 2,000 query proteins; three GEM draws plus overlap and identity sensitivities | Weighted dark-query recovery is approximately 10–11 percentage points higher for the sampled SMAG pool. This is conditional on the convenience pool, caller releases and composition, not a causal biome effect. |
| H4: structure | Does cropping confident regions recover PDB similarity differently for whole-chain-unmatched models with versus without baseline local matches? | All 5,851 models cropped (7,773 regions); paired whole/crop searches on a weighted 300-model pilot (413 regions) | Weighted recovery among fresh whole-chain-unmatched models is 83.0% versus 36.0%. Recovery is reference- and algorithm-dependent; nonretrieval does not establish a novel fold. |

## Approved designs and reviews

- H1: [dossier](../../../reproducibility/darkmatter/round6_hypotheses/h1_lineage/DOSSIER.md), [independent review](h1_review2/review.txt), [final confirmation](h1_review2/verdicts_final.json).
- H2: [dossier](../../../reproducibility/darkmatter/round6_hypotheses/h2_context/DOSSIER.md), [independent review](h2_review2/review.txt), [final confirmation](h2_review2/verdicts_final.json).
- H3: [dossier](../../../reproducibility/darkmatter/round6_hypotheses/h3_references/DOSSIER.md), [independent review](h3_review2/review.txt), [final confirmation](h3_review2/verdicts_final.json).
- H4: [dossier](../../../reproducibility/darkmatter/round6_hypotheses/h4_structure/DOSSIER.md), [independent review](h4_review2/review.txt), [final confirmation](h4_review2/verdicts_final.json).

Each analysis directory contains DESIGN.md, DOSSIER.md, INPUTS.tsv, scripts, compact results, logs and the response to the previous review. The final approved copies are frozen under snapshots/h1_final through snapshots/h4_final, each with a SHA-256 manifest. Prior versions are preserved separately. The [state record](STATE.json) records individual acceptance and completion state; the [review contract](REVIEW_CONTRACT.md) defines viability.

## IBEX data and reproducibility

The underlying inputs, extracted genomes, matched protein databases, cropped structures, raw searches and Slurm outputs are retained under:

`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/`

Subdirectories h1_lineage, h2_context, h3_references and h4_structure contain the corresponding analyses. The review package is mirrored under review_package/round6/. INPUTS.tsv and the result-specific manifests identify upstream data paths, versions, hashes and unavailable inputs. All 32 actual PDB reference components (4,530,919,251 bytes) passed the independent checksum audit; the historical checksum list's self-entry is excluded and its current digest pinned separately.

This is an internally reproducible feasibility package, not a claim that every original manuscript input is publicly accessible. Previously documented deposition/access gaps remain. The full H2 catalogue analysis, full H3 query/catalogue expansion and full H4 alignments are runnable proposed extensions, not completed census analyses. H4's full-crop run was stopped before full alignment after measuring runtime scaling; the corrected weighted pilot is the executed search evidence.

## Verification

- [H1 numerical and covariance checks](../../../reproducibility/darkmatter/round6_hypotheses/h1_lineage/results/v3_closure_verification.json).
- [H2 independent table audit](root_checks/h2_table_audit.json): eligibility before aggregation, unique species loci, sequence-verified bins, recurrence bounds.
- [H3 independent arithmetic audit](root_checks/h3_table_audit.json): ten weighted contrasts and standard errors; 212 complete archive members and rejected truncated member 213.
- [H4 independent invariants](root_checks/h4_pilot_invariants.json): 6,152 model/region assertions, seven sampling strata and finite-population weights.
- [PDB reference audit](root_checks/pdb_reference_verification.json).
- [Final local verification](root_checks/final_package_verification.json) and [IBEX verification](root_checks/ibex_final_package_verification.json): all 265 package files match approved snapshots and live analysis files; IBEX Slurm job 53716667 completed with exit 0. [Analysis job accounting](root_checks/final_slurm_accounting.tsv) confirms final analysis runs completed successfully.

The manuscript and bibliography remain at the round-six baseline; no exploratory pilot result has been silently promoted into the manuscript. Original workflow and acceptance criteria are retained in [WORKFLOW_ORIGINAL.md](WORKFLOW_ORIGINAL.md).
