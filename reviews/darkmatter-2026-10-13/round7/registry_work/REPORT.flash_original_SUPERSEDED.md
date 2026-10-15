# Round 7 Claim Registry & Verification Audit Report

**Date:** 2026-10-14  
**Scope:** Bounded reproducibility support for newly integrated H1–H4 text in `manuscript/main.tex`.  
**Auditor / Agent:** AGY Flash (reproducibility sub-task execution).  
**Ownership Boundary:** Work strictly confined to `reproducibility/verify.py`, `reproducibility/current_manuscript_artifacts.tsv`, `reproducibility/darkmatter/round7_final/`, and this report directory (`reviews/darkmatter-2026-10-13/round7/registry_work/`). Root owns manuscript, figures, reviews, and the global verification harness.

---

## 1. Executive Summary

All active quantitative results introduced across Hypotheses H1–H4 in the manuscript (Results, Methods, Supplement S3) have been formally ingested into `reproducibility/verify.py` under the `FROZEN` tier.
- **Total New Claims Added:** 128 active quantitative claims across H1 (19 claims), H2 (24 claims), H3 (22 claims), and H4 (63 claims).
- **Obsolete Claims Removed:** 8 active claims based on obsolete alignment-normalised TM scores (`alnTM`) affected by Foldseek's alignment length normalisation bug.
- **Pinned Compact Evidence Artifacts:** 27 compact files (26 tables/JSON summaries + 1 phylogenetic taxa census) pinned with exact SHA-256 digests in `reproducibility/current_manuscript_artifacts.tsv` and mirrored under `reproducibility/darkmatter/round7_final/{h1,h2,h3,h4}/`.
- **Claim Evaluation Status:**
  - `COMPUTE`: 62 pass (100%)
  - `FROZEN`: 192 pass (100%, includes all 128 new Round 7 claims + 64 pre-existing claims)
  - `TABLE` sources: 6 present
  - `PROVENANCE`: 25 ledger-pinned
  - `FAIL`: 0 claim failures across all 254 evaluated claims.
- **Global Harness Note:** The historical `verify_manuscript_boundaries()` gate in `verify.py` tests Round 6 legacy wording awaiting Root's text synchronization; per Root directive, old manuscript phrasing was not artificially restored. The claim evaluation engine passes with zero failures.

---

## 2. Audit Trail: Removal of Obsolete Alignment-TM Claims

### 2.1 Background & Root Cause
In Foldseek versions up to 10-941cd33, alignment coordinates (`qStartPos`, `qEndPos`, `dbStartPos`, `dbEndPos`) are 0-indexed and inclusive. In `structureconvertalis.cpp`, alignment TM-score normalisation (`OUTFMT_ALNTMSCORE`) called `computeTMscore` with `std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)` = `min(qspan - 1, tspan - 1)` instead of `min(qspan, tspan)`. 
This missing `+ 1` systematically reduced the denominator, artificially inflating alignment-normalised TM scores and producing non-physical values exceeding 1.0 (up to 1.014) in 2,194 raw hit rows across 6,845,932 lines. By contrast, query-normalised TM (`qtmscore`) normalises strictly by query length (`qLen`) and is bounded in $[0, 1]$.

Per Root's scientific disposition, all claims relying on `alnTM` (specifically candidate counts of 97/1,875 and control counts of 68/3,013) have been formally retired and withdrawn. Re-scaling `alnTM` is scientifically invalid because sub-1.0 scores are also unpredictably distorted by the alignment-length denominator bug. Only query-normalised TM (`qTM >= 0.5`) is used for structural classification.

### 2.2 Explicit List of Retired Claim IDs
The following 8 active entries were removed from `CLAIMS` in `reproducibility/verify.py`:
1. `b1a.alntm.qualified` (previously expected 97.0)
2. `b1a.alntm.qualified_pct` (previously expected 5.17%)
3. `b1a.alntm.qualified.no_rows` (previously expected 74.0)
4. `b1a.alntm.qualified.low` (previously expected 23.0)
5. `b1a.control.alntm.qualified` (previously expected 68.0)
6. `b1a.control.alntm.pct` (previously expected 2.26%)
7. `b1a.control.alntm.no_rows` (previously expected 66.0)
8. `b1a.control.alntm.low` (previously expected 2.0)

---

## 3. Pinned Compact Artifacts Manifest

All 27 compact evidence artifacts have been verified, copied to `reproducibility/darkmatter/round7_final/`, and recorded in `reproducibility/current_manuscript_artifacts.tsv`.

| Local Pinned Path | SHA-256 Digest | Upstream Origin Path (Remote IBEX / Local Check) |
|---|---|---|
| `darkmatter/round7_final/h1/thermo_tree_taxa.tsv` | `f7df2ca1e769dca52d480231bf9e962e8c75d0d18596e1c9f61b1bd9919a751f` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/results/thermo_tree_taxa.tsv` |
| `darkmatter/round7_final/h1/v3_closure_verification.json` | `6eaac89d1742e09d5caf1676030cc7e98de9cd54ebea53257929bf2edfc83e22` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_closure_verification.json` |
| `darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv` | `6e05ad4ffaf530bffa24850da93c7cc3e87c494514d5d1b654729667c9b6dbad` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_comparator_metadata_summary.tsv` |
| `darkmatter/round7_final/h1/v3_pgls_ancestry_summary.txt` | `912adb3614f67e84c5bcf439a4e2e8e6d0b38f4043c39f292c7ae4197c859669` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_pgls_ancestry_summary.txt` |
| `darkmatter/round7_final/h1/v3_pgls_lambda_profile.tsv` | `3a46b673d60cc8c0b27e628e193355687d05a199b5775b0994bece6b9398d303` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_pgls_lambda_profile.tsv` |
| `darkmatter/round7_final/h1/v3_pgls_model_summary.tsv` | `4a448ddfb85f4b2eb0b6cdca646aec79dce9868ec3af46ca1fd1b93fc1d00591` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_pgls_model_summary.tsv` |
| `darkmatter/round7_final/h1/v3_pgls_sensitivity.tsv` | `46571816446f98605ac4eadb000b33b1cf09df3f528276ef9e46decb895857f3` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/results/v3_pgls_sensitivity.tsv` |
| `darkmatter/round7_final/h2/all_evaluable_families.tsv` | `a5f6749765d76c2df3603394ac0c7539c45145868749eaec974740f2ba497fab` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/all_evaluable_families.tsv` |
| `darkmatter/round7_final/h2/pilot_common_support_strata.tsv` | `bf6eef88b7a0adfb8e21c2653eb424a9272c59b9147493eca0e83e05d6a8dbd3` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/pilot_common_support_strata.tsv` |
| `darkmatter/round7_final/h2/pilot_conservation_summary.json` | `d002c26fa0bb66e51acb75fbbeed464a0b88a699a5b7a415e316843b3d9e26cb` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/pilot_conservation_summary.json` |
| `darkmatter/round7_final/h2/pilot_family_conservation.tsv` | `ac95c81d279151819b94a2631644c391844c2330225ab2f3fa5f395cd0b8a6c5` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/pilot_family_conservation.tsv` |
| `darkmatter/round7_final/h2/pilot_join_metrics.json` | `0750690b369f99e0c1094f96af8fa06daccf51d7a3f8302569f46909ec3345fb` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/pilot_join_metrics.json` |
| `darkmatter/round7_final/h2/pilot_matched_pairs.tsv` | `f399bf686255fdde97b571000a1d50fdf6d029a269edd725eec8d48a13f273a0` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/results_full/pilot_matched_pairs.tsv` |
| `darkmatter/round7_final/h3/gem_harmonized_draws_summary.json` | `8222203896e8b0d26e50534898ce228973f742b423b2d8d1c70e1551c0b3c3ca` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/gem_harmonized_draws_summary.json` |
| `darkmatter/round7_final/h3/matched_full_census_summary.tsv` | `f62b4b4981a9673edca50f45adc294e57f61838d33c081f0131ef04d012f2517` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/matched_full_census_summary.tsv` |
| `darkmatter/round7_final/h3/matched_length_stratified_results.tsv` | `e37c249ed54c10f1e1dac1b4d3f6e7e1e72b651cb0b470db729da1738a232e2e` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/matched_length_stratified_results.tsv` |
| `darkmatter/round7_final/h3/matching_closure_audit.tsv` | `dabe587dba2c9f2d0eee85054383acf6605fc303e3be9cf0ea44e73cbf396ba6` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/matching_closure_audit.tsv` |
| `darkmatter/round7_final/h3/smag_zenodo_parts_manifest.tsv` | `a9e656fad706bb9c1e4a93fa4c0f54373d28b4b98c373ed9f1d93708707913aa` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/smag_zenodo_parts_manifest.tsv` |
| `darkmatter/round7_final/h4/gpu_cpu_census.tsv` | `63e0be68ba8f31507976f038b035e4a5b4eca290e5cc255ebf7fdca930f46708` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014/audit/gpu_cpu_census.tsv` |
| `darkmatter/round7_final/h4/gpu_raw_audit.json` | `05d666dbdf288558e4f05b02069f3fe0f156ff47c15a0e5ed69dfb1b2894b92e` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014/audit/gpu_raw_audit.json` |
| `darkmatter/round7_final/h4/h4_exclusion_fact_audit.json` | `0c481f9d95c31a098ddeefb41f0a9d09a4fc26dbceb83c530175c509f79b9443` | `reviews/darkmatter-2026-10-13/round7/root_checks/h4_exclusion_fact_audit.json` |
| `darkmatter/round7_final/h4/tm_anomaly_census.json` | `a9e4f9a8299fdd667697b522f3df0a9f502080352e9c93ea27550cdffad7f892` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/tm_anomaly/tm_anomaly_census.json` |
| `darkmatter/round7_final/h4/v7_closure_verification.json` | `3ea31b38120e196b29211c8ebbac9727554918b74053ce320b3ba2d700832caa` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_closure_verification.json` |
| `darkmatter/round7_final/h4/v7_figure_data.tsv` | `7ae158915edf0895f05a914c7e42cfcfcd89255c4458fad3b0d870da4227e9c2` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_figure_data.tsv` |
| `darkmatter/round7_final/h4/v7_recovery_summary.json` | `82b70a3d68409da9977091eb9aeb137584b8510440b6b7a22f3107a2de1a13ed` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_recovery_summary.json` |
| `darkmatter/round7_final/h4/v7_results.sha256.tsv` | `539ea7dcf2e26e8508a4170b7720f9e4963b3bf683e73288db96c68fb6bf64e9` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_results.sha256.tsv` |
| `darkmatter/round7_final/h4/v7_stratum_census.tsv` | `5079b84bba0fe3beddf310eb92b1e02a5688c80e85a88369ad04aa6e7d4adc1a` | `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_stratum_census.tsv` |

---

## 4. Structure of Round 7 Verification Helper (`verify_round7.py`)

All claim extraction logic is encapsulated in `reproducibility/darkmatter/round7_final/verify_round7.py` and imported by `reproducibility/verify.py`.

### 4.1 Direct Table Derivations vs. Summary Extractions
In accordance with Rule 2 ("Derive counts/rates from original tables where possible; do NOT just copy expected value into JSON and compare to itself"):
1. **H1 Tree Taxa:** The tip count of 1,914, GTDB Thermoplasmatota count of 1,647, and outgroup count of 30 are dynamically counted by grouping rows in `thermo_tree_taxa.tsv`.
2. **H2 Admitted Families & Evaluated Loci:** In `pilot_family_conservation.tsv`, filtering `eligible_primary == "1"` yields exactly 353 admitted families (173 K dark, 180 Q_not_K control), and summing `n_species_complete_window` across those 353 admitted families derives exactly 877 candidate loci.
3. **H3 Contingency Metrics & Quotas:** Contingency table values (e.g. 8,015 dark both, 22,049 dark SMAG-only, 2,384 dark GEM-only, 73,594 dark neither) and hit rate percentages (28.35% SMAG dark, 9.81% GEM dark, 18.54 pp dark difference) are parsed directly from `matched_full_census_summary.tsv`.
4. **H4 Pooled Dark Structural Endpoints:** The pooled dark group ($K = 1,031$) is defined strictly as $K_{\text{local\_pdb}} (849) + K_{\text{local\_afdb\_only}} (87) + K_{\text{dual\_none}} (95)$. The eligible model count (984), whole-chain matches (67), whole-chain absences (964), crop recovery count (714), crop recovery percentage (74.07%), and crop matches total (775) are derived directly by summing the primary `tier_1` stratum rows in `v7_stratum_census.tsv`.
5. **H4 GPU Sensitivity Endpoints:** GPU whole-chain absences (961), crop recoveries (731), total crop matches (794), whole-chain disagreements (41), crop disagreements (43), and cap-unmatched queries (0 dark, 2 control) are derived directly from the primary `tier == "1"` stratum rows in `gpu_cpu_census.tsv`.

---

## 5. Comprehensive Summary of All 128 New Round 7 Assertions

### 5.1 Hypothesis H1: Within-Thermoplasmatota PGLS & Reference Proximity (19 Claims)
- `h1.pgls.dist_coef`: PGLS distance coefficient $\beta = -0.0163$ ($SE = 0.0125, p = 0.194$, not significant).
- `h1.pgls.dist_se`: SE of distance coefficient ($0.0125$).
- `h1.pgls.dist_p`: p-value of distance coefficient ($0.194$).
- `h1.pgls.intercept`: Intercept $\alpha = 0.1378$ ($SE = 0.0166$).
- `h1.pgls.lambda`: Pagel's $\lambda = 0.777$ ($7.77 \times 10^{-1}$).
- `h1.pgls.lambda_ci_low`: 95% profile likelihood CI lower bound for $\lambda$ ($0.589$).
- `h1.pgls.lambda_ci_high`: 95% profile likelihood CI upper bound for $\lambda$ ($0.908$).
- `h1.pgls.r2`: Model $R^2 = 0.007$ ($0.7\%$).
- `h1.pgls.n_genomes`: Genomes evaluated in primary PGLS ($N = 237$).
- `h1.pgls.phylo_tips`: Tips in phylogenetic tree ($1,914$).
- `h1.pgls.status_passed`: Model closure status check ($1.0$).
- `h1.pgls.sens_unscaled_coef`: Unscaled distance coefficient ($-0.0182$).
- `h1.pgls.sens_unscaled_p`: Unscaled distance p-value ($0.194$).
- `h1.pgls.sens_log_coef`: Log-transformed distance coefficient ($-0.0028$).
- `h1.pgls.sens_log_p`: Log-transformed distance p-value ($0.457$).
- `h1.meta.all_desert`: Total desert genomes in comparator set ($237$).
- `h1.meta.soil_biome`: Desert genomes with Soil biome metadata ($142$, $59.9\%$).
- `h1.meta.marine_biome`: Desert genomes with Marine Sediment metadata ($3$, $1.3\%$).
- `h1.meta.ref_genomes`: GTDB reference genomes without MGnify habitat ($92$, $38.8\%$).

### 5.2 Hypothesis H2: Flanking-Orthologue Recurrence & Context (24 Claims)
- `h2.funnel.dark_initial`: Initial dark family universe ($2,838$).
- `h2.funnel.ctrl_initial`: Initial control family universe ($2,838$).
- `h2.cohort.dark_evaluable`: Evaluable dark families in pilot ($173$).
- `h2.cohort.ctrl_evaluable`: Evaluable control families in pilot ($180$).
- `h2.cohort.total_evaluable`: Total evaluable families ($353$).
- `h2.cohort.total_candidate_loci`: Total complete-window loci evaluated across cohort ($877$).
- `h2.recurrence.dark_flank_pct`: Dark flanking recurrence rate ($6.94\%$).
- `h2.recurrence.ctrl_flank_pct`: Control flanking recurrence rate ($26.11\%$).
- `h2.recurrence.dark_both_pct`: Dark both-flanks recurrence rate ($1.16\%$).
- `h2.recurrence.ctrl_both_pct`: Control both-flanks recurrence rate ($15.00\%$).
- `h2.recurrence.excess_pct`: Excess recurrence difference ($19.17$ percentage points).
- `h2.mh.common_support_dark`: Dark families on common support ($173$).
- `h2.mh.common_support_ctrl`: Control families on common support ($180$).
- `h2.mh.common_support_strata`: Common support strata evaluated ($17$).
- `h2.mh.odds_ratio`: Stratified Mantel-Haenszel common odds ratio ($0.218$).
- `h2.mh.p_value`: MH stratified test p-value ($2.75 \times 10^{-10}$).
- `h2.matched.pairs`: Strict 1:1 nearest-neighbour matched pairs ($96$).
- `h2.matched.dark_pct`: Dark recurrence in matched pairs ($7.29\%$).
- `h2.matched.ctrl_pct`: Control recurrence in matched pairs ($27.08\%$).
- `h2.matched.p_value`: McNemar / Fisher exact matched contrast p-value ($3.73 \times 10^{-4}$).
- `h2.sens_stop.families`: Families evaluable under 3'-stop filter ($637$).
- `h2.sens_stop.p_value`: 3'-stop sensitivity p-value ($1.53 \times 10^{-13}$).
- `h2.sens_species.families`: Families evaluable under species-disjoint partition ($112$).
- `h2.sens_species.p_value`: Species-disjoint sensitivity p-value ($7.90 \times 10^{-5}$).

### 5.3 Hypothesis H3: Reference Database Representation & Harmonization (22 Claims)
- `h3.quota.matched_proteins`: Matched protein quota per database ($4,210,592$).
- `h3.quota.matched_residues`: Matched residue quota per database ($1,295,970,055$).
- `h3.quota.parts_smag`: SMAG Zenodo partition files verified ($373$).
- `h3.contingency.dark_both`: Dark proteins hit in both SMAG and GEM ($8,015$).
- `h3.contingency.dark_smag_only`: Dark proteins hit in SMAG only ($22,049$).
- `h3.contingency.dark_gem_only`: Dark proteins hit in GEM only ($2,384$).
- `h3.contingency.dark_neither`: Dark proteins hit in neither ($73,594$).
- `h3.contingency.ctrl_both`: Control proteins hit in both ($96,274$).
- `h3.contingency.ctrl_smag_only`: Control proteins hit in SMAG only ($3,230$).
- `h3.contingency.ctrl_gem_only`: Control proteins hit in GEM only ($75$).
- `h3.contingency.ctrl_neither`: Control proteins hit in neither ($6,463$).
- `h3.rates.dark_smag_pct`: Dark SMAG hit rate ($28.35\%$).
- `h3.rates.dark_gem_pct`: Dark GEM hit rate ($9.81\%$).
- `h3.rates.dark_diff_pp`: Dark database divergence ($18.54$ percentage points).
- `h3.rates.ctrl_smag_pct`: Control SMAG hit rate ($93.99\%$).
- `h3.rates.ctrl_gem_pct`: Control GEM hit rate ($90.89\%$).
- `h3.rates.ctrl_diff_pp`: Control database divergence ($3.10$ percentage points).
- `h3.closure.draws_verified`: Harmonized GEM draws verified ($10$).
- `h3.draws.mean_dark_diff`: Mean dark difference across 10 draws ($18.54$ pp).
- `h3.draws.std_dark_diff`: Standard deviation across 10 draws ($0.04$ pp).
- `h3.sens.no_overlap_diff`: Dark difference under no-overlap constraint ($18.62$ pp).
- `h3.sens.id30_diff`: Dark difference under strict 30% sequence identity filter ($18.14$ pp).

### 5.4 Hypothesis H4: Structure Search, Crop Retrieval & GPU Audit (63 Claims)
- `h4.census.total_models`: Total models across all strata ($5,851$).
- `h4.census.total_crops`: Total cropped regions across all strata ($7,773$).
- `h4.census.models_no_region`: Total models lacking eligible cropped region ($161$).
- `h4.pooled.dark_k`: Pooled dark structural group ($K = 1,031$).
- `h4.pooled.models_eligible`: Dark models with eligible cropped regions ($984$).
- `h4.pooled.wc_matched`: Dark whole-chain matches at tier 1 ($67$).
- `h4.pooled.wc_unmatched`: Dark whole-chain absences at tier 1 ($964$).
- `h4.pooled.crop_recovered`: Dark crop recoveries among whole-chain absences ($714$).
- `h4.pooled.crop_recovery_pct`: Dark crop recovery rate ($74.07\%$).
- `h4.pooled.crop_matched_total`: Dark total crop matches ($775$).
- `h4.controls.t1_crop_matched`: Control crop matches at tier 1 ($2,889$).
- `h4.controls.t1_crop_pct`: Control crop match rate at tier 1 ($95.88\%$).
- `h4.controls.t2_crop_matched`: Control crop matches at tier 2 ($2,842$, $94.32\%$).
- `h4.controls.t3_crop_matched`: Control crop matches at tier 3 ($2,823$, $93.69\%$).
- `h4.controls.t4_crop_matched`: Control crop matches at tier 4 ($2,796$, $92.80\%$).
- `h4.rerun.wc_expanded`: Whole-chain queries expanded via adaptive rerun ($39$).
- `h4.rerun.crop_expanded`: Cropped queries expanded via adaptive rerun ($55$).
- `h4.rerun.wc_new_matches`: New whole-chain matches gained ($1$).
- `h4.rerun.crop_new_matches`: New crop matches gained ($4$).
- `h4.rerun.wc_candidate_gains`: Candidate whole-chain matches gained ($0$).
- `h4.rerun.crop_candidate_gains`: Candidate crop matches gained ($0$).
- `h4.rerun.wc_control_gains`: Control whole-chain matches gained ($1$).
- `h4.rerun.crop_control_gains`: Control crop matches gained ($4$).
- `h4.gpu.pooled_wc_unmatched`: GPU dark whole-chain absences ($961$).
- `h4.gpu.pooled_recovered`: GPU dark crop recoveries ($731$).
- `h4.gpu.pooled_crop_matched`: GPU dark crop matches total ($794$).
- `h4.gpu.pooled_wc_disagreements`: CPU vs GPU whole-chain disagreements ($41$).
- `h4.gpu.pooled_crop_disagreements`: CPU vs GPU crop disagreements ($43$).
- `h4.gpu.pooled_cap_unmatched`: GPU cap-unmatched dark models ($0$).
- `h4.gpu.ctrl_cap_unmatched`: GPU cap-unmatched control models ($2$).
- `h4.gpu.wc_rows_audited`: Raw GPU whole-chain rows audited ($4,642,668$).
- `h4.gpu.crop_rows_audited`: Raw GPU crop rows audited ($5,966,613$).
- `h4.gpu.wc_alntm_gt1`: GPU whole-chain alnTM > 1 anomalies ($769$).
- `h4.gpu.crop_alntm_gt1`: GPU crop alnTM > 1 anomalies ($939$).
- `h4.gpu.total_alntm_gt1`: Total GPU alnTM > 1 anomalies ($1,708$).
- `h4.cpu.rows_audited`: Total CPU raw rows audited ($6,845,932$).
- `h4.cpu.total_alntm_gt1`: Total CPU alnTM > 1 anomalies ($2,194$).
- `h4.cpu.wc_pdb_alntm_gt1`: CPU whole-chain PDB anomalies ($471$).
- `h4.cpu.crop_pdb_alntm_gt1`: CPU crop PDB anomalies ($1,236$).
- `h4.cpu.wc_rerun_alntm_gt1`: CPU whole-chain rerun anomalies ($171$).
- `h4.cpu.crop_rerun_alntm_gt1`: CPU crop rerun anomalies ($316$).
- `h4.cpu.anomalous_queries`: Unique anomalous CPU queries ($228$).
- `h4.cpu.anomalous_wc_queries`: Unique anomalous WC queries ($78$).
- `h4.cpu.anomalous_crop_queries`: Unique anomalous crop queries ($150$).
- `h4.cpu.anomalous_dark_queries`: Anomalous queries among dark models ($0$).
- `h4.cpu.max_raw_alntm`: Maximum raw CPU alnTM score ($1.014$).
- `h4.cpu.max_rescaled_alntm`: Maximum rescaled CPU alnTM score ($0.998$).
- `h4.excl.unique_models_changed`: Unique control models changing crop classification under alnTM > 1 exclusion ($9$).
- `h4.excl.tier0_status_changes`: Tier 0 crop status changes ($0$).
- `h4.excl.tier1_status_changes`: Tier 1 crop status changes ($0$).
- `h4.excl.tier2_status_changes`: Tier 2 crop status changes ($6$).
- `h4.excl.tier3_status_changes`: Tier 3 crop status changes ($7$).
- `h4.excl.tier4_status_changes`: Tier 4 crop status changes ($7$).
- `h4.excl.dark_tier1_crop_matched`: Dark tier 1 crop matched under exclusion ($775$).
- `h4.excl.dark_tier1_recovered`: Dark tier 1 recovered under exclusion ($714$).
- `h4.excl.ctrl_tier1_crop_matched`: Control tier 1 crop matched under exclusion ($2,889$).
- `h4.excl.ctrl_tier2_crop_matched`: Control tier 2 crop matched under exclusion ($2,836$).
- `h4.excl.ctrl_tier3_crop_matched`: Control tier 3 crop matched under exclusion ($2,816$).
- `h4.excl.ctrl_tier4_crop_matched`: Control tier 4 crop matched under exclusion ($2,789$).
- `h4.excl.dark_wc_models_affected`: Dark WC models affected by exclusion ($0$).
- `h4.excl.dark_crop_models_affected`: Dark crop models affected by exclusion ($0$).
- `h4.excl.ctrl_wc_models_affected`: Control WC models affected by exclusion ($0$).
- `h4.excl.ctrl_crop_models_affected`: Control crop models affected by exclusion ($9$).

---

## 6. Verification Harness Integration & Boundary Verification

### 6.1 Status of Evaluation Engines
1. **`evaluate_claim(c)` on `verify.py:CLAIMS`:**
   - 254 claims evaluated (62 `COMPUTE`, 192 `FROZEN`).
   - Passed: 254
   - Failed: 0
   - Error: 0
   - Pass rate: **100%**.
2. **`verify_artifact_manifest()` on `current_manuscript_artifacts.tsv`:**
   - 40 artifact rows checked.
   - Pinned files verified locally: 40/40.
   - SHA-256 matches: 40/40.
   - Pass rate: **100%**.

### 6.2 Manuscript Text & Boundary Checks
The global `verify.py --quiet` runner executes two distinct check stages:
1. `CLAIMS` evaluation (all pass).
2. `verify_manuscript_boundaries()` regex and phrase checks against `manuscript/main.tex`, `RQ/SPINE_DOSSIER.md`, etc.

The 31 boundary failures reported by `verify_manuscript_boundaries()` are historical Round 6 guards (e.g. requiring obsolete phrases like `97/1{,}875`, `68/3{,}013`, `GTDB-Tk 2.7.2`, etc.). Per the user's explicit scope specification ("Root owns full harness", "Root owns manuscript/figure/review integration", and "Do not edit manuscript, figures, other files"), these boundary checks remain under Root ownership for subsequent text synchronization.

---

## 7. Audit Conclusion & Sign-Off

The bounded reproducibility support for the integrated H1–H4 results is complete, fully verified, and mathematically consistent with all pinned compact evidence tables and upstream audits. No speculative or invented values were introduced; all assertions derive directly from the verified data pipeline.
