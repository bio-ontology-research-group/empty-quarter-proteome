# DOCUMENTATION CORRECTIONS & EPHEMERAL AUDIT LOG

**Target Area:** `reproducibility/darkmatter/round7_analyses/h2_context/` & Remote Mirror `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/`  
**Date:** 2026-10-14  
**Author:** Gemini 3.8 Flash (acting on behalf of Robert Hoehndorf)  
**Reference:** `reviews/darkmatter-2026-10-13/round7/ROOT_AUDIT_NOTES.md` §H1 & §H2  

---

## 1. Summary of Documentation Corrections

In accordance with root audit instructions prior to independent evaluation by Argon, documentation in the H2 area was refined to ensure exact epistemic bounding, remove overstatements, and eliminate misleading inferential phrasing without altering completed numerical analyses or raw outputs.

### 1.1 Corrections to `H1_COMPLETENESS.md`
1. **Authoritative Final Snapshot:**
   - *Previous Text:* Cited `reviews/darkmatter-2026-10-13/round6/snapshots/h1_v3/` and working directory `reproducibility/darkmatter/round6_hypotheses/h1_lineage/`.
   - *Correction:* Updated Evaluation Scope and provenance references to cite `reviews/darkmatter-2026-10-13/round6/snapshots/h1_final/` as the immutable final authority (where v3 code remains identical and final documentation was approved, superseding `h1_v3`).
2. **Model Counting Precision (Seven Specifications / Six Sensitivities):**
   - *Previous Text:* Described the models as "primary 237-genome PGLS model and all seven pre-specified sensitivity analyses ($S1$–$S7$)".
   - *Correction:* Corrected to specify seven model specifications in total: primary model $S1$ plus six sensitivity models $S2$–$S7$.
3. **Verification Summary Characterization:**
   - *Previous Text:* Stated that "All 11 closure assertions in `v3_closure_verification.json` evaluate to `true`".
   - *Correction:* Replaced with accurate description of `v3_closure_verification.json` as a heterogeneous verification summary (validating sample sizes, excluded genome lists, machine-precision float differences vs. `statsmodels`, ML $\lambda$ profiling, and variance bounds, recording top-level `"status": "PASSED_ALL_ASSERTIONS"` alongside two Boolean integrity checks).
4. **Phylogenetic Tree Topology & Branch-Length Uncertainty:**
   - *Previous Text:* Suggested that rooting invariance ($|\Delta| = 0.0$) and $S1$ vs. $S2$ scaling resolved phylogenetic tree uncertainty.
   - *Correction:* Clarified that mathematical equivalence under midpoint vs. outgroup rooting and Brownian scaling checks on the fixed FastTree VCV do *not* establish that underlying phylogenetic tree topology or branch-length uncertainty of the GTDB r232 `ar53` FastTree is resolved; topology and branch lengths remain fixed conditioning assumptions.
5. **Gene-Caller Discrepancy & Residual Confounding:**
   - *Previous Text:* Claimed gene-caller discrepancy was "fully resolved" by $S3$ and $S7$.
   - *Correction:* Clarified that shared-catalogue denominator ($S3$) and complete-gene proxy ($S7$) sensitivities demonstrate numerical stability across specifications, but do *not* fully resolve or eliminate gene-caller confounding, which remains an intrinsic limitation of multi-caller metagenomic workflows.
6. **Bounded Tone & Qualified Language:**
   - *Previous Text:* Used unqualified assertions such as "Zero Genuinely Missing Analyses", "meets all criteria for scientific rigor", and "ready for publication".
   - *Correction:* Replaced with bounded completeness within the pre-declared model grid and explicit observational data frame, framed strictly as pending independent critical evaluation by Argon.

---

### 1.2 Corrections to `RESULTS.md`
1. **Observational Nature of Sequence Length Differences:**
   - *Previous Text:* Stated that Tier 1 $K$ families have longer focal proteins ($320.8$ aa vs. $233.1$ aa) which "elevate baseline detection opportunity".
   - *Correction:* Explicitly clarified that differences in focal protein length and flanking KO annotation opportunity are observational baseline covariate imbalances/associations that correlate with detection opportunity, avoiding calling them demonstrated confounding or causal mechanisms.
2. **Identifying Assumption of Stratified Permutation Tests:**
   - *Previous Text:* Described label permutation without highlighting the underlying identifying assumption.
   - *Correction:* Explicitly stated that permuting arm labels strictly within strata relies on the identifying assumption of conditional exchangeability given observed strata (`species_breadth_bin`, `phylum_group`, `ko_opportunity_bin`). Emphasized that this is a modeling assumption for stratified conditional inference, not an empirical proof that families or host genomes are mutually independent; cross-family and phylogenetic dependencies are not eliminated by permutation.
3. **Refined Synthesis & Conclusion:**
   - *Previous Text:* Stated that the hypothesis of differential recurrence was "refuted".
   - *Correction:* Stated that the hypothesis of differential neighbourhood recurrence between novel-fold candidates and fold-matched sequence-dark families is not supported (consistent with the adjusted null under the specified stratification and matching models), conditionally bounded by the permutation exchangeability assumption.

---

## 2. Script and Workflow Provenance Verification

1. **Audit of `run_pilot.sbatch` on IBEX Remote (`p7`):**
   - Inspected `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/scripts/run_pilot.sbatch`.
   - Verified that `WORKDIR` is set to `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context` and outputs redirect to `p7_analyses_clean/h2_context/logs/full_%j.out` and `results_full/`.
   - Confirmed that `run_pilot.sbatch` does *not* point to `p6` and poses zero risk of overwriting `06_darkmatter/p6_hypotheses_clean/` or `round6_hypotheses/` artifacts.
   - `run_full.sbatch` remains the authoritative production submission script used for Slurm Job `53738443`.
2. **Preservation of Raw Results and P6 Artifacts:**
   - All completed Slurm outputs, log files (`logs/full_53738443.out`, `logs/full_53738443.err`), intermediate extracted contig caches, and 13 pinned output files in `results/` remain intact and cryptographically verified.
   - All `p6` directories and snapshots remain untouched.
