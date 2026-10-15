# EVIDENCE DOSSIER: Flanking Gene-Neighbourhood Recurrence of Low-Query-TM vs. Fold-Matched Sequence-Dark Structural Families (v3 Protocol)

**Hypothesis 2:** Bounded flanking gene-neighbourhood KEGG KO recurrence of low-query-TM ($K$, $\text{qTM} < 0.5$ vs. PDB and AFDB) versus fold-matched sequence-dark ($Q \setminus K$, $\text{qTM} \geq 0.5$) monomer structural families across distinct refined-MAG dRep 95% ANI species clusters.  
**Author:** Gemini 4 Argon (acting on behalf of Robert Hoehndorf)  
**Date:** 2026-10-13  
**Status:** Version 3 (Empirical Pilot Executed under IBEX Slurm Jobs `53710523`, `53712511`, and `53713357`; Initial Diagnostic Run `53710190` Preserved)  
**Local Artifact Directory:** `reproducibility/darkmatter/round6_hypotheses/h2_context/`  
**Remote Workbench:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/`  
**Superseded Versions:** Version 1 (`v1_superseded/`, frozen in `reviews/darkmatter-2026-10-13/round6/snapshots/h2_initial/`) and Version 2 (`v2_superseded/`, frozen in `reviews/darkmatter-2026-10-13/round6/snapshots/h2_v2/`).

---

## 1. Executive Summary & Resolution of All v2 Review Findings (`C1`–`C5`)

Following the independent v2 audit (`reviews/darkmatter-2026-10-13/round6/h2_review1/review.txt`) and root audit notes (`reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`), Version 3 (v3) redesigns and re-executes the Hypothesis 2 pipeline on an expanded 300-family stratified cohort ($150 \, K$, $150 \, Q \setminus K$):

1. **`C1` Resolved — Locus-First Eligibility Gate Before KO Aggregation (`R_fam in [0, 1]` Guaranteed):**  
   In v2, KOs were aggregated across non-truncated flanks before checking locus window completeness, causing 4 families to exceed $R_{\text{fam}} > 1.0$ and 2 families to be false-positive $C_{\text{fam}} = 1$ calls, while `source_file_unavailable` and `coordinate_mismatch` flanks were left as callable negatives. In v3 (`scripts/02_extract_neighbourhoods.py` and `scripts/03_analyse_conservation.py`), every flank is classified into 7 mutually exclusive categories, locus eligibility is determined **prior** to any family-level aggregation, and primary numerators, denominators, and `ko_annotation_opportunity` are computed over the **exact same set** of Tier 1 eligible primary loci (`eligible_locus_primary == 1`, requiring all 4 flanks in `{"mapped_KO", "mapped_noKO"}`). All $R_{\text{fam}} \in [0, 1]$ across all 300 families (`max = 1.0000`), and all ineligible families ($S_F < 2$) are coded `-1` (`INELIGIBLE / MISSING`), never `0`.
2. **`C2` Resolved — Real Contig Lengths, Verified GTDB-Tk R232 Taxonomy (`100%` Coverage), and Converged Full-Rank Regressions:**  
   In v2, `contig_len` was `0` for 100% of rows (`534 / 534`) and 99.6% of rows had `taxonomy == ""` (`phylum = "unassigned"`), causing quasi-complete separation and logistic regression non-convergence (`converged = False`). In v3, `scripts/01_select_cohort.py` resolves species-level taxonomy via `Cdb.csv -> Wdb.csv -> 05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` (`5,438 / 5,438` species clusters and `1,461 / 1,461` multi-species families resolved; **0% unassigned**), and `scripts/02_extract_neighbourhoods.py` measures true base-pair contig lengths (`contig_len_bp > 0` for 100% of binned rows) and window completeness **before** selecting one deterministic primary locus per `(family, species)`. Multivariable OLS/HC3 and MLE/L2-penalized logistic regressions converge (`converged = True`, full column rank `9 / 9`, all standard errors finite and bounded `< 3.1`).
3. **`C3` Resolved — Separated Within-Arm Positional Enrichment vs. Between-Arm Covariate-Adjusted Contrast:**  
   In v2, the permutation null sampled edge-truncated loci/positions and compared the observed between-arm risk difference against the lower-prevalence random-position difference distribution, conflating within-arm synteny with between-arm differences. In v3 (`scripts/03_analyse_conservation.py`), inference is separated into:
   - **Part (a) Within-Arm Recurrence vs. Same-Contig Eligible-Position Reference ($B = 5{,}000$):** Pseudo-focal positions are drawn strictly from indices on the same contig whose $\pm 2$ flanks all satisfy the exact same Tier 1 complete-window gate. Both arms show strong within-arm recurrence enrichment over random eligible contig positions: Arm $K$ observed $\bar{C}_K^{\text{obs}} = 48.84\%$ (`21/43`) vs. positional expectation $\bar{\mu}_K^{\text{pos}} = 12.86\%$ (excess $+35.97$ percentage points, $p = 0.0002$); Arm $Q \setminus K$ observed $\bar{C}_{Q\setminus K}^{\text{obs}} = 31.03\%$ (`9/29`) vs. positional expectation $\bar{\mu}_{Q\setminus K}^{\text{pos}} = 4.96\%$ (excess $+26.08$ percentage points, $p = 0.0002$).
   - **Part (b) Between-Arm Contrast ($K$ vs. $Q \setminus K$):** Tested on family-level background-adjusted excess $E_{\text{fam}} = C_{\text{fam}}^{\text{obs}} - \mu_{\text{fam}}^{\text{pos}}$ and observed $C_{\text{fam}}^{\text{obs}}$ via (1) primary covariate-stratified label permutation across `(species_breadth_bin, phylum_group, ko_opportunity_bin)` ($B = 10{,}000$; Tier 1 $\Delta_{\text{MH}}(E_{\text{fam}}) = -0.0648, p = 0.4907$; Tier 2 $\Delta_{\text{MH}}(E_{\text{fam}}) = +0.0119, p = 0.8678$), (2) 1:1 caliper matching with explicit soft- vs. strict-caliper balance diagnostics, and (3) converged multivariable OLS/HC3 (`is_k` $\hat{\beta} = +0.0083, p = 0.9315$) and logistic regression (`is_k` $\text{OR} = 1.4916, p = 0.6236$).
4. **`C4` Resolved — Cryptographic Contig Sequence SHA-256 + Bakta Log Verification & Complete 3,800-File Checksum Manifest:**  
   In v2, contigs were paired by FASTA order and length only (`209 / 293` bins had duplicate lengths), `14` missing bins were silently dropped, and `INPUTS.tsv` contained `DYNAMIC` placeholders. In v3 (`scripts/02_extract_neighbourhoods.py`), all `709` unique refined bins are audited in `results/bin_contig_verification.tsv`: **`685 / 685` (`100.0%`) available bins pass both exact forward-strand uppercase nucleotide SHA-256 matching and Bakta `.log` `orig-id` verification (`0` hash mismatches, `0` reverse-complement mismatches, despite `521` bins containing duplicate contig lengths)**, while **`24 / 709` bins lacking Bakta/eggNOG files are explicitly recorded (`annotation_files_missing: 24`)**. Every single input file read (`3,800` files total: `11` fixed tables, `340` sample `.faa.gz`, `709` bin `.fa.gz`, and `2,740` Bakta/eggNOG files) is recorded with exact byte size and SHA-256 in `INPUTS.tsv` and `results/executed_inputs_manifest.tsv` (`0` `DYNAMIC` placeholders).
5. **`C5` Resolved — Explicit 90%-Sequence-Cluster Representative Population Definition & Strict Interpretive Bounds:**  
   Upstream provenance (`00_structure_clean/scripts/01_prep.sh`, `04_cluster.sh`, `darkclean_5_mag.sh`) confirms that `04_cluster/clusters.tsv` contains **`12,556` 90%-sequence-cluster representative proteins (`mag_dark_rep_ids.txt.gz`)** across `12,382` unique contigs (`6,138 / 12,382` [`49.57%`] on refined MAGs). All v3 estimands are explicitly conditioned on this verified 90%-representative member population without unverified extrapolation to unbinned or unexpanded `all_cluster90.tsv` instances, and all references to "novel folds", "operons", "expression", or "causal adaptation" are removed.

---

## 2. Bounded Population, Attrition Funnel, & Formal Estimands

### 2.1 Explicit Member Population & Taxonomic Lineage Resolution (`C2`, `C5`)
* **Confidence-Qualified Structural Universe (`Gate 0`):** `2,838` monomer structural families (`mean_plddt >= 0.70`, `length >= 100` aa) from `05_validate/per_rep_structures.tsv` and `05_validate/fam_reps.tsv`:
  * **Arm $K$ (Low-query-TM structural clusters):** $n = 1{,}031$ families with $\text{qTM} < 0.5$ against both PDB and AFDB (`in_K == 1`).
  * **Arm $Q \setminus K$ (Fold-matched sequence-dark structural clusters):** $n = 1{,}807$ families with $\text{qTM} \geq 0.5$ against PDB or AFDB (`in_K == 0`).
* **Member Basis (`04_cluster/clusters.tsv`):** Contains `12,556` 90%-sequence-cluster representative member proteins (`mag_dark_rep_ids.txt.gz`) across `12,382` unique contigs, of which `6,138` contigs (`49.57%`) map to refined MAG bins in `contig_to_bin_refined.tsv.gz`.
* **Species Identity & Verified GTDB-Tk R232 Taxonomy:** Refined bins are mapped to dRep 95% ANI secondary clusters (`species`) via `05_refine/_drep/out/data_tables/Cdb.csv` (`12,433` binned genomes across `5,438` species clusters). Refined bins absent from `Cdb.csv` receive `species = ""` (`binned_unclustered_no_drep_species`) and never count as species. Species-level GTDB-Tk R232 taxonomy is joined via `Wdb.csv` and `05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` (`5,438 / 5,438` [`100.0%`] species clusters resolved).

### 2.2 Complete Attrition Funnel (`2,838` Universe to Evaluable Tiers)

| Funnel Stage / Gate | Total Families | Arm $K$ ($\text{qTM} < 0.5$) | Arm $Q \setminus K$ ($\text{qTM} \geq 0.5$) | Classification of Excluded Families |
|---|---:|---:|---:|---|
| **Gate 0:** Confidence-qualified structural families (`per_rep_structures.tsv`) | 2,838 | 1,031 | 1,807 | — |
| **Gate 1:** $\geq 1$ 90%-rep member on a refined MAG contig (`contig_to_bin_refined.tsv.gz`) | 2,426 | 917 | 1,509 | 412 families with only unbinned 90%-rep contigs |
| **Gate 2:** $\geq 1$ 90%-rep member in a valid dRep 95% ANI species cluster (`Cdb.csv`) | 2,378 | 898 | 1,480 | 48 families only in unclustered refined bins |
| **Gate 3:** $\geq 2$ distinct dRep 95% ANI species clusters (`all_evaluable_families.tsv`) | 1,461 | 570 | 891 | 917 single-species families (`n_species == 1`) |
| **Pilot Cohort:** Stratified sample across `(phylum_group, species_breadth_bin)` (`pilot_cohort.tsv`) | 300 | 150 | 150 |Deterministic seed `20261013` across all 8 strata |
| **Tier 1 (Primary Gate):** $\geq 2$ distinct species with 4 exact-coordinate Bakta-mapped flanks | **72** | **43** | **29** | **228 Ineligible / Missing (`c_fam_primary = -1`)** |
| **Tier 2 (3'-Stop Sensitivity Gate):** $\geq 2$ distinct species with 4 exact or unique 3'-stop+strand+frame+overlap flanks | **133** | **67** | **66** | **167 Ineligible / Missing (`c_fam_stop3_sens = -1`)** |
| **Tier 3 (Partial-Window Sensitivity Gate):** $\geq 2$ distinct species with $\geq 1$ exact-mapped flank | **274** | **138** | **136** | **26 Ineligible / Missing (`c_fam_sensitivity = -1`)** |

### 2.3 Primary Locus Selection & Formal Estimands (`C1`, `C3`)
* **Two-Stage Outcome-Independent Primary Locus Selection:** Across the 300 pilot families, `scripts/01_select_cohort.py` emits all `1,518` candidate 90%-representative member rows (`964` in dRep species clusters, `32` in unclustered refined bins, `522` unbinned). After measuring true nucleotide contig lengths (`contig_len_bp`) and flank mapping status in `scripts/02_extract_neighbourhoods.py`, exactly one deterministic primary locus (`is_primary_locus = 1`; `900` total `(family, species)` pairs) is selected per `(family, species)` using the exact 8-element outcome-independent sort key (`scripts/02_extract_neighbourhoods.py:838-847`, ignoring KO annotations):
  $$\text{key}(m) = \left(-\text{eligible\_locus\_primary}(m), \; -\text{eligible\_locus\_stop3\_sens}(m), \; -\text{n\_exact\_mapped\_flanks}(m), \; -\text{n\_stop3\_or\_exact\_mapped\_flanks}(m), \; -\text{n\_present\_flanks}(m), \; -\text{contig\_len\_bp}(m), \; -\text{total\_genes\_on\_contig}(m), \; \text{member}(m)\right)$$
* **Locus-First Aggregation:** Let $\mathcal{L}_{\text{prim}}(F)$ denote the set of designated primary species loci for family $F$ that satisfy the Tier 1 gate (`eligible_locus_primary == 1`), and $S_F^{\text{prim}} = |\mathcal{L}_{\text{prim}}(F)|$. For $S_F^{\text{prim}} \geq 2$, let $n_{F,k}^{\text{prim}} = \sum_{\ell \in \mathcal{L}_{\text{prim}}(F)} \mathbb{I}(k \in \text{KOs}(\ell))$ be the number of distinct eligible species carrying KEGG KO $k$ in their 4-flank window:
  $$C_{\text{fam}}^{\text{obs}} = \mathbb{I}\left(\max_{k} n_{F,k}^{\text{prim}} \geq 2\right) \in \{0, 1\}, \qquad R_{\text{fam}}^{\text{obs}} = \begin{cases} \dfrac{\max_{k} n_{F,k}^{\text{prim}}}{S_F^{\text{prim}}} \in \left[\dfrac{2}{S_F^{\text{prim}}}, 1\right], & \text{if } \max_{k} n_{F,k}^{\text{prim}} \geq 2, \\[8pt] 0.0, & \text{otherwise.} \end{cases}$$

---

## 3. Cryptographic Contig Identity & Cross-Caller Flank Join Audit (`C1`, `C4`)

### 3.1 Contig Sequence SHA-256 + Bakta Log Verification (`results/bin_contig_verification.tsv`)
Across the `1,518` candidate member rows, `996` reside on `709` unique refined MAG bins (`522` are unbinned contigs):
* **Verified by Exact Forward-Strand Nucleotide SHA-256 + Bakta `.log` `orig-id`:** **`685 / 685` (`100.0%`)** of refined bins with Bakta annotation files passed both cryptographic uppercase nucleotide SHA-256 identity between `03_binning/.../output_bins/{bin}.fa.gz` and `07_annotation/.../bakta/{bin}/{bin}.fna` and explicit `qc: revised sequence: id=contig_N, orig-id=<raw_contig>` verification in `{bin}.log`.
* **Duplicate Contig Lengths Handled Safely:** Even though **`521 / 685` (`76.1%`)** available bins contain multiple contigs with identical base-pair lengths (`0` bins had duplicate sequence SHA-256 hashes), combining sequence SHA-256 with `.log` `orig-id` achieved a **100.0% unique 1-to-1 contig mapping** with **`0` hash mismatches** and **`0` reverse-complement orientation mismatches**.
* **Explicit Missing-Annotation Bins:** Exactly **`24 / 709` (`3.39%`)** refined bins lack Bakta/eggNOG annotation files on IBEX (`annotation_files_missing: 24`, enumerated in `results/pilot_join_metrics.json:10-35` and `results/bin_contig_verification.tsv`). All `35` primary-locus flanks on these `24` bins are classified as `source_file_unavailable` and excluded from callable/complete windows.
* **Complete Input Checksum Pinning:** `INPUTS.tsv` and `results/executed_inputs_manifest.tsv` record exact byte sizes and SHA-256 checksums for all **`3,800` input files** read (`11` fixed upstream tables + `340` sample `.faa.gz` + `709` bin `.fa.gz` + `2,740` Bakta `.fna`/`.log`/`.tsv` and eggNOG `.emapper.annotations` files), with **`0` `DYNAMIC` placeholders**.

### 3.2 Seven-Category Flank Audit & Exact Closure (`results/pilot_join_metrics.json`)
Every flanking position ($\{-2, -1, +1, +2\}$) across all `1,518` candidate rows (`6,072` flanks) is classified into 7 mutually exclusive categories, partitioned transparently across primary species loci (`900` rows, `3,600` flanks), non-primary rows (`618` rows, `2,472` flanks), Tier 1 eligible primary loci (`278` rows, `1,112` flanks), and Tier 1 ineligible primary loci (`622` rows, `2,488` flanks):

| Flank Audit Category | All Candidate Rows (`1,518` rows; `6,072` flanks) | Designated Primary Loci (`900` rows; `3,600` flanks) | Non-Primary Rows (`618` rows; `2,472` flanks) | Tier 1 Eligible Primary Loci (`278` rows; `1,112` flanks) | Tier 1 Ineligible Primary Loci (`622` rows; `2,488` flanks) | Tier 2 Eligible Primary Loci (`458` rows; `1,832` flanks) |
|---|---:|---:|---:|---:|---:|---:|
| `mapped_KO` (exact `(contig, start, stop, strand)` + KO) | 878 (14.46%) | 824 (22.89%) | 54 (2.18%) | **411 (36.96%)** | 413 (16.60%) | **584 (31.88%)** |
| `mapped_noKO` (exact `(contig, start, stop, strand)`, no KO) | 1,548 (25.49%) | 1,456 (40.44%) | 92 (3.72%) | **701 (63.04%)** | 755 (30.35%) | **1,027 (56.06%)** |
| `stop3prime_matched_KO` (unique 3'-stop+strand+frame+overlap + KO) | 145 (2.39%) | 138 (3.83%) | 7 (0.28%) | 0 (0.00%) | 138 (5.55%) | **76 (4.15%)** |
| `stop3prime_matched_noKO` (unique 3'-stop+strand+frame+overlap, no KO) | 287 (4.73%) | 274 (7.61%) | 13 (0.53%) | 0 (0.00%) | 274 (11.01%) | **145 (7.91%)** |
| `coordinate_mismatch_unmapped` (no exact or unique 3'-stop match) | 261 (4.30%) | 235 (6.53%) | 26 (1.05%) | 0 (0.00%) | 235 (9.45%) | 0 (0.00%) |
| `source_file_unavailable` (unbinned contig or missing Bakta bin) | 1,285 (21.16%) | 35 (0.97%) | 1,250 (50.57%) | 0 (0.00%) | 35 (1.41%) | 0 (0.00%) |
| `contig_edge_truncated` (position falls outside $[0, N_{\text{genes}}-1]$) | 1,668 (27.47%) | 638 (17.72%) | 1,030 (41.67%) | 0 (0.00%) | 638 (25.64%) | 0 (0.00%) |
| **Total Flanks Evaluated** | **6,072 (100.0%)** | **3,600 (100.0%)** | **2,472 (100.0%)** | **1,112 (100.0%)** | **2,488 (100.0%)** | **1,832 (100.0%)** |

* **Cross-Caller 3'-Stop Discordance Resolution (`C1`):** Among `647` non-exact binned primary flanks (`138 + 274 + 235`), **`412 / 647` (`63.68%`)** share the exact 3' stop coordinate, strand, reading frame (`start % 3`), and positive nucleotide overlap with a unique Bakta CDS on the same verified contig (5' start-codon differences between meta-Pyrodigal and bin-trained Bakta), rescuing `138` KO-bearing flanks in Tier 2 (`133` evaluable families: `67 K, 66 QnotK`).

---

## 4. Empirical Results (`results/pilot_conservation_summary.json`)

The complete v3 pipeline ran on IBEX under Slurm job `53712511` (`gpu510-32`, 16 CPUs, 3 min 21 s) with post-audit summary verification under Slurm job `53713357`.

### 4.1 Part (a): Within-Arm Recurrence Relative to Same-Contig Eligible-Position Reference ($B = 5{,}000$)
For each Tier 1 eligible primary locus $\ell \in \mathcal{L}_{\text{prim}}(F)$, `scripts/03_analyse_conservation.py` identifies the exact set of gene indices $\mathcal{I}_{\text{elig}}(\ell)$ on the same contig whose 4 flanking genes $(j-2, j-1, j+1, j+2)$ are all exact-coordinate Bakta-mapped (`mapped_KO` or `mapped_noKO`):
* **Eligible Position Counts Across `278` Tier 1 Primary Loci:** `274 / 278` (`98.56%`) reside on contigs with $M_c^{\text{prim}} = |\mathcal{I}_{\text{elig}}(\ell)| \geq 2$ eligible 4-flank windows (`249 / 278` [`89.57%`] have $M_c^{\text{prim}} \geq 5$; $\min = 1$, median $= 20.0$, mean $= 44.71$, $\max = 472$). Only `4 / 278` (`1.44%`) loci have $M_c^{\text{prim}} = 1$.
* **Within-Arm Enrichment Results ($B = 5{,}000$ Monte Carlo draws):**

| Arm | Evaluable Families ($N$) | Observed Conserved Rate $\bar{C}_{\text{arm}}^{\text{obs}}$ | Same-Contig Positional Expected Rate $\bar{\mu}_{\text{arm}}^{\text{pos}}$ (Null SD) | Within-Arm Excess $\bar{E}_{\text{arm}} = \bar{C}^{\text{obs}} - \bar{\mu}^{\text{pos}}$ | Within-Arm $p_{\text{pos}}$ ($C_{\text{fam}}$) | Observed Mean $\bar{R}_{\text{arm}}^{\text{obs}}$ | Positional Expected $\bar{\rho}_{\text{arm}}^{\text{pos}}$ | Within-Arm Excess $\bar{E}_{\text{arm}}^{(R)}$ ($p_{\text{pos}}$) | Multi-Position Subset ($M_c^{\text{prim}} \geq 2$ in $\geq 2$ sp) Excess |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Low-qTM ($K$)** | 43 | **48.84%** (21/43) | 12.86% (3.93%) | **+35.97 pp** | **0.0002** | 0.4729 | 0.1056 | +0.3672 (**0.0002**) | +35.97 pp ($n=43$) |
| **Fold-Matched ($Q \setminus K$)** | 29 | **31.03%** (9/29) | 4.96% (3.63%) | **+26.08 pp** | **0.0002** | 0.3103 | 0.0384 | +0.2720 (**0.0002**) | +29.09 pp ($n=26$) |

* **Finding (Part a):** Both low-query-TM ($K$) and fold-matched ($Q \setminus K$) structural families exhibit highly significant within-arm flanking KO recurrence relative to random eligible 4-flank windows on the exact same host contigs ($+35.97$ pp in $K$ and $+26.08$ pp in $Q \setminus K$, both $p = 0.0002$).

### 4.2 Unadjusted Comparison Across All Three Eligibility Tiers & Confounder Profile
Before covariate adjustment, raw recurrence rates are higher in Arm $K$ than in Arm $Q \setminus K$ across all three eligibility tiers, but without reaching statistical significance:

| Eligibility Tier | Evaluable Families ($K$ / $Q \setminus K$) | $K$ Conserved ($C_{\text{fam}}=1$) | $Q \setminus K$ Conserved ($C_{\text{fam}}=1$) | Raw Risk Difference $\Delta_{\text{raw}}(C_{\text{fam}})$ | Fisher Exact $p$ | Mean $R_{\text{fam}}$ ($K$ vs. $Q \setminus K$) | Mann-Whitney $p$ ($R_{\text{fam}}$) | Mean Excess $E_{\text{fam}}$ ($K$ vs. $Q \setminus K$) | Mann-Whitney $p$ ($E_{\text{fam}}$) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Tier 1 (Primary: 4 exact flanks)** | 72 (43 / 29) | 21 / 43 (48.84%) | 9 / 29 (31.03%) | +17.80 pp | 0.1515 | 0.4729 vs. 0.3103 | 0.1763 | +0.3597 vs. +0.2608 ($\Delta = +0.0989$) | 0.2779 |
| **Tier 2 (Sensitivity: 4 exact/3'-stop flanks)** | 133 (67 / 66) | 31 / 67 (46.27%) | 24 / 66 (36.36%) | +9.91 pp | 0.2920 | 0.4084 vs. 0.3212 | 0.2467 | +0.3391 vs. +0.2976 ($\Delta = +0.0416$) | 0.6852 |
| **Tier 3 (Sensitivity: $\geq 1$ exact flank)** | 274 (138 / 136) | 45 / 138 (32.61%) | 33 / 136 (24.26%) | +8.34 pp | 0.1417 | 0.2526 vs. 0.1982 | 0.1622 | — | — |

* **Why Unadjusted $K$ Appears Higher in Tier 1:** Inspection of the Tier 1 covariate profile (`results/pilot_conservation_summary.json:596-626`) reveals that while species breadth is matched (`2.349` in $K$ vs. `2.379` in $Q \setminus K$, $\text{SMD} = -0.0385$), Tier 1 $K$ families reside in windows with substantially higher **flank KO annotation opportunity** (`mean = 0.4219` in $K$ vs. `0.3273` in $Q \setminus K$, $\text{SMD} = +0.4431$) and longer focal protein length (`309.4` aa vs. `206.6` aa, $\text{SMD} = +0.9298$). Because higher KO density directly increases both random positional recurrence ($\bar{\mu}_K^{\text{pos}} = 12.86\%$ vs. $\bar{\mu}_{Q\setminus K}^{\text{pos}} = 4.96\%$) and observed recurrence detection, valid between-arm inference requires conditioning on annotation opportunity, species breadth, and coarse host lineage.

### 4.3 Part (b): Between-Arm Covariate-Adjusted Contrast ($K$ vs. $Q \setminus K$)

#### 1. Primary Estimator — Covariate-Stratified Label Permutation ($B = 10{,}000$, `results/pilot_common_support_strata.tsv`)
Evaluable families are stratified jointly by `species_breadth_bin` (`sp2` vs. `sp3plus`), verified R232 `phylum_group` (`Actinomycetota`, `Thermoplasmatota`, `other_bacteria_archaea`, `multi`), and `ko_opportunity_bin` (`<= median` vs. `> median` `ko_annotation_opportunity`). Under the identifying assumption of conditional exchangeability of families within joint strata (noting that 4-bin `phylum_group` stratification provides coarse taxonomic conditioning rather than full phylogenetic covariance control), arm labels are permuted strictly within common-support strata ($n_{s,K} \geq 1$ and $n_{s,Q\setminus K} \geq 1$):

| Stratified Permutation Metric ($B = 10{,}000$) | Tier 1 Primary Gate (`72` families) | Tier 2 3'-Stop Sensitivity Gate (`133` families) |
|---|---|---|
| **Common-Support Strata Coverage** | **13 / 15 strata** (`42/43 K` [97.67%], `28/29 QnotK` [96.55%]) | **16 / 16 strata** (`67/67 K` [100.0%], `66/66 QnotK` [100.0%]) |
| **MH-Weighted Excess Contrast $\Delta_{\text{MH}}(E_{\text{fam}})$** | **$-0.0648$** (Null SD $= 0.0937$, **$p = 0.4907$**) | **$+0.0119$** (Null SD $= 0.0712$, **$p = 0.8678$**) |
| **Pooled Common-Support Excess Contrast $\Delta_{\text{pool}}(E_{\text{fam}})$** | $+0.0821$ (Null SD $= 0.0788$, $p = 0.7629$) | $+0.0416$ (Null SD $= 0.0620$, $p = 0.5616$) |
| **MH-Weighted Raw Contrast $\Delta_{\text{MH}}(C_{\text{fam}})$** | **$-0.0034$** (Null SD $= 0.0965$, **$p = 1.0000$**) | **$+0.0717$** (Null SD $= 0.0768$, **$p = 0.3768$**) |
| **MH-Weighted Continuous Excess $\Delta_{\text{MH}}(E_{\text{fam}}^{(R)})$** | **$-0.0687$** ($p = 0.4384$) | **$+0.0152$** ($p = 0.8115$) |

* **Key Result:** Once stratified by species breadth, coarse R232 host lineage, and flank KO annotation opportunity, the apparent unadjusted Tier 1 difference ($+17.80$ pp) is not supported ($\Delta_{\text{MH}}(E_{\text{fam}}) = -0.0648, p = 0.4907$; $\Delta_{\text{MH}}(C_{\text{fam}}) = -0.0034, p = 1.0000$).

#### 2. Diagnostic Estimator — 1:1 Covariate Caliper Matching (`results/pilot_matched_pairs.tsv`)
Because $K$ and $Q \setminus K$ differ structurally in sequence length (`309.4` aa vs. `206.6` aa, pre-match $\text{SMD} = +0.9298$), 1:1 matching without replacement in $n = 72$ faces a direct trade-off between pair count and multidimensional balance across all 5 covariates. We report both soft-caliper and strict exact-phylum/tight-length caliper matching transparently:
* **Soft-Caliper Matching (`25` pairs; exact `species_breadth_bin`, $|\Delta z_{\text{opp}}| \leq 0.75$, $|\Delta z_{\text{len}}| \leq 1.25$, soft phylum penalty):** Balances the two primary drivers of recurrence — `n_species_complete_window` ($\text{SMD} = -0.0645$) and `ko_annotation_opportunity` ($\text{SMD} = +0.0059$, reduced from $+0.4431$) — and yields **null between-arm differences**: $\Delta_{\text{match}}(E_{\text{fam}}) = +0.0195$ (sign-flip $p = 0.8966$) and $\Delta_{\text{match}}(C_{\text{fam}}) = +0.0400$ ($p = 1.0000$). However, this 25-pair soft match is **not balanced on sequence length ($\text{SMD} = +0.6209$) or phylum (`44.0%` exact match)** and must not be described as globally balanced across all covariates.
* **Strict Exact-Phylum + Tight-Length Caliper Sensitivity (`9` pairs; exact `phylum_group` [`100.0%`], exact `species_breadth_bin`, $|\Delta z_{\text{opp}}| \leq 0.60$, $|\Delta z_{\text{len}}| \leq 0.50$):** Enforces exact coarse `phylum_group` matching (`100%` same phylum group; coarse taxonomic conditioning rather than full phylogenetic covariance control) and reduces post-match SMDs to `n_species_complete_window = 0.0000`, `ko_annotation_opportunity = +0.0723`, and `length = +0.2879`, while similarly showing **no significant difference between arms**: $\Delta_{\text{strict}}(E_{\text{fam}}) = +0.0563$ (sign-flip $p = 0.5864$) and $\Delta_{\text{strict}}(C_{\text{fam}}) = +0.1111$ ($p = 1.0000$).

#### 3. Parametric & Penalized Multivariable Regression Models ($n = 72$, `rank = 9 / 9`, All `converged = True`)
Simultaneous adjustment for all covariates (`is_k`, `n_species_complete_window`, `ko_annotation_opportunity`, `length_100aa`, `mean_plddt_scaled`, and R232 `phylum_group` indicators with reference `Thermoplasmatota`) was verified across four model specifications:

| Covariate / Parameter | OLS / HC3 on Excess $E_{\text{fam}}$ ($\hat{\beta} \pm \text{SE}_{\text{HC3}}, p$) | OLS / HC3 on Raw $C_{\text{fam}}$ ($\hat{\beta} \pm \text{SE}_{\text{HC3}}, p$) | MLE Logistic on Raw $C_{\text{fam}}$ ($\hat{\beta} \pm \text{SE}$, OR, $p$) | L2-Penalized Ridge Logistic ($\lambda = 0.5$) ($\hat{\beta} \pm \text{SE}$, OR, $p$) |
|---|---|---|---|---|
| **Arm Indicator (`is_k`: $K$ vs. $Q \setminus K$)** | **$+0.0083 \pm 0.0961$ ($p = 0.9315$)** | **$+0.0601 \pm 0.1084$ ($p = 0.5791$)** | **$+0.3999 \pm 0.8147$ ($\text{OR} = 1.49, p = 0.6236$)** | **$+0.5492 \pm 0.5771$ ($\text{OR} = 1.73, p = 0.3413$)** |
| **Species Breadth ($S_F^{\text{prim}}$)** | $-0.0659 \pm 0.0636$ ($p = 0.3001$) | $+0.0055 \pm 0.0622$ ($p = 0.9289$) | $-0.0866 \pm 0.4283$ ($\text{OR} = 0.92, p = 0.8398$) | $+0.0924 \pm 0.3406$ ($\text{OR} = 1.10, p = 0.7862$) |
| **KO Annotation Opportunity ($A_{\text{flank}}^{\text{prim}}$)** | **$+0.9513 \pm 0.2241$ ($p = 2.2 \times 10^{-5}$)** | **$+1.2729 \pm 0.1883$ ($p < 10^{-6}$)** | **$+11.1593 \pm 3.0922$ ($p = 0.0003$)** | **$+3.5413 \pm 1.0047$ ($p = 0.0004$)** |
| **Sequence Length (per 100 aa)** | $-0.0237 \pm 0.0374$ ($p = 0.5275$) | $-0.0248 \pm 0.0465$ ($p = 0.5938$) | $+0.0148 \pm 0.2984$ ($\text{OR} = 1.01, p = 0.9604$) | $-0.0802 \pm 0.2524$ ($\text{OR} = 0.92, p = 0.7506$) |
| **Mean pLDDT (scaled per 0.10)** | $-0.0360 \pm 0.1065$ ($p = 0.7351$) | $+0.0023 \pm 0.0975$ ($p = 0.9815$) | $+0.4857 \pm 0.7547$ ($\text{OR} = 1.63, p = 0.5199$) | $+0.0565 \pm 0.5257$ ($\text{OR} = 1.06, p = 0.9144$) |
| **Lineage: `Actinomycetota` vs. `Thermoplasmatota`** | $+0.0688 \pm 0.1595$ ($p = 0.6664$) | $+0.0665 \pm 0.1742$ ($p = 0.7025$) | $+0.0734 \pm 1.0186$ ($\text{OR} = 1.08, p = 0.9425$) | $+0.6164 \pm 0.6778$ ($\text{OR} = 1.85, p = 0.3631$) |
| **Lineage: `multi` vs. `Thermoplasmatota`** | $-0.2248 \pm 0.1448$ ($p = 0.1204$) | $-0.3229 \pm 0.1608$ ($p = 0.0446$) | $-3.4777 \pm 1.4067$ ($\text{OR} = 0.03, p = 0.0134$) | $-1.1622 \pm 0.7212$ ($\text{OR} = 0.31, p = 0.1071$) |
| **Lineage: `other_bacteria_archaea` vs. `Thermoplasmatota`** | $-0.0065 \pm 0.1551$ ($p = 0.9664$) | $+0.0082 \pm 0.1806$ ($p = 0.9639$) | $+0.1117 \pm 1.0201$ ($\text{OR} = 1.12, p = 0.9128$) | $+0.4363 \pm 0.6855$ ($\text{OR} = 1.55, p = 0.5245$) |
| **Model Diagnostics** | `converged=True`, $R^2 = 0.3141$ | `converged=True`, $R^2 = 0.4371$ | `converged=True`, pseudo-$R^2 = 0.4387$ | `converged=True`, all $\text{SE} \leq 1.01$ |

### 4.4 Cross-Family Contig, Bin, & Species Dependence Audit & Exchangeability Assumptions (`results/pilot_conservation_summary.json:647-667`)
Conditional exchangeability within strata/matches and cross-family independence are **identifying assumptions** of the between-arm tests rather than empirically verified facts, and 4-bin R232 `phylum_group` conditioning or exact phylum-group matching (`9` strict pairs) represents **coarse taxonomic conditioning rather than full phylogenetic covariance control**. To audit physical locus overlap and shared host taxa across the `72` Tier 1 evaluable families (`170` eligible primary loci):
* **Zero Shared Contigs or Overlapping Windows:** All **`170 / 170` (`100.0%`)** eligible primary loci across the `72` Tier 1 evaluable families reside on **distinct `(bin, contig)` pairs (`170` unique contigs; `0` shared contigs across families; `0` overlapping $\pm 2$-gene windows)**.
* **Shared Host Lineage Audit & Fixed-Outcome Species-Reweighting Sensitivity:** The `170` loci span `151` unique refined bins (`17` bins shared by $\geq 2$ families) and `144` unique dRep 95% ANI species clusters (`21` species clusters shared by $\geq 2$ families, including `8` shared across both arms; at broader taxonomic ranks, `13 / 57` named GTDB genera and `19 / 56` named GTDB families appear in $\geq 2$ Tier 1 structural families). In `scripts/03_analyse_conservation.py:859-883`, resampling the `144` dRep species clusters with replacement ($B = 2{,}000$) weights each family's **fixed pre-computed** $E_{\text{fam}}$ and $C_{\text{fam}}$ by its mean resampled species multiplicity to evaluate the **unstratified** arm difference (`+0.0989` baseline), yielding $\text{SE}_{\text{boot}}(\Delta E_{\text{fam}}) = 0.0730$ and $\text{SE}_{\text{boot}}(\Delta C_{\text{fam}}) = 0.0852$, while restricting to the `40` families with zero shared species (`26 K, 14 QnotK`) yields an unstratified difference of $\Delta(E_{\text{fam}}) = +0.1091$ (Mann-Whitney $p = 0.5452$, Welch $t$ $p = 0.4036$). Because this species-reweighting check holds family-level outcomes fixed and does not recompute multi-species KO recurrence ($\geq 2$ distinct resampled species) or the primary Mantel-Haenszel stratified estimator ($\Delta_{\text{MH}}(E_{\text{fam}}) = -0.0648$), it is reported strictly as a descriptive sensitivity check and cannot prove the absence of masking or phylogenetic confounding.

---

## 5. Descriptive Neighbourhood Context (`O1`) & Recurrent Flanking KOs

### 5.1 Strand Co-Orientation & Intergenic Distances (`O1`)
Across the `278` Tier 1 eligible complete-window primary loci (`151` in Arm $K$, `127` in Arm $Q \setminus K$):
* **Same-Strand Fraction (`same_strand_fraction`):** Family-level mean co-orientation with the focal CDS across the 4 flanks is **`0.6441`** (median `0.6250`) in Arm $K$ (`n = 43`) vs. **`0.6174`** (median `0.6667`) in Arm $Q \setminus K$ (`n = 29`); locus-level mean is `0.6325` in $K$ vs. `0.6024` in $Q \setminus K$.
* **Intergenic Distances (`_dist`):** Median distance from the focal CDS to immediately adjacent flanks (`-1, +1`) is **`70.5` bp** ($\text{IQR}: 14.2\text{--}171.8$ bp) in Arm $K$ and **`65.5` bp** ($\text{IQR}: 10.5\text{--}137.5$ bp) in Arm $Q \setminus K$ (across all 4 flanks `m2, m1, p1, p2`: median `400.5` bp in $K$ vs. `317.5` bp in $Q \setminus K$). These metrics are reported strictly as descriptive physical neighbourhood statistics without inferring co-transcription.

### 5.2 Complete List of 30 Tier 1 Recurrent Flanking KO Families & Top-KO Keyword Flag Audit
Across the `30` Tier 1 conserved families (`21 K`, `9 QnotK`), **`0 / 30` (`0.0%`)** of the selected `best_recurrent_ko` annotations match the 11-keyword substring screen (`MOBILOME_TERMS`: transposase, integrase, recombinase, insertion element, phage, prophage, terminase, capsid, tail, portal, plasmid; `best_ko_mobilome_flag == 0` in `scripts/03_analyse_conservation.py:99-116,189-197`). This `0 / 30` flag applies strictly to keyword screening of the selected top recurrent KO in each family: it does **not** establish the absence of mobile or defense elements across the windows (across all `296` Tier 1 flanks in these `30` families, `141 / 296` [`47.6%`] have no KEGG KO [`mapped_noKO`], `172 / 296` [`58.1%`] have Bakta product `"hypothetical protein"`, `K00590` in `F31Dr2::F31Dr2_428982_6` encodes a `DNA methylase`, and an adjacent `mapped_noKO` flank in `F38Sr1::F38Sr1_386592_2` carries eggNOG description `"COG2931 RTX toxins and related Ca2 -binding proteins"`), nor does it demonstrate that all 30 loci are chromosomal or core housekeeping/metabolic genes (`8 / 30` top recurrent flanks have Bakta product `"hypothetical protein"` and `6 / 30` have eggNOG COG category `S` [function unknown], alongside `K06888` `DUF255 domain-containing protein` and `K02221` `YggT family protein`), or establish any functional mechanism:

| Family ID | Arm | Verified R232 Phylum | Best Recurrent KO | Verified Bakta Product of Recurrent Flank | Species ($n / S_F^{\text{prim}}$) | $R_{\text{fam}}^{\text{obs}}$ | Excess $E_{\text{fam}}$ | Mobilome Flag |
|---|---|---|---|---|---:|---:|---:|---:|
| `1Dr2::1Dr2_476456_56` | `K` | `Actinomycetota` | `K06888` | DUF255 domain-containing protein | 7 / 7 | 1.0000 | +0.5426 | 0 |
| `33PRr2::33PRr2_1315111_5` | `K` | `Actinomycetota` | `K03545` | hypothetical protein | 2 / 2 | 1.0000 | +0.2250 | 0 |
| `36PRr1::36PRr1_462970_6` | `QnotK` | `Actinomycetota` | `K18455` | Mycothiol conjugate amidase Mca | 2 / 2 | 1.0000 | +0.7806 | 0 |
| `3Dr2::3Dr2_768169_20` | `K` | `Thermoplasmatota` | `K00260` | hypothetical protein | 2 / 2 | 1.0000 | +0.9170 | 0 |
| `44Dr2::44Dr2_3882404_4` | `K` | `Thermoplasmatota` | `K00261` | Glutamate dehydrogenase | 2 / 2 | 1.0000 | +0.7216 | 0 |
| `54PRr2::54PRr2_192771_23` | `K` | `Pseudomonadota` | `K01420` | hypothetical protein | 2 / 3 | 0.6667 | +0.9128 | 0 |
| `F23PRr1::F23PRr1_142702_24` | `K` | `Actinomycetota` | `K13280` | hypothetical protein | 2 / 2 | 1.0000 | +0.9466 | 0 |
| `F25Dr2::F25Dr2_61875_4` | `K` | `Thermoplasmatota` | `K06883` | hypothetical protein | 4 / 4 | 1.0000 | +0.6236 | 0 |
| `F2Dr1::F2Dr1_2447321_2` | `K` | `Acidobacteriota` | `K02045` | ABC transporter domain-containing protein | 2 / 2 | 1.0000 | +0.6724 | 0 |
| `F31Dr2::F31Dr2_428982_6` | `QnotK` | `Thermoproteota` | `K00590` | hypothetical protein | 2 / 2 | 1.0000 | +0.9036 | 0 |
| `F34PRr2::F34PRr2_1575772_4` | `K` | `multi` | `K02888` | Large ribosomal subunit protein bL21 | 2 / 2 | 1.0000 | +0.9362 | 0 |
| `F38Sr1::F38Sr1_386592_2` | `K` | `Thermoplasmatota` | `K03924` | ATPase | 2 / 2 | 1.0000 | +0.8266 | 0 |
| `F41Dr3::F41Dr3_524800_3` | `K` | `Actinomycetota` | `K01870` | Isoleucine--tRNA ligase | 2 / 2 | 1.0000 | +0.9262 | 0 |
| `F41Dr3::F41Dr3_561308_4` | `QnotK` | `Gemmatimonadota` | `K03924` | ATPase | 3 / 3 | 1.0000 | +0.8910 | 0 |
| `F41Dr3::F41Dr3_7398_3` | `QnotK` | `multi` | `K02108` | ATP synthase subunit a | 2 / 2 | 1.0000 | +0.9650 | 0 |
| `F6Dr2::F6Dr2_614823_13` | `QnotK` | `Myxococcota_A` | `K00360` | Bacterioferritin-associated ferredoxin | 2 / 2 | 1.0000 | +0.9622 | 0 |
| `S22PRr1::S22PRr1_266671_1` | `QnotK` | `Actinomycetota` | `K07161` | hypothetical protein | 5 / 5 | 1.0000 | +0.5312 | 0 |
| `S2PRr1::S2PRr1_94455_15` | `K` | `Acidobacteriota` | `K06867` | hypothetical protein | 2 / 3 | 0.6667 | +0.4684 | 0 |
| `V12Dr3::V12Dr3_456902_4` | `QnotK` | `Actinomycetota` | `K02902` | Large ribosomal subunit protein bL28 | 2 / 2 | 1.0000 | +0.8824 | 0 |
| `V12Sr1::V12Sr1_1454330_2` | `K` | `Actinomycetota` | `K02221` | YggT family protein | 2 / 2 | 1.0000 | +0.5838 | 0 |
| `V15PRr2::V15PRr2_848361_1` | `QnotK` | `Myxococcota_A` | `K00075` | UDP-N-acetylenolpyruvoylglucosamine reductase | 2 / 2 | 1.0000 | +0.8392 | 0 |
| `V30PRr1::V30PRr1_1219613_7` | `K` | `Actinomycetota` | `K07566` | Threonylcarbamoyl-AMP synthase | 2 / 2 | 1.0000 | +0.9676 | 0 |
| `V30PRr1::V30PRr1_62945_5` | `K` | `Actinomycetota` | `K00760` | tRNA(Ile)-lysidine synthase | 2 / 2 | 1.0000 | +0.9126 | 0 |
| `V37PRr2::V37PRr2_488468_1` | `K` | `Cyanobacteriota` | `K02033` | Glutathione transport system permease protein GsiC | 3 / 3 | 1.0000 | +0.0096 | 0 |
| `V39Sr2::V39Sr2_267206_2` | `K` | `multi` | `K06902` | Major facilitator superfamily (MFS) profile domain-containing protein | 2 / 2 | 1.0000 | +0.8148 | 0 |
| `V58PRr1::V58PRr1_187836_44` | `K` | `multi` | `K01802` | Peptidyl-prolyl cis-trans isomerase | 2 / 2 | 1.0000 | +0.9696 | 0 |
| `V58PRr1::V58PRr1_97362_17` | `K` | `Cyanobacteriota` | `K06929` | CoA-binding domain-containing protein | 2 / 2 | 1.0000 | +0.8914 | 0 |
| `V5Dr1::V5Dr1_585320_11` | `K` | `Actinomycetota` | `K07177` | endopeptidase La | 2 / 2 | 1.0000 | +0.9796 | 0 |
| `V60PRr2::V60PRr2_435968_3` | `K` | `multi` | `K01647` | Citrate synthase | 2 / 2 | 1.0000 | +0.9446 | 0 |
| `V8PRr1::V8PRr1_1011131_1` | `QnotK` | `Actinomycetota` | `K00575` | protein-glutamate O-methyltransferase | 2 / 2 | 1.0000 | +0.9816 | 0 |

---

## 6. Independent Closure Assertions & Full-Cohort Execution Recipe

### 6.1 Programmatic Closure Verification (`results/pilot_conservation_summary.json:668-680`)
All 11 programmatic closure assertions in `scripts/03_analyse_conservation.py` passed (`True`):
1. `c1_all_tier1_eligible_loci_have_4_exact_mapped_flanks: true` (`1,112 / 1,112` flanks in `{"mapped_KO", "mapped_noKO"}`).
2. `c1_all_ineligible_families_coded_minus_one: true` (`228` Tier 1, `167` Tier 2, `26` Tier 3 families coded `-1`).
3. `c1_all_r_fam_within_0_1_bounds: true` (`0.0 <= R_fam <= 1.0` across all 300 families in all 3 tiers).
4. `c1_numerator_denominator_species_closure: true` (`best_ko_species_count <= n_species_complete_window` for all 300 families).
5. `c2_zero_unassigned_r232_phylum_in_cohort: true` (`0` unassigned phyla across all `1,461` Gate 3 and `300` pilot families).
6. `c2_all_binned_loci_have_positive_contig_len_bp: true` (`996 / 996` binned rows have `contig_len_bp > 0`).
7. `c3_all_positional_null_draws_satisfy_primary_gate: true` (`278 / 278` Tier 1 loci have $M_c^{\text{prim}} \geq 1$).
8. `c3_between_arm_stratified_common_support_verified: true` (`13` strata covering `97.67%` $K$ and `96.55%` $Q \setminus K$ in Tier 1; `16 / 16` strata covering `100%` in Tier 2).
9. `c3_ols_hc3_and_logit_models_converged_full_rank_finite_se: true` (`rank = 9/9`, `converged = True`, all $\text{SE} < 3.1$).
10. `c4_flank_category_totals_close_exactly: true` (`3,600 + 2,472 = 6,072`; `1,112 + 2,488 = 3,600`).
11. `c5_cohort_attrition_totals_close_exactly: true` (`72 + 228 = 300`; `133 + 167 = 300`; `274 + 26 = 300`).

### 6.2 Runnable Full-Cohort Recipe (`1,461` Gate 3 Families) & Interpretive Boundaries
* **Full-Cohort Command (`<= 16` CPUs, `64` GB RAM, `<= 30` min on IBEX `batch`/`debug`):**
  ```bash
  cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context
  sbatch scripts/run_pilot.sbatch full 1461 results_full false
  ```
  Extrapolating the pilot eligibility fractions (`24.0%` Tier 1 [`28.7%` $K$, `19.3%` $Q \setminus K$], `44.3%` Tier 2) across all `1,461` Gate 3 families (`570 K, 891 QnotK`) provides a pilot-based sample-size planning estimate of $\approx 336\text{--}351$ Tier 1 and $\approx 648$ Tier 2 evaluable families (reported strictly as a finite-catalogue planning projection rather than a guaranteed statistical power claim).
* **Scientific Conclusion & Scope Bounds:** Both low-query-TM ($K$) and fold-matched ($Q \setminus K$) structural dark families exhibit strong, non-random cross-species flanking KO recurrence relative to random eligible windows on the same host contigs ($p = 0.0002$ in both arms), while `0 / 30` selected top recurrent KOs match the 11-term phage/transposase/integrase keyword screen (`best_ko_mobilome_flag == 0`; this keyword flag applies only to the selected recurrent KO and does not establish the absence of mobile elements, chromosomal location, or core housekeeping/metabolic status, as several recurrent flanks are annotated as `"hypothetical protein"` or `DUF255`/`YggT` domains of unknown function). Under the identifying assumptions of conditional exchangeability and cross-family independence — conditioning on species breadth, coarse 4-group R232 host phylum classification, and flank KO annotation opportunity — recurrence excess does not differ significantly between $K$ and $Q \setminus K$ ($\Delta_{\text{MH}}(E_{\text{fam}}) = -0.0648, p = 0.4907$). These results support bounded computational gene-context hypotheses for specific high-recurrence families (e.g., `1Dr2::1Dr2_476456_56` adjacent to `K06888` [`DUF255`] across 7/7 species; `F25Dr2::F25Dr2_61875_4` adjacent to `K06883` across 4/4 species) without claiming functional mechanisms, novel protein folds, operon co-transcription, gene expression, or adaptive evolution.
