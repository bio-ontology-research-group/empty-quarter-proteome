# DESIGN: Flanking Gene Neighbourhood Recurrence of Low-Query-TM vs Fold-Matched Sequence-Dark Structural Families (v3 Protocol)

**Hypothesis 2:** Cross-species flanking KEGG KO recurrence of low-query-TM ($K$, $\text{qTM} < 0.5$) versus fold-matched sequence-dark ($Q \setminus K$, $\text{qTM} \geq 0.5$) structural families.  
**Date:** 2026-10-13  
**Status:** Version 3 (Pre-computation protocol resolving `round6/h2_review1` `C1`–`C5` and `round6/ROOT_AUDIT_NOTES.md`)  
**Local Workspace:** `reproducibility/darkmatter/round6_hypotheses/h2_context/`  
**Remote Workbench:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/`  
**Superseded Versions:** Version 1 (`v1_superseded/`, Slurm `53701050`) and Version 2 (`v2_superseded/` and `reviews/darkmatter-2026-10-13/round6/snapshots/h2_v2/`, Slurm `53703301`).

---

## 1. Audit Reconciliation & Protocol Changes from v2 (`C1`–`C5`)

| Finding | Defect in v2 | Mandatory Protocol & Code Specification in v3 |
|---|---|---|
| **C1: Locus Eligibility & Callability** | `ko_by_sp_primary` aggregated KOs from non-truncated flanks of incomplete-window loci ($R_{\text{fam}} > 1.0$ in 4 families, 2 false-positive `QnotK` calls); `source_file_unavailable` and `coordinate_mismatch` flanks set `is_complete = True` and `has_any_callable = True`. | Determine locus-level eligibility **before** KO aggregation. A primary-eligible locus (`eligible_locus_primary == 1`) requires all 4 flanks ($m2, m1, p1, p2$) to be exact-coordinate Bakta-mapped (`mapped_KO` or `mapped_noKO`). Unmapped, unavailable, or truncated flanks are strictly missing, never observed `noKO`. Primary numerator, denominator ($n_{\text{sp,prim}}$), and $A_{\text{flank}}$ use the exact same eligible loci/species, enforcing $R_{\text{fam}} \in [0, 1]$. Audit exact 3'-stop+strand cross-caller matches (`stop3prime_matched_KO`, `stop3prime_matched_noKO`) as a separate sensitivity tier (`Tier 2`), never substituting silently for exact coordinates (`Tier 1`). |
| **C2: Multiplicity Contig Lengths, R232 Taxonomy, & Model Convergence** | `contig_len` parsed column 2 (`"refined"`) of `contig_to_bin_refined.tsv.gz`, yielding `0` on 534/534 rows; taxonomy read unrefined `02_taxonomy/gtdbtk_out/` IDs, leaving 99/100 pilot families as `phylum = "unassigned"`, inducing quasi-complete separation and `"converged": false`. | Measure true nucleotide `contig_len_bp` from assembly/bin FASTA and evaluate outcome-independent locus completeness **before** selecting the single deterministic primary locus per `(family, species)`. Join verified GTDB-Tk R232 taxonomy through `Cdb.csv` (`secondary_cluster`) $\to$ `Wdb.csv` (`winner_genome`) $\to$ `05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` (100% coverage across all 12,433 dRep-clustered MAGs and 1,461 multi-species families) plus direct bin cross-checks in `05_refine/_gtdbtk_refined/`. Enforce full-rank, convergence, and finite-SE checks on all regression fits. |
| **C3: Same-Contig Positional Reference & Between-Arm Excess Inference** | Positional null sampled incomplete loci and edge positions (33.5% truncated draws) and compared observed $\Delta^{\text{obs}}$ against a low-prevalence positional-null difference distribution (`SD = 0.0681` vs. `0.1336`), creating an artifactual $p = 0.0010$. | Separate **(a) within-arm recurrence relative to a same-contig positional reference** (sampling pseudo-focal positions strictly from indices $j \in \mathcal{I}_c$ on the same contig that satisfy the **identical** 4-flank exact-mapped gate as observed loci, recording $M_c = |\mathcal{I}_c|$ per locus) from **(b) between-arm comparison** on family-level background-adjusted excess $E_{\text{fam}} = C_{\text{fam}}^{\text{obs}} - \mu_{\text{fam}}^{\text{pos}}$. Test between-arm differences using covariate-stratified label permutation across joint `(species_breadth_bin, phylum_group, ko_opportunity_bin)` strata, 1:1 Mahalanobis caliper matching within species-breadth/lineage blocks, and converged regression models. Expand deterministic pilot cohort to **300 families (150 $K$, 150 $Q \setminus K$)** (and execute/provide the full 1,461-family cohort recipe) to ensure common support across strata. |
| **C4: Proven Contig Identity & Complete Input Checksum Manifest** | Raw-to-Bakta contig mapping used FASTA order/length despite 209/293 bins having duplicate contig lengths; 14/307 bins lacked Bakta files on IBEX; `INPUTS.tsv` had `DYNAMIC` placeholders and omitted `.fna`/`.log`. | Prove 1:1 raw-bin-to-Bakta contig identity via exact uppercase nucleotide SHA-256 matching between `output_bins/{parent}.fa.gz` and `bakta/{parent}/{parent}.fna` **plus** `id=contig_N, orig-id=...` verification in `bakta/{parent}/{parent}.log`, failing closed on duplicate sequence hashes, reverse-complement mismatches, or log discrepancies. Record missing-annotation bins explicitly. Enumerate and SHA-256-hash every read file in `INPUTS.tsv` and `results/executed_inputs_manifest.tsv` with zero `DYNAMIC` placeholders. |
| **C5: Exact Member Population & Scientific Bounds** | `04_cluster/clusters.tsv` contains 90%-sequence-cluster representatives (`mag_dark_rep_ids.txt.gz`) rather than all gene instances in `all_cluster90.tsv`, which was undocumented; `DOSSIER.md` retained `"candidate novel folds"` and `"operons"` claims. | Explicitly define and bound the target population on the verified 90%-sequence-cluster representative members in `04_cluster/clusters.tsv` mapped to dRep species-clustered refined MAGs, documenting why unbinned 90%-cluster reps arise from `darkclean_5_mag.sh`. Eliminate all novel-fold, operon, co-transcription, expression, and adaptive claims. |

---

## 2. Formal Target Population, Membership Basis, & Attrition Funnel

### 2.1 Upstream Provenance & Exact Member Population (`C5`)
* **Structural Universe ($Q$):** $N = 2{,}838$ confidence-qualified monomer ESMFold models (mean $\text{pLDDT} \geq 0.70$, length $\geq 100$ aa) from `04_structure/00_structure_clean/05_validate/per_rep_structures.tsv`:
  * **Low-query-TM arm ($K$):** $n = 1{,}031$ structural families with query-normalised Foldseek TM-score $\text{qTM} < 0.5$ against both PDB and AFDB (`conf_qualified_dual_novel == 1`).
  * **Fold-matched sequence-dark arm ($Q \setminus K$):** $n = 1{,}807$ structural families with $\text{qTM} \geq 0.5$ against at least one reference (`conf_qualified_dual_novel == 0`).
* **Membership Basis (`04_cluster/clusters.tsv` vs. `all_cluster90.tsv`):**
  Upstream scripts (`darkclean_5_mag.sh`, `01_prep.sh`, `04_cluster.sh`) constructed `04_cluster/clusters.tsv` by 3Di-clustering the $\geq 100$ aa triage-no-hit subset of `mag_dark_rep_ids.txt.gz` ($4{,}292{,}446$ **90%-sequence-identity cluster representatives** from `all_cluster90.tsv`, keeping structural clusters with $\geq 3$ 90%-cluster representatives). In `darkclean_5_mag.sh`, a 90%-sequence cluster was included in `mag_dark_rep_ids.txt.gz` if *any* member protein in `all_cluster90.tsv` sat on an HQ+MQ MAG contig, even when the 90%-cluster representative ID itself sat on an unbinned contig.
* **Explicit Conditioning Rule:** This bounded experiment conditions explicitly on the **12,556 verified 90%-sequence-cluster representative members in `04_cluster/clusters.tsv`** across the $2{,}838$ structural families. Inferences apply strictly to this defined representative population and are not generalized to unexpanded `all_cluster90.tsv` copies or unbinned contigs.

### 2.2 Explicit Admission & Eligibility Gates
1. **Gate 0 (Confidence-Qualified Structural Universe):** $2{,}838$ families ($1{,}031$ $K$, $1{,}807$ $Q \setminus K$; $12{,}556$ 90%-cluster representative members on $12{,}382$ unique contigs).
2. **Gate 1 (Refined MAG Binned Subset):** $2{,}426$ families ($917$ $K$, $1{,}509$ $Q \setminus K$) with $\geq 1$ member contig assigned to a refined MAG in `contig_to_bin_refined.tsv.gz` ($6{,}138 / 12{,}382$ contigs binned).
3. **Gate 2 (Valid dRep 95% ANI Species-Cluster Presence):** $2{,}378$ families ($898$ $K$, $1{,}480$ $Q \setminus K$) with $\geq 1$ member in a refined MAG assigned to a non-blank dRep 95% ANI secondary cluster in `Cdb.csv` ($12{,}433$ HQ+MQ refined MAGs across $5{,}438$ species clusters, 100% mapped to GTDB-Tk R232 taxonomy via `Wdb.csv`). Refined bins absent from `Cdb.csv` receive `species = ""` (`binned_unclustered_no_drep_species`) and never form species units.
4. **Gate 3 (Multi-Species Admission Gate):** $1{,}461$ families ($570$ $K$, $891$ $Q \setminus K$) with members in $\geq 2$ distinct dRep 95% ANI species clusters.
   * **Deterministic Expanded Pilot Cohort:** $300$ families ($150$ $K$, $150$ $Q \setminus K$) sampled deterministically from Gate 3 stratified by `(r232_phylum, species_breadth_bin, family)` (`seed = 20261013`).
   * **Full Multi-Species Cohort:** All $1{,}461$ Gate 3 families (`mode = full`).
5. **Gate 4 (Locus-Level & Family-Level Callability Tiers — Evaluated Before KO Aggregation):**
   * **Tier 1 (Primary Gate — Complete 4-Flank Exact-Coordinate Bakta-Mapped Windows):**
     A species representative locus has `eligible_locus_primary == 1` iff its bin has SHA-256 + log proven contig identity and all 4 flanking CDS positions ($m2, m1, p1, p2$) are present on the contig and matched to Bakta CDSs by exact `(norm_contig, start, stop, strand)` coordinates (`mapped_KO` or `mapped_noKO`). A family is primary-evaluable (`eligible_primary == 1`) iff it has $n_{\text{sp,prim}} \geq 2$ distinct species with `eligible_locus_primary == 1`. Families with $n_{\text{sp,prim}} < 2$ are recorded as `INELIGIBLE / MISSING` (`c_fam_primary = -1`), never counted as `0`.
   * **Tier 2 (Audited Cross-Caller Sensitivity Gate — Complete 4-Flank Exact or 3'-Stop+Strand Bakta-Mapped Windows):**
     A locus has `eligible_locus_stop3_sens == 1` iff all 4 flanks are mapped to Bakta either by exact `(norm_contig, start, stop, strand)` or by a unique same-strand, same-3'-stop-codon Bakta CDS `(norm_contig, stop_3prime, strand)` (`mapped_KO`, `mapped_noKO`, `stop3prime_matched_KO`, `stop3prime_matched_noKO`). Family eligibility requires $n_{\text{sp,stop3}} \geq 2$.
   * **Tier 3 (Partial-Window Mapped Sensitivity Gate):**
     A locus has `eligible_locus_partial_sens == 1` iff $\geq 1$ flank is exact-coordinate Bakta-mapped (`mapped_KO` or `mapped_noKO`), treating only exact-mapped flanks as callable. Family eligibility requires $n_{\text{sp,part}} \geq 2$.

---

## 3. Provenance, Contig Verification, & Outcome-Independent Multiplicity Control

### 3.1 Fail-Closed Raw-Bin-to-Bakta Contig Identity (`C2`/`C4`)
For every refined MAG bin `b`:
1. Parse `output_bins/{parent}.fa.gz`: record `(norm_raw_contig, len_bp, sha256_upper, rc_sha256_upper)`.
2. Parse `bakta/{parent}/{parent}.fna`: record `(bakta_contig_id, len_bp, sha256_upper)`.
3. Parse `bakta/{parent}/{parent}.log`: extract `id=<bakta_contig_id>, orig-id=<norm_raw_contig>`.
4. **Fail-Closed Identity Verification:** A Bakta contig `b_hdr` is mapped to `r_hdr` iff:
   * `r_sha256[r_hdr] == b_sha256[b_hdr]` and `r_len[r_hdr] == b_len[b_hdr]` (exact forward-strand nucleotide SHA-256 equality), AND
   * Either `b_log` explicitly maps `b_hdr -> r_hdr`, or (if `.log` is absent) `r_sha256[r_hdr]` is unique (`count == 1`) in both the raw bin FASTA and Bakta `.fna`.
   Any contig with a hash mismatch, duplicate sequence hash without log disambiguation, or reverse-complement orientation fails closed (`contig_identity_failed`).

### 3.2 Exhaustive 7-Category Flank Audit
Every flanking position $w \in \{m2, m1, p1, p2\}$ of every evaluated member row is classified into exactly one category:
1. `mapped_KO`: Proven contig, exact `(norm_contig, start, stop, strand)` Bakta CDS match, eggNOG KEGG KO present.
2. `mapped_noKO`: Proven contig, exact `(norm_contig, start, stop, strand)` Bakta CDS match, no eggNOG KEGG KO.
3. `stop3prime_matched_KO`: Proven contig, no exact start+stop match, unique `(norm_contig, stop_3prime, strand)` Bakta CDS match, eggNOG KEGG KO present.
4. `stop3prime_matched_noKO`: Proven contig, no exact start+stop match, unique `(norm_contig, stop_3prime, strand)` Bakta CDS match, no eggNOG KEGG KO.
5. `coordinate_mismatch_unmapped`: Proven contig, flanking Pyrodigal CDS exists on contig, no exact or unique 3'-stop Bakta CDS match.
6. `source_file_unavailable`: Unbinned contig (`bin == ""`), bin not in `Cdb.csv`, or bin missing verified Bakta/eggNOG files.
7. `contig_edge_truncated`: Flanking index $j + w \notin [0, L_c - 1]$.

### 3.3 Outcome-Independent Multiplicity Control (`C1`/`C2`)
When a family has multiple member loci in the same valid dRep species cluster (`species != ""`), all candidate loci in that `(family, species)` are evaluated for coordinate/caller completeness and real contig length (`contig_len_bp`) **without inspecting KO annotations**, and ranked deterministically in ascending order of the exact 8-tuple (`scripts/02_extract_neighbourhoods.py:838-847`):
$$\text{rank\_key}(r) = \left(-r.\text{eligible\_locus\_primary}, \; -r.\text{eligible\_locus\_stop3\_sens}, \; -r.\text{n\_exact\_mapped\_flanks}, \; -r.\text{n\_stop3\_or\_exact\_mapped\_flanks}, \; -r.\text{n\_present\_flanks}, \; -r.\text{contig\_len\_bp}, \; -r.\text{total\_genes\_on\_contig}, \; r.\text{member}\right).$$
Exactly one top-ranked locus per `(family, species)` receives `is_primary_locus = 1`.

---

## 4. Estimands, Same-Contig Positional Reference, & Between-Arm Inference (`C1`/`C3`)

### 4.1 Family-Level Estimands (Strictly Closed Over Eligible Loci)
For any eligibility tier $t \in \{\text{prim}, \text{stop3}, \text{part}\}$ with eligible species set $S_F^{(t)}$ ($n_{\text{sp}}^{(t)}(F) = |S_F^{(t)}| \geq 2$):
* Let $\text{KO}_F^{(t)}(k) \subseteq S_F^{(t)}$ be the subset of **tier-$t$-eligible species** in which KEGG KO $k$ appears on a tier-$t$-callable flank of that species' designated primary locus.
* **Binary Recurrence Indicator:** $C_{\text{fam}}^{(t)} = \mathbb{I}\left(\max_{k} |\text{KO}_F^{(t)}(k)| \geq 2\right) \in \{0, 1\}$.
* **Continuous Recurrence Proportion:** $R_{\text{fam}}^{(t)} = \frac{1}{n_{\text{sp}}^{(t)}(F)} \max_{k:\,|\text{KO}_F^{(t)}(k)| \geq 2} |\text{KO}_F^{(t)}(k)| \in [0, 1]$ (and $0$ when $\max_{k} |\text{KO}_F^{(t)}(k)| < 2$).

### 4.2 Part (a): Within-Arm Recurrence vs. Same-Contig Eligible-Position Reference
For each contig $c$ of a Tier-$t$-eligible locus ($t \in \{\text{prim}, \text{stop3}\}$), define the eligible pseudo-focal index set:
$$\mathcal{I}_c^{(t)} = \left\{j \in \{2, \dots, L_c - 3\} : \text{all 4 flanks } \{j-2, j-1, j+1, j+2\} \text{ on contig } c \text{ satisfy Tier-}t \text{ Bakta mapping}\right\}$$
with saved count $M_c^{(t)} = |\mathcal{I}_c^{(t)}| \geq 1$.
In each Monte Carlo iteration $b \in \{1, \dots, B_{\text{pos}}\}$ ($B_{\text{pos}} = 5{,}000$), for every Tier-$t$-eligible locus $r \in \mathcal{L}_F^{(t)}$, sample $j^{(b)} \sim \text{Uniform}(\mathcal{I}_c^{(t)})$ and compute $C_{\text{fam}}^{\text{pos},(b)}$ and $R_{\text{fam}}^{\text{pos},(b)}$.
* **Family-level positional expectation:** $\mu_F^{\text{pos}} = \frac{1}{B_{\text{pos}}} \sum_{b=1}^{B_{\text{pos}}} C_{\text{fam}}^{\text{pos},(b)}$, $\rho_F^{\text{pos}} = \frac{1}{B_{\text{pos}}} \sum_{b=1}^{B_{\text{pos}}} R_{\text{fam}}^{\text{pos},(b)}$.
* **Within-arm positional enrichment ($A \in \{K, Q \setminus K\}$):**
  $$\Delta_{\text{pos}}(A) = \bar{C}_A^{\text{obs}} - \bar{\mu}_A^{\text{pos}}, \qquad p_{\text{pos}}(A) = \frac{1 + \sum_{b=1}^{B_{\text{pos}}} \mathbb{I}\left(|\bar{C}_A^{\text{pos},(b)} - \bar{\mu}_A^{\text{pos}}| \geq |\bar{C}_A^{\text{obs}} - \bar{\mu}_A^{\text{pos}}|\right)}{1 + B_{\text{pos}}}$$

### 4.3 Part (b): Between-Arm Contrast on Observed & Background-Adjusted Excess Recurrence
Define each eligible family's **background-adjusted excess recurrence**:
$$E_{\text{fam}} = C_{\text{fam}}^{\text{obs}} - \mu_F^{\text{pos}} \in [-1, 1], \qquad E_{\text{fam}}^{(R)} = R_{\text{fam}}^{\text{obs}} - \rho_F^{\text{pos}} \in [-1, 1]$$
To test whether arm $K$ differs from arm $Q \setminus K$ while conditioning on coarse host lineage, species breadth, and annotation opportunity (under the identifying assumptions of conditional exchangeability within strata/matches and cross-family independence; 4-group R232 phylum stratification and exact phylum-group matching provide coarse taxonomic conditioning rather than full phylogenetic covariance control):
1. **Covariate-Stratified Permutation Test ($B_{\text{perm}} = 10{,}000$, plus-one corrected):**
   * Partition eligible families into joint strata $h = (\text{species\_breadth\_bin}, \text{phylum\_group}, \text{ko\_opportunity\_bin})$, where $\text{species\_breadth\_bin} \in \{2, \geq 3\}$, $\text{phylum\_group} \in \{\text{Thermoplasmatota}, \text{Actinomycetota}, \text{multi}, \text{other}\}$, and $\text{ko\_opportunity\_bin} \in \{\text{low}, \text{high}\}$ (median split of $A_{\text{flank}}$).
   * Audit and report common support ($n_{h,K} \geq 1$ and $n_{h,Q\setminus K} \geq 1$).
   * Within each common-support stratum $h$, permute arm labels $(K, Q \setminus K)$ preserving $(n_{h,K}, n_{h,Q\setminus K})$ to test the Mantel-Haenszel stratum-weighted excess difference $\Delta_{\text{strat}}(E_{\text{fam}})$ and raw risk difference $\Delta_{\text{strat}}(C_{\text{fam}})$.
2. **1:1 Covariate-Matched Pair Analysis:**
   * Match $K$ and $Q \setminus K$ families 1:1 without replacement within species-breadth and lineage blocks minimizing standardized Mahalanobis distance on $(n_{\text{sp}}, A_{\text{flank}}, \text{length}, \text{mean\_pLDDT})$; report pre/post-matching SMDs and paired sign-flip permutation $p$-values ($B = 10{,}000$, plus-one corrected) on $E_{\text{fam}}$ and $C_{\text{fam}}$.
3. **Validated Multivariable Regression Models:**
   * Fit (i) OLS/HC3 regression of $E_{\text{fam}}$ and $C_{\text{fam}}$ on `is_k`, `n_species`, `ko_annotation_opportunity`, `length_100aa`, `mean_plddt`, and supported phylum indicators ($\geq 5$ families with outcome variation), and (ii) multivariable logistic regression (plus L2/Firth-penalized sensitivity check) for $C_{\text{fam}}$, enforcing full column rank, `converged == True`, and finite bounded standard errors (`0 < SE < 100`).
4. **Unadjusted Baselines:** Two-sided Fisher's exact test on $C_{\text{fam}}$ and Mann-Whitney $U$ test on $R_{\text{fam}}$ and $E_{\text{fam}}$.
