# DESIGN: Conserved Gene Neighbourhoods of Low-Query-TM vs Fold-Matched Sequence-Dark Families (v2 Revised Protocol)

**Hypothesis 2:** Conserved gene neighbourhoods of low-query-TM versus fold-matched sequence-dark families.  
**Author:** Gemini Flash High (acting on behalf of Robert Hoehndorf)  
**Date:** 2026-10-13  
**Status:** Version 2 (Revised Protocol addressing Root Audit Notes; pre-computation specification)  
**Local Workspace:** `reproducibility/darkmatter/round6_hypotheses/h2_context/`  
**Remote Workbench:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/`  
**Superceded:** Version 1 (archived in `v1_superseded/` and frozen in `reviews/darkmatter-2026-10-13/round6/snapshots/h2_initial/`).

---

## 1. Audit Reconciliation & Required Corrections

Following review of `reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`, the study design is revised to eliminate four specific defects in the v1 pilot:

| Defect in v1 | Root Cause in v1 Code | Mandated Correction in v2 |
|---|---|---|
| **Contig-ID Dropping** | `bakta_by_coord[(start, stop, strand)]` dropped contig identity; coordinates collision across contigs in same bin. | Key annotations strictly by `(sample, raw_contig, start, stop, strand)`. Prove 1:1 contig mapping from raw bin to Bakta FASTA by length matching; fail closed on mismatch. |
| **Conflated Join Counter** | `exact_coordinate_matches=554` only tracked finding the focal protein in assembly FAA; did not audit cross-caller joins. | Formally audit all flanking CDS joins into 5 mutually exclusive categories: (1) `mapped_KO`, (2) `mapped_noKO`, (3) `coordinate_mismatch`, (4) `source_file_unavailable`, and (5) `contig_edge_truncated`. |
| **Species Callability & Blank Species** | Unbinned contigs (`species == ""`) were not strictly prevented from acting as distinct species; partial windows counted toward multi-species callability. | Only mapped, non-blank dRep species clusters count toward multi-species presence. Define primary callability on complete four-flank observable windows; partial windows reported only as sensitivity. |
| **Multiplicity within Species** | Multiple members in same species were aggregated without explicit multiplicity control. | For families with multiple copies per species, deterministically select the primary locus with the longest contig / complete window, and record copy number explicitly. |
| **Global Shuffling in Null** | `03_analyse_conservation.py` pooled all KOs globally across all families and shuffled, destroying contig/host structure. | Implement a strictly host/contig-preserving random focal-position permutation null: pseudo-focal loci are sampled on the *same contig*, inheriting exact contig geometry and real annotations. Use Monte Carlo plus-one correction. |
| **Overstated Claims** | Prose discussed novel folds, operons, and $p > 0.05$ equivalence. | Remove all novel-fold, operon, expression, and equivalence claims. Non-significance is reported as non-rejection of the null under fixed power. |

---

## 2. Formal Hypothesis Specification

### 2.1 Research Question & Contrast
* **Biological Question:** Do low-query-TM families repeatedly occur beside the same annotated neighbours across distinct species, compared with fold-matched sequence-dark families, after controlling for sequence length, member count, species breadth, contig fragmentation, and host lineage?
* **Population:** Confidence-qualified monomer ESMFold predictions from the Empty Quarter MAG catalogue (mean pLDDT $\geq 0.70$, length $\geq 100$ aa; $n = 2{,}838$ families).
  * **Low-query-TM arm ($K$):** $n = 1{,}031$ families with query-normalised Foldseek TM-scores $< 0.5$ against both PDB and AFDB (dual no-hits; candidate structurally uncharacterized folds).
  * **Fold-matched sequence-dark arm ($Q \setminus K$):** $n = 1{,}807$ families with query-normalised TM-scores $\geq 0.5$ against at least one reference (known/familiar fold despite sequence darkness).

### 2.2 Primary Unit of Analysis
The **gene family** ($F$) is the primary unit of analysis. Member proteins within a family share evolutionary history and cannot be treated as independent units (anti-pseudoreplication rule).

### 2.3 Attrition & Admission Gate
Families are narrowed through explicit, recorded gates:
1. **Total Structural Universe:** $n = 2{,}838$ (1,031 $K$, 1,807 $Q \setminus K$).
2. **Refined MAG Binned Subset:** Families with $\geq 1$ member on a contig assigned to a post-refinement MAG (`contig_to_bin_refined.tsv.gz`).
3. **Multi-Species Presence:** Families with admitted members binned into $\geq 2$ distinct, non-blank dRep 95% ANI species clusters (`Cdb.csv`).
4. **Primary Observable Window Gate:** Families with $\geq 2$ distinct species exhibiting a **complete four-flank observable window** (no contig-edge truncation among flanks $m2, m1, p1, p2$).
   * *Ineligible rule:* Families with $< 2$ species meeting the complete window gate are classified as **INELIGIBLE / UNRESOLVED**, not counted as biological negatives.
   * *Sensitivity Gate:* Families with $\geq 2$ species having at least one callable flank on each side (partial window).

### 2.4 Primary Estimand
Let $F$ be an admitted family with eligible species $S_F \geq 2$. For each species $s \in S_F$, exactly one deterministic representative locus is evaluated (the locus on the longest contig).
* **Flanking Window:** $W = \{-2, -1, +1, +2\}$ relative CDS positions along the contig.
* **Family-Level Conserved Neighbourhood Indicator ($C_{\text{fam}}$):**
  $$C_{\text{fam}} = \mathbb{I}\left(\exists \, k \in \text{KEGG KOs} \text{ such that } k \text{ occurs in } W \text{ across } \geq 2 \text{ distinct species}\right)$$
* **Continuous Recurrence Score ($R_{\text{fam}}$):**
  $$R_{\text{fam}} = \max_{k \in \text{KOs}} \frac{\text{Number of distinct species with } k \text{ in flank } W}{S_F}$$
* **Primary Contrast:** Risk difference $\Delta = \mathbb{E}[C_{\text{fam}} \mid K] - \mathbb{E}[C_{\text{fam}} \mid Q \setminus K]$, accompanied by adjusted odds ratios from multivariable logistic regression.

---

## 3. Provenance & Cross-Caller Coordinate Joins

### 3.1 Source-to-Bakta Contig Mapping
Bakta preserves input contig ordering but renames contigs to `contig_1, contig_2, \dots`.
1. Stream `output_bins/{parent}.fa.gz` to record `R = [(raw_contig_name, length_bp)]` in file order.
2. Stream `annotation/{mode}/bakta/{parent}/{parent}.fna` to record `K = [(bakta_contig_name, length_bp)]` in file order.
3. **Fail-Closed Verification:** Assert $\text{len}(R) == \text{len}(K)$ and $R[i][1] == K[i][1]$ for all $i$. If any length differs or count mismatches, refuse the bin.
4. Mapping: $\text{bakta\_to\_raw}[K[i][0]] = R[i][0]$.

### 3.2 Annotation Join & Audit Categories
Bakta annotations from `annotation/{mode}/bakta/{parent}/{parent}.tsv` are mapped to original contig names and keyed strictly as:
$$\text{key} = (\text{sample}, \text{raw\_contig}, \text{start}, \text{stop}, \text{strand})$$
Every flanking CDS is evaluated and placed into exactly one audit category:
1. `mapped_KO`: Contig in bin, exact start/stop/strand match to Bakta CDS, eggNOG assigned a KEGG KO.
2. `mapped_noKO`: Exact match to Bakta CDS, but eggNOG assigned no KO (e.g. unannotated or hypothetical).
3. `coordinate_mismatch`: Contig in bin, but no Bakta CDS shares exact start/stop/strand.
4. `source_file_unavailable`: Contig unbinned, or bin annotation files missing.
5. `contig_edge_truncated`: Contig terminates before this flanking position can be occupied.

---

## 4. Null Model: Host/Contig-Preserving Random Position Permutation

To test whether observed neighbourhood recurrence exceeds background expectation without breaking genomic architecture:
* For each contig harbouring an admitted member, determine the total number of called CDSs $L$.
* In each permutation $b \in \{1, \dots, B\}$:
  * Draw a random pseudo-focal locus $j \in \{0, \dots, L-1\}$ uniformly on the **same contig**.
  * The pseudo-flank window is $\{j-2, j-1, j+1, j+2\}$ on that same contig.
  * Flanks inherit real annotations and real edge truncations from their actual positions on that contig.
  * Evaluate whether any KO is shared across $\geq 2$ distinct species under this permuted configuration.
* Compute the empirical Monte Carlo $p$-value with the **plus-one correction**:
  $$p = \frac{1 + \sum_{b=1}^B \mathbb{I}(|\Delta^{\text{null}, b}| \geq |\Delta^{\text{obs}}|)}{1 + B}$$

---

## 5. Covariate Conditioning & Statistical Power

### 5.1 Covariate Adjustment
In addition to stratified 2x2 comparison, multivariable logistic regression models:
$$\text{logit}(P(C_{\text{fam}} = 1)) = \beta_0 + \beta_1 \mathbb{I}(\text{arm} = K) + \beta_2 S_F + \beta_3 N_F + \beta_4 L_{\text{rep}} + \beta_5 \text{pLDDT} + \sum_p \gamma_p \mathbb{I}(\text{phylum} = p) + \beta_6 A_{\text{flank}}$$
where $A_{\text{flank}}$ is the annotation opportunity (fraction of callable flanks with an assigned KO).

### 5.2 Scientific Boundaries
* **Candidate Hypotheses Only:** Recurrent neighbours generate specific candidate pathway associations for targeted follow-up; they do not prove operon co-transcription, biochemical activity, or cellular function.
* **No Causal Adaptation:** Proximity to stress or transport genes cannot establish environmental adaptation or positive selection.
* **Informative Negative Outcome:** If $\Delta \approx 0$ with narrow confidence intervals, low-query-TM families are confirmed to occupy standard chromosomal neighbourhoods at rates identical to familiar-fold sequence-dark families, ruling out preferential displacement into mobile islands or spurious non-coding artifacts.
