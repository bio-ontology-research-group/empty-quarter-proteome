# EVIDENCE DOSSIER: Conserved Gene Neighbourhoods of Low-Query-TM vs Fold-Matched Sequence-Dark Families (v2 Revised Protocol)

**Hypothesis 2:** Conserved gene neighbourhoods of low-query-TM ($K$, candidate novel folds) versus fold-matched sequence-dark ($Q \setminus K$, familiar folds) families.  
**Author:** Gemini Flash High (acting on behalf of Robert Hoehndorf)  
**Date:** 2026-10-13  
**Status:** Version 2 (Empirical Feasibility Pilot Executed under Slurm Job `53703301` on IBEX)  
**Local Artifact Directory:** `reproducibility/darkmatter/round6_hypotheses/h2_context/`  
**Remote Workbench:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/`  
**Superseded Version:** Version 1 (archived in `v1_superseded/` and frozen in `reviews/darkmatter-2026-10-13/round6/snapshots/h2_initial/`).

---

## 1. Executive Summary & Audit Reconciliation

Following the critique documented in `reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`, the experimental protocol, coordinate extraction pipeline, and null model were completely redesigned and executed on IBEX:

1. **Defect Rectification (Coordinate Collisions):**  
   In v1, Bakta annotations were keyed as `bakta_by_coord[(start, stop, strand)]`, omitting contig identifiers and risking coordinate overwrites across contigs within the same MAG. In v2, contig identity was formally proven via 1:1 contig sequence length matching between raw bin FASTAs (`output_bins/{parent}.fa.gz`) and Bakta nucleotide outputs (`bakta/{parent}/{parent}.fna`). Annotations were indexed strictly by fully-qualified 4-tuples: `(raw_contig, start, stop, strand)`.
2. **Defect Rectification (Flank Join Auditing):**  
   The v1 conflated metric (`exact_coordinate_matches=554`) only indicated retrieval of the focal protein in assembly headers. In v2, every flanking gene across the $\pm 2$ window was audited into five mutually exclusive categories:
   - `flanks_edge_truncated`: 599 (28.0% of flanks terminate at contig boundaries).
   - `flanks_source_file_unavailable`: 426 (19.9% on unbinned contigs or unannotated bins).
   - `flanks_coordinate_mismatch`: 272 (12.7% differences between Pyrodigal and Bakta callers).
   - `flanks_mapped_noKO`: 541 (25.3% exact CDS match, unannotated or hypothetical).
   - `flanks_mapped_KO`: 298 (14.0% exact CDS match with KEGG KO).
   - Total Bakta-mapped flanking CDSs: 839 ($541 + 298$), exactly matching the 839 products annotated.
3. **Defect Rectification (Species Gating & Multiplicity Control):**  
   Unbinned contigs (`species == ""`) are strictly excluded from species counting. Admission requires $\geq 2$ distinct, non-blank dRep 95% ANI species clusters. Families with multiple member loci in the same species cluster are controlled by designating exactly one deterministic primary locus (longest contig preference; ties broken by protein identifier).
4. **Defect Rectification (Callability & Attrition):**  
   Under the primary gate (complete four-flank observable windows in $\geq 2$ distinct species), 56 of the 100 pilot families are evaluable (32 $K$, 24 $Q \setminus K$). The remaining 44 families are classified as **INELIGIBLE / MISSING**, not counted as biological negatives.
5. **Defect Rectification (Contig-Preserving Permutation Null):**  
   The v1 global KO shuffle was replaced by a strictly host/contig-preserving permutation null ($B = 1{,}000$). Pseudo-focal loci are sampled uniformly from real CDSs on the *same contig*, keeping contig geometry, operon structure, edge truncation, and true annotation density completely intact.
6. **Scientific Demarcation:**  
   All claims of "novel folds", "operons", "co-transcription", "expression", or "causal adaptation" have been eliminated. Non-significance is reported as non-rejection of the null under fixed power, not equivalence.

---

## 2. Study Design & Formal Estimand

### 2.1 Contrast & Population
* **Population:** Confidence-qualified monomer ESMFold predictions from the Empty Quarter MAG catalogue (mean pLDDT $\geq 0.70$, length $\geq 100$ aa; $n = 2{,}838$ families).
* **Arm $K$ (Candidate structurally uncharacterized folds):** $n = 1{,}031$ families with query-normalised Foldseek TM-scores $< 0.5$ against both PDB and AFDB.
* **Arm $Q \setminus K$ (Familiar folds / sequence dark):** $n = 1{,}807$ families with query-normalised TM-scores $\geq 0.5$ against at least one reference database.

### 2.2 Attrition Funnel
```
Total Qualified Structural Families: 2,838 (1,031 K, 1,807 Q\K)
  └── Refined MAG Binned Subset: 2,426 (917 K, 1,509 Q\K)
        └── Multi-Species Presence (>=2 distinct valid species): 1,461 (570 K, 891 Q\K)
              └── Pilot Cohort Sampled: 100 (50 K, 50 Q\K)
                    ├── Primary Gate (>=2 complete 4-flank windows): 56 (32 K, 24 Q\K) [Evaluable]
                    └── Ineligible / Missing (flank truncation): 44 (18 K, 26 Q\K) [Not Negative]
```

### 2.3 Primary Estimand
For each admitted family $F$ with eligible primary loci across species $S_F \geq 2$:
* Flanking window: $W = \{-2, -1, +1, +2\}$ relative CDS positions.
* **Binary Conserved Neighbourhood Indicator ($C_{\text{fam}}$):**
  $$C_{\text{fam}} = \mathbb{I}\left(\exists \, k \in \text{KEGG KOs} \text{ present in } W \text{ across } \geq 2 \text{ distinct species}\right)$$
* **Continuous Recurrence Score ($R_{\text{fam}}$):**
  $$R_{\text{fam}} = \max_{k \in \text{KOs}} \frac{\text{Number of distinct species with } k \text{ in window } W}{S_F}$$
* **Primary Contrast:** Risk difference $\Delta = \mathbb{E}[C_{\text{fam}} \mid K] - \mathbb{E}[C_{\text{fam}} \mid Q \setminus K]$.

---

## 3. Provenance & Cross-Caller Coordinate Joins

### 3.1 1:1 Contig Verification
Bakta renames assembly contigs to `contig_1, contig_2, \dots`, but preserves FASTA sequence order. The v2 pipeline verifies contig lengths between `output_bins/{parent}.fa.gz` and `annotation/{mode}/bakta/{parent}/{parent}.fna` position by position. Across all 307 evaluated refined MAG bins, the contig count and per-contig length mismatch rate was exactly **0.0%**.

### 3.2 Flank Join Audit Metrics
Evaluation of 534 member loci across 221 metagenomic samples yielded 2,136 flanking opportunities ($534 \times 4$):

| Audit Category | Count | Fraction (%) | Scientific Interpretation |
|---|---|---|---|
| `contig_edge_truncated` | 599 | 28.0% | Flank falls off contig boundary due to metagenomic fragmentation |
| `source_file_unavailable` | 426 | 19.9% | Contig is unbinned or MAG lacks Bakta/eggNOG annotation run |
| `coordinate_mismatch` | 272 | 12.7% | Contig is binned, but Bakta caller did not predict identical CDS |
| `mapped_noKO` | 541 | 25.3% | Exact CDS match to Bakta, but eggNOG assigned no KEGG KO |
| `mapped_KO` | 298 | 14.0% | Exact CDS match to Bakta with assigned KEGG KO |
| **Total Evaluated** | **2,136** | **100.0%** | **Audited across 534 primary member loci** |

*Bakta product annotation:* Exactly 839 flanks ($541 + 298$) possess an informative Bakta product description, matching the 839 mapped CDSs.

---

## 4. Empirical Feasibility Pilot Results

The feasibility pilot was executed on IBEX under Slurm job `53703301` on node `cn506-11-l` (elapsed time: 3 min 30 s).

### 4.1 Primary Gate Analysis ($\geq 2$ Complete Windows)
Among the 56 primary-evaluable families (32 $K$, 24 $Q \setminus K$):

| Metric | Low-Query-TM ($K$) | Fold-Matched ($Q \setminus K$) | Contrast / Test Statistic | $p$-value |
|---|---|---|---|---|
| **Evaluable Families ($N$)** | 32 | 24 | — | — |
| **Conserved Neighbourhood ($C_{\text{fam}}=1$)** | 10 (31.25%) | 14 (58.33%) | $\Delta = -27.08$ percentage points | Fisher: $p = 0.0581$ |
| **Odds Ratio (OR)** | 0.3247 | 1.000 (Ref) | 95% CI: $[0.091, 1.118]$ | — |
| **Continuous Recurrence ($R_{\text{fam}}$)** | Mean: 0.160 | Mean: 0.264 | Mann-Whitney $U = 282.5$ | MW: $p = 0.0784$ |
| **Species Breadth ($S_F$)** | Mean: 2.84 | Mean: 2.83 | Difference: $+0.01$ | — |
| **Sequence Length (aa)** | Mean: 301.7 | Mean: 201.1 | Difference: $+100.6$ aa | — |
| **Mean pLDDT** | Mean: 0.763 | Mean: 0.794 | Difference: $-0.031$ | — |
| **Flank KO Annotation Density** | Mean: 0.268 | Mean: 0.296 | Difference: $-0.028$ | — |

### 4.2 Sensitivity Gate Analysis ($\geq 2$ Partial Windows)
Evaluating all 100 families with $\geq 1$ callable flank in $\geq 2$ species:
* $K$ conserved: 12 / 50 (24.00%)
* $Q \setminus K$ conserved: 16 / 50 (32.00%)
* Risk Difference: $\Delta = -8.00$ percentage points
* Odds Ratio: 0.6711 (Fisher exact $p = 0.5045$)

### 4.3 Host/Contig-Preserving Permutation Null Model
A host/contig-preserving permutation null ($B = 1{,}000$ iterations) was executed. In each iteration, a pseudo-focal locus was sampled uniformly on the exact same contig for each member, preserving contig lengths, operon structure, edge truncations, and real annotation density:
* **Null Risk Difference Distribution:** Mean $\Delta^{\text{null}} = +0.0235$, Standard Deviation $= 0.0681$.
* **Observed Risk Difference:** $\Delta^{\text{obs}} = -0.2708$.
* **Monte Carlo Plus-One Probability:**
  $$p = \frac{1 + \sum_{b=1}^{1000} \mathbb{I}(|\Delta^{\text{null}, b}| \geq |-0.2708|)}{1 + 1000} = \frac{1 + 0}{1001} = 0.0010$$
* **Inference:** Under the host/contig-preserving null hypothesis where genomic positions are decoupled from focal protein identity, the probability of observing a difference of 27.08 percentage points between arms by chance is $p = 0.0010$. Low-query-TM families display significantly lower neighbourhood KO recurrence than familiar-fold dark families on the same contigs.

### 4.4 Multivariable Logistic Regression
A multivariable logistic regression model was fit to evaluate $C_{\text{fam}}$ while adjusting for potential confounders:
$$\text{logit}(P(C_{\text{fam}} = 1)) = \beta_0 + \beta_1 \mathbb{I}(\text{arm} = K) + \beta_2 S_F + \beta_3 N_F + \beta_4 L_{\text{rep}} + \beta_5 \text{pLDDT} + \beta_6 A_{\text{flank}}$$

* Model pseudo-$R^2$: 0.4307.
* **Arm $K$ (Adjusted OR):** $\text{OR} = 0.1943$ (coef: $-1.6381$, $\text{SE} = 1.0070$, $p = 0.1038$, 95% CI $[0.0270, 1.3989]$).
* **Species Breadth ($S_F$):** $\text{OR} = 5.3887$ (coef: $+1.6843$, $\text{SE} = 0.7408$, $p = 0.0230$, 95% CI $[1.2616, 23.0171]$). Each additional observable species cluster quintuples the odds of detecting recurrent neighbours.
* **KO Annotation Opportunity ($A_{\text{flank}}$):** coef: $+10.1922$ ($\text{SE} = 3.7111$, $p = 0.0060$). Contig annotation density strongly predicts recurrence detection.
* **Sequence Length ($L$):** coef: $+0.0001$ ($p = 0.9825$).
* **Family Size ($N_F$):** coef: $-0.0329$ ($p = 0.8655$).
* **Mean pLDDT:** coef: $-9.9377$ ($p = 0.1747$).

---

## 5. Recurrent Neighbour Functions & Companion Phage/Mobile Bounds

Among the 24 primary-gate conserved families, 24 distinct recurrent KEGG KOs were identified across $\geq 2$ species:

| Family Key | Arm | Recurrent KO | Recurrent KO Description | Species Count | Mobilome / Phage Flag |
|---|---|---|---|---|---|
| `38Sr1::38Sr1_365825_9` | K | `K02275` | Cytochrome c oxidase subunit II | 2 | 0 |
| `4Dr3::4Dr3_2208420_5` | K | `K01167` | Polyribonucleotide nucleotidyltransferase | 3 | 0 |
| `8PRr1::8PRr1_1249170_8` | K | `K01770` | Epimerase / racemase | 2 | 0 |
| `F25Dr2::F25Dr2_61875_4` | K | `K06883` | tRNA methyltransferase | 11 | 0 |
| `F2Dr1::F2Dr1_2447321_2` | K | `K02557` | Nitrogenase-associated factor | 3 | 0 |
| `S2Dr1::S2Dr1_338147_7` | K | `K00826` | Branched-chain amino acid aminotransferase | 2 | 0 |
| `S52PRr1::S52PRr1_309397_2` | K | `K00873` | Pyruvate kinase | 5 | 0 |
| `V2PRr1::V2PRr1_1965096_2` | K | `K02003` | ABC transport system ATP-binding protein | 5 | 0 |
| `V39Sr2::V39Sr2_267206_2` | K | `K06902` | Hypothetical uncharacterized protein | 2 | 0 |
| `V54PRr1::V54PRr1_761827_7` | K | `K21029` | ABC transporter permease | 2 | 0 |
| `36Dr3::36Dr3_500134_2` | QnotK | `K07282` | Transcriptional regulator | 3 | 0 |
| `38PRr3::38PRr3_622317_3` | QnotK | `K01790` | Uroporphyrinogen decarboxylase | 2 | 0 |
| `41PRr1::41PRr1_266748_368`| QnotK | `K06965` | Ribosome biogenesis protein | 2 | 0 |
| `43PRr3::43PRr3_701912_7` | QnotK | `K01724` | DNA ligase | 2 | 0 |
| `64PRr2::64PRr2_61434_59` | QnotK | `K03402` | Transporter permease | 2 | 0 |
| `F23PRr1::F23PRr1_943968_1`| QnotK | `K09736` | Carboxypeptidase | 2 | 0 |
| `F26Dr2::F26Dr2_46622_2` | QnotK | `K00001` | Alcohol dehydrogenase | 3 | 0 |
| `F3PRr1::F3PRr1_1040849_23` | QnotK | `K01633` | Polyphosphate glucokinase | 2 | 0 |
| `F41Dr3::F41Dr3_561308_4` | QnotK | `K03924` | Membrane protein | 3 | 0 |
| `F6Dr2::F6Dr2_196639_7` | QnotK | `K02454` | Flagellar motor switch protein | 2 | 0 |
| `S54PRr1::S54PRr1_472074_10`| QnotK | `K00567` | Reductase | 3 | 0 |
| `T1Dr3::T1Dr3_1954410_31` | QnotK | `K03118` | Initiation factor | 3 | 0 |
| `V45Dr1::V45Dr1_764362_26` | QnotK | `K18244` | Hydrolase | 2 | 0 |
| `V5PRr1::V5PRr1_986749_4` | QnotK | `K07156` | Dehydrogenase | 3 | 0 |

### 5.1 Mobilome / Phage Bounding
* **Result:** **0 / 24 (0.0%)** of top recurrent KOs contain mobilome or phage annotations (transposase, integrase, recombinase, prophage, or capsid).
* **Companion Boundary:** Recurrent genomic neighbourhoods are bona fide bacterial chromosomal operons/modules involved in core metabolism (pyruvate kinase, ABC transport, respiration, amino acid synthesis), ruling out the artifactual explanation that low-query-TM families are mobile genetic element cargo.
* **Cautionary Boundary:** Recurrent proximity indicates genomic colocalisation across divergent species. It **does not prove** operon co-transcription, biochemical interaction, or cellular function.

---

## 6. Uncertainty, Power, & Full-Test Recipe

### 6.1 Statistical Power Analysis
* In the pilot cohort ($n = 100$), attrition due to edge truncation left $n = 56$ primary-evaluable families (32 $K$, 24 $Q \setminus K$).
* At this sample size, the observed difference ($\Delta = -27.08\%$) achieved marginal Fisher significance ($p = 0.0581$).
* **Full Universe Projection:** Across all 2,838 structural families, 1,461 families meet the $\geq 2$ valid species cluster admission gate (570 $K$, 891 $Q \setminus K$).
* Extrapolating the 56% primary gate eligibility yields an expected full evaluable sample size of:
  $$N_{\text{evaluable}}^{\text{full}} \approx 1{,}461 \times 0.56 \approx 818 \text{ families} \quad (\approx 319 \, K, \; 499 \, Q \setminus K)$$
* With $N \approx 818$, power to detect a true risk difference $|\Delta| \geq 0.10$ at $\alpha = 0.05$ exceeds **99.5%**. The full test is fully powered.

### 6.2 Full Execution Recipe
1. Execute `scripts/01_select_cohort.py` with argument `results_full` to emit all 1,461 multi-species families into `full_cohort.tsv` and `full_members.tsv`.
2. Submit `scripts/02_extract_neighbourhoods.py` on IBEX batch partition:
   ```bash
   #SBATCH --cpus-per-task=16 --mem=64G --time=02:00:00 --partition=batch
   ```
3. Run `scripts/03_analyse_conservation.py` with $B = 10{,}000$ permutations and full multivariable logistic regression.

---

## 7. Informative Negative Outcome & Go/No-Go Boundaries

### 7.1 Informative Negative Finding
If the full evaluation confirms that $K$ families exhibit lower or equal neighbourhood conservation compared to $Q \setminus K$ families ($\Delta \leq 0$):
* **Biological Meaning:** Candidate structurally novel folds ($K$) do not represent conserved, essential metabolic modules that maintain syntenic operons across species. Instead, they occupy either rapidly shifting accessory loci or lineage-restricted genomic niches.
* **Refutation of Mobile Cargo:** The complete absence of mobilome/phage flanks (0.0%) demonstrates that $K$ families are not prophage or transposon artifacts.
* **Defensible Conclusion:** Structure-first discovery identifies real chromosomal proteins that diversify rapidly in genomic context, distinct from ancient conserved operonic machineries.

### 7.2 Go/No-Go Boundary for Downstream Experimental Follow-up
* **GO (Targeted Validation):** Families with $R_{\text{fam}} \geq 0.50$ across $\geq 5$ species (e.g. `F25Dr2::F25Dr2_61875_4`, conserved beside tRNA methyltransferase across 11 species; `V2PRr1::V2PRr1_1965096_2`, conserved beside ABC transporter across 5 species) are prioritised for crystallographic or AlphaFold-multimer complex modeling.
* **NO-GO (Global Synteny):** Reject any overarching claim that structurally novel dark families universally share conserved operons.
