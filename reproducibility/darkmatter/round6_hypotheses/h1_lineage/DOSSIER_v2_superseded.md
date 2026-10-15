# EVIDENCE DOSSIER: Hypothesis 1 — What Explains Within-Thermoplasmatota Darkness?

**Author:** Gemini 3.8 Flash (High) via agy pair programming  
**Date:** 2026-10-13  
**Status:** Completed v2 Rigorous Phylogenetic Covariance & Reconciled Association Analysis  
**Local Artifact Directory:** `reproducibility/darkmatter/round6_hypotheses/h1_lineage/`  
**Remote Workbench Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h1_lineage/`  
**Remote Slurm Job Execution:**
- **Job ID `53701040` (v2 Confirmatory Analysis):** Partition `debug`, node `cn506-11-l`, 8 CPUs, 32 GB RAM, runtime 24.19s, exit code 0.
- **Job ID `53700584` (v1 Exploratory Pilot):** Superseded initial pilot, preserved in logs.

---

## 1. Exact Hypothesis Formulation

**Context & Original Question (Round 5 EXTENSIONS.md):**  
Within the Rub' al-Khali desert metagenome-assembled genome (MAG) catalogue, the archaeal phylum *Thermoplasmatota* displays an exceptionally high dark proteome fraction (median 30.96%, ~4.5-fold enriched in dark families relative to species count). What explains within-*Thermoplasmatota* darkness? Does darkness track external amino acid divergence to reference genomes, shared-gene coverage, or lineage-specific gene content after accounting for genome quality and phylogenetic ancestry?

**Specific Refined Hypothesis (H1):**  
Across the 237 desert *Thermoplasmatota* species representatives, variation in the dark proteome fraction ($Y_i = N_{\mathrm{dark}, i} / N_{\mathrm{prot}, i}$) is associated with shared external gene coverage ($N_{\mathrm{RBH}} / N_{\mathrm{desert}}$), assembly length, and contig fragmentation, but exhibits no negative dilution with continuous external reference proximity (AAI) once conditioned on phylogenetic branch-length covariance on the ar53 marker tree.

**Explicit Scope and Epistemic Boundaries:**  
- This is a **finite-cohort association study** conditioned on the observed 237 species representatives and their placement in the GTDB-Tk r232 ar53 tree.
- It does **not** assert universal macroevolutionary rules across the tree of life.
- It does **not** assert macroevolutionary directionality (e.g., accessory "expansion" in `c__SW-10-69-26` vs genome "reduction" in `c__E2`), which requires outgroup-rooted ancestral character state reconstruction.
- It treats gene-caller variation and contig fragmentation as measured quality covariates, not as definitive proof that all unclustered genes are caller artifacts.

---

## 2. Primary Unit of Analysis and Identifiable Estimands

### Primary Unit of Analysis
The species-representative MAG ($i = 1, \dots, 237$) in *Thermoplasmatota*, analyzed within its continuous phylogenetic covariance matrix $\mathbf{C} \in \mathbb{R}^{237 \times 237}$ parameterized from the 1,914-taxon ar53 FastTree (`thermo_tree.nwk`) and conditioned on assembly quality.

### Primary Identifiable Estimands
In the primary Phylogenetic Generalized Least Squares (PGLS) model:

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim \mathcal{N}\left(\mathbf{0}, \sigma^2 [\lambda \mathbf{C} + (1 - \lambda)\mathbf{I}_{237}]\right)$$

with response $\mathbf{y} = \mathrm{dark\_fraction} \in [0, 1]$, the primary estimands are:
1. **External Divergence Slope ($\beta_{\mathrm{AAI}}$):** The expected change in dark fraction per 1 percentage-point increase in nearest-relative AAI, conditional on shared gene coverage, completeness, contamination, $\log_{10}(\mathrm{length})$, $\log_{10}(\mathrm{N50})$, and phylogenetic correlation at profile MLE $\hat{\lambda}_{\mathrm{ML}}$.
2. **Shared RBH Coverage Slope ($\beta_{\mathrm{frac\_rbh\_desert}}$):** The expected change in dark fraction per unit increase in shared RBH coverage ($N_{\mathrm{RBH}} / N_{\mathrm{desert}}$), conditional on AAI, quality, and phylogenetic covariance.

### Secondary Descriptive Estimands
- **Reconciled Ratio vs Count Correlations:** Contrasting the observed ratio correlation against the algebraic common-denominator permutation null and the raw decoupled count correlation ($N_{\mathrm{dark}}$ vs $N_{\mathrm{RBH}}$).
- **Accessory Dark Family Retention Profile:** The descriptive distribution of 14,645 admitted dark families across genome sharing tiers (singletons, rare, intermediate, frequent, core).

---

## 3. Formal Null and Alternative Hypotheses

| Estimand | Null Hypothesis ($H_0$) | Alternative Hypothesis ($H_1$) | Statistical Test & Model | Primary Result |
|---|---|---|---|---|
| **External Reference Proximity ($\beta_{\mathrm{AAI}}$)** | $\beta_{\mathrm{AAI}} \le 0$ under PGLS; higher reference AAI either reduces or has no association with dark fraction. | $\beta_{\mathrm{AAI}} > 0$; conditional on shared coverage and quality, higher AAI within this range tracks slightly higher dark fraction. | PGLS Wald $t$-test ($df = 230$), profile $\hat{\lambda}_{\mathrm{ML}} = 0.9158$. | **$\beta = +0.00170$, $t = +4.25$, $p = 3.12 \times 10^{-5}$** (Rejects negative dilution $H_0$). |
| **Shared RBH Coverage ($\beta_{\mathrm{frac\_rbh}}$)** | $\beta_{\mathrm{frac\_rbh}} = 0$; shared external gene coverage does not predict dark fraction beyond quality and tree structure. | $\beta_{\mathrm{frac\_rbh}} < 0$; genomes sharing a higher fraction of their proteome have lower dark fraction. | PGLS Wald $t$-test ($df = 230$), profile $\hat{\lambda}_{\mathrm{ML}} = 0.9158$. | **$\beta = -0.13748$, $t = -7.17$, $p = 7.50 \times 10^{-13}$** (Significant negative association). |
| **Phylogenetic Signal ($\lambda$)** | $\lambda = 0$; residual variation in dark fraction is phylogenetically independent. | $\lambda > 0$; dark fraction displays significant phylogenetic autocorrelation along branch lengths. | Profile Likelihood Ratio Test: $2(\ln L_{\mathrm{ML}} - \ln L_{\lambda=0})$. | **$\hat{\lambda}_{\mathrm{ML}} = 0.9158$, $\mathrm{LRT} = 229.56$, $p < 10^{-50}$** (Strong phylogenetic signal). |
| **Common-Denominator Independence** | Observed ratio correlation matches the algebraic null of independent numerators divided by common denominator ($N_{\mathrm{prot}}$). | Observed ratio correlation deviates from the independent-count common-denominator null. | Permutation null test with plus-one tail correction: $(k+1)/(N+1)$. | **Observed $\rho = -0.3942$ vs Null mean $\rho = +0.3017$** ($p = 9.999 \times 10^{-5}$; rejects independent-count null). |

---

## 4. Methodological Audit & Scientific Rectifications

### A. Resolution of the Ratio Coupling Sign Error
- **The Initial Draft Mistake:** The exploratory draft incorrectly claimed that the common-denominator permutation null was $\rho = -0.3017$ and asserted that the observed negative ratio correlation was "100% mathematical artifact".
- **The Corrected Mathematics:** Dividing two independent positive random variables by a common variable denominator ($Z$) mechanically creates a **positive** correlation: when $Z$ is large, both $X/Z$ and $Y/Z$ decrease; when $Z$ is small, both ratios increase. In our audited 10,000-iteration permutation test across the observed $N_{\mathrm{prot}}$ distribution, the null expected ratio correlation is **positive: mean $\rho = +0.3017$ (95% CI: $[+0.1789, +0.4194]$)**.
- **The Biological & Structural Reality:** The observed correlation within `SW-10-69-26` is strongly **negative** ($\rho = -0.3942, p = 2.83 \times 10^{-8}$), which decisively departs from the positive independent-count null ($p = 9.999 \times 10^{-5}$). 
- **Why is the observed correlation negative?** Because dark genes and RBH genes are mutually exclusive subsets constrained by the total proteome:
  $$N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$$
  Across these genomes, proteome size variation is heavily driven by dark genes ($\rho = +0.9457, p = 2.82 \times 10^{-91}$ with $N_{\mathrm{prot}}$), whereas external RBH gene counts scale much more slowly ($\rho = +0.3294, p = 4.68 \times 10^{-6}$). As genomes expand, dark genes occupy an increasing fraction of the proteome, naturally reducing the fraction available for conserved RBH genes.
- **Decoupled Raw Counts:** At the raw count level, $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ are **positively correlated** ($\rho = +0.2419, p = 9.09 \times 10^{-4}$). Larger genomes have more of *both* gene classes, but dark genes expand at a much higher relative rate.

### B. RBH Denominator Audit
- In `v2_thermo_relatives_aai.py`, `aai_frac_rbh` was computed as $N_{\mathrm{RBH}} / \min(N_a, N_b)$.
- In **74.7% of pairs (177 / 237)**, the desert MAG has more proteins than its nearest external comparator. Consequently, the published `aai_frac_rbh` used the *external comparator's* proteome size as the denominator.
- Audited values across the 237 representatives:
  * Published $\min$-denominator ratio: median 0.4690 (IQR 0.4050–0.5300).
  * True desert-normalized ratio ($N_{\mathrm{RBH}} / N_{\mathrm{desert}}$): median 0.4092 (IQR 0.3409–0.4934).
  * External-normalized ratio ($N_{\mathrm{RBH}} / N_{\mathrm{nearest}}$): median 0.4564 (IQR 0.3935–0.5245).
- In all v2 models, we strictly use the true desert-normalized ratio `frac_rbh_desert`.

### C. Tree-Based Continuous Covariance vs Inadequate Clade Clustering
- Earlier attempts to adjust for lineage using asymptotic cluster-sandwich standard errors at the Order ($G=3$) or Family ($G=7$) level are statistically invalid because cluster asymptotics fail when $G < 30$ and sample sizes are severely unbalanced (one family, `f__JAJPHT01`, contains 140 / 237 representatives).
- We resolved this by extracting the full $237 \times 237$ cophenetic branch-length covariance matrix $\mathbf{C}$ from the 1,914-taxon FastTree (`thermo_tree.nwk`) using DendroPy 5.0.8, verifying exact tip mapping (237 / 237 tips mapped, 100.0%).
- All ancestry adjustments are parameterized via continuous PGLS with Pagel's $\lambda$ profiling.

---

## 5. Executed Empirical Results (IBEX Slurm Job 53701040)

All results are audited against exact saved outputs in `results/`.

### Table 1: Primary PGLS Model Summary ($n = 237, \hat{\lambda}_{\mathrm{ML}} = 0.9158$)
Model fit: $\ln L = 640.933$, $\mathrm{AIC} = -1265.866$, $\mathrm{BIC} = -1238.121$, $R^2 = 0.99916$, $\sigma^2 = 0.001028$, $df_{\mathrm{resid}} = 230$.

| Covariate | PGLS Coef ($\beta$) | PGLS SE | $t$-statistic | $p$-value | 95% Model-Conditional CI | OLS Coef ($\lambda=0$) | OLS SE | OLS $p$-value |
|---|---|---|---|---|---|---|---|---|
| **Intercept** | -0.45321 | 0.11020 | -4.11 | $5.53 \times 10^{-5}$ | [-0.66921, -0.23721] | -1.25228 | 0.12995 | $< 10^{-15}$ |
| **completeness** | -0.00034 | 0.00012 | -2.70 | $0.00748$ | [-0.00058, -0.00009] | -0.00051 | 0.00020 | $0.01065$ |
| **contamination** | -0.00032 | 0.00048 | -0.66 | $0.50959$ | [-0.00127, +0.00063] | -0.00011 | 0.00087 | $0.90346$ |
| **log10(length)** | +0.11543 | 0.01812 | +6.37 | $1.01 \times 10^{-9}$ | [+0.07992, +0.15094] | +0.24985 | 0.02042 | $< 10^{-15}$ |
| **log10(N50)** | -0.00900 | 0.00299 | -3.02 | $0.00285$ | [-0.01486, -0.00315] | -0.02264 | 0.00476 | $3.46 \times 10^{-6}$ |
| **aai** | +0.00170 | 0.00040 | +4.25 | $3.12 \times 10^{-5}$ | [+0.00092, +0.00249] | +0.00193 | 0.00038 | $9.88 \times 10^{-7}$ |
| **frac_rbh_desert** | -0.13748 | 0.01917 | -7.17 | $7.50 \times 10^{-13}$ | [-0.17505, -0.09992] | -0.10541 | 0.02398 | $1.67 \times 10^{-5}$ |

*Key Inference:* After accounting for continuous phylogenetic covariance and quality covariates:
1. Nearest-relative AAI displays a small, statistically positive slope ($\beta = +0.00170$ per 1% AAI), decisively showing that within this evolutionary distance band (55–72% AAI), higher AAI does *not* dilute dark fraction.
2. Shared RBH coverage maintains a strong negative conditional association ($\beta = -0.13748, p = 7.50 \times 10^{-13}$).
3. Assembly length ($\beta = +0.11543, p = 1.01 \times 10^{-9}$) and contig fragmentation ($\beta = -0.00900, p = 0.00285$) remain highly significant independent predictors.

---

### Table 2: Pagel's $\lambda$ Profile and Signal Sensitivity ($n = 237$)

| Pagel's $\lambda$ | Model Status | Log-Likelihood ($\ln L$) | AIC | BIC | $R^2$ | $\beta_{\mathrm{AAI}}$ (SE) | $p_{\mathrm{AAI}}$ | $\beta_{\mathrm{frac\_rbh}}$ (SE) | $p_{\mathrm{frac\_rbh}}$ |
|---|---|---|---|---|---|---|---|---|---|
| **0.0000** | OLS (Independent) | 526.154 | -1036.31 | -1008.56 | 0.99787 | +0.00193 (0.00038) | $9.88 \times 10^{-7}$ | -0.10541 (0.02398) | $1.67 \times 10^{-5}$ |
| **0.2500** | Fixed $\lambda$ | 595.697 | -1175.39 | -1147.65 | 0.99973 | +0.00198 (0.00034) | $4.59 \times 10^{-8}$ | -0.13175 (0.02172) | $8.77 \times 10^{-9}$ |
| **0.5000** | Fixed $\lambda$ | 617.167 | -1218.33 | -1190.59 | 0.99961 | +0.00193 (0.00035) | $2.37 \times 10^{-7}$ | -0.13776 (0.02077) | $3.29 \times 10^{-10}$ |
| **0.7500** | Fixed $\lambda$ | 633.846 | -1251.69 | -1223.95 | 0.99942 | +0.00185 (0.00037) | $1.49 \times 10^{-6}$ | -0.14036 (0.01986) | $2.37 \times 10^{-11}$ |
| **0.9158** | **Profile MLE ($\hat{\lambda}_{\mathrm{ML}}$)** | **640.933** | **-1265.87** | **-1238.12** | **0.99916** | **+0.00170 (0.00040)** | **$3.12 \times 10^{-5}$** | **-0.13748 (0.01917)** | **$7.50 \times 10^{-13}$** |
| **1.0000** | Strict Brownian Motion | 608.673 | -1201.35 | -1173.60 | 0.99875 | +0.00148 (0.00053) | $0.00565$ | -0.12356 (0.01986) | $2.84 \times 10^{-9}$ |

*Observation:* Both $\beta_{\mathrm{AAI}}$ and $\beta_{\mathrm{frac\_rbh}}$ are remarkably stable across the entire spectrum of phylogenetic autocorrelation from $\lambda = 0.25$ to $\lambda = 1.00$. The profile MLE $\hat{\lambda} = 0.9158$ demonstrates strong phylogenetic inertia, improving model AIC by 229.56 units over non-phylogenetic OLS.

---

### Table 3: Robustness Across Lineage and Quality Sensitivities

| Analysis Cohort | $n$ | Profile $\hat{\lambda}_{\mathrm{ML}}$ | Model $R^2$ | $\beta_{\mathrm{AAI}}$ (SE) | $p_{\mathrm{AAI}}$ | $\beta_{\mathrm{frac\_rbh}}$ (SE) | $p_{\mathrm{frac\_rbh}}$ | $\beta_{\log_{10}\mathrm{len}}$ (SE) | $p_{\log_{10}\mathrm{len}}$ | $\beta_{\log_{10}\mathrm{N50}}$ (SE) | $p_{\log_{10}\mathrm{N50}}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Primary (All 237 Reps)** | 237 | 0.9158 | 0.9992 | +0.00170 (0.00040) | $3.12 \times 10^{-5}$ | -0.13748 (0.01917) | $7.50 \times 10^{-13}$ | +0.11543 (0.01812) | $1.01 \times 10^{-9}$ | -0.00900 (0.00299) | $0.00285$ |
| **Restricted: `c__SW-10-69-26` Only** | 185 | 0.9008 | 0.9988 | +0.00192 (0.00048) | $0.00010$ | -0.18044 (0.02219) | $1.34 \times 10^{-13}$ | +0.10202 (0.01918) | $3.93 \times 10^{-7}$ | -0.00200 (0.00371) | $0.59076$ |
| **Dominant Family Omitted (No `JAJPHT01`)** | 97 | 0.9011 | 0.9988 | +0.00149 (0.00049) | $0.00326$ | -0.10238 (0.02942) | $0.00078$ | +0.10431 (0.03666) | $0.00566$ | -0.01355 (0.00427) | $0.00213$ |
| **Caller Quality Proxy (`frac_complete`)** | 232 | 0.9135 | 0.9992 | +0.00155 (0.00044) | $0.00059$ | -0.10476 (0.02642) | $0.00010$ | +0.13324 (0.01955) | $2.32 \times 10^{-10}$ | -0.00158 (0.00475) | $0.74031$ |

*Sensitivity Inferences:*
1. **Lineage Restriction:** In `c__SW-10-69-26` alone ($n = 185$), the results hold firmly ($\beta_{\mathrm{AAI}} = +0.00192, p = 0.00010; \beta_{\mathrm{frac\_rbh}} = -0.18044, p = 1.34 \times 10^{-13}$). The findings are not an artifact of pooling `SW-10-69-26` with `E2`.
2. **Lineage Replication:** Removing the dominant family `f__JAJPHT01` ($n = 97$ remaining across 6 families) preserves both significant associations ($\beta_{\mathrm{AAI}} = +0.00149, p = 0.00326; \beta_{\mathrm{frac\_rbh}} = -0.10238, p = 0.00078$).
3. **Gene-Caller Sensitivity:** Controlling for complete gene call rate (`frac_complete`, $n = 232$, 5 MAGs excluded due to unparsed CheckM2 FAA calls) does not alter the sign or significance of AAI or shared RBH coverage.

---

### Table 4: Audited Reconciled Coupling & Allometric Scaling in `c__SW-10-69-26` ($n = 185$)

| Metric / Comparison | Variable Pair | Value | $p$-value | Interpretation & Meaning |
|---|---|---|---|---|
| **Argon Ratio Correlation** | $\mathrm{dark\_fraction}$ vs $\mathrm{aai\_frac\_rbh}$ | $\rho = -0.30102$ | $3.13 \times 10^{-5}$ | Negative, but used distorted $\min(N_a, N_b)$ denominator. |
| **Corrected Ratio Correlation** | $\mathrm{dark\_fraction}$ vs $\mathrm{frac\_rbh\_desert}$ | $\rho = -0.39418$ | $2.83 \times 10^{-8}$ | Negative, true desert MAG proteome denominator. |
| **Raw Decoupled Counts** | $N_{\mathrm{dark}}$ vs $N_{\mathrm{RBH}}$ | $\rho = +0.24189$ | $9.09 \times 10^{-4}$ | **Positive**: larger genomes harbor more dark and more RBH genes. |
| **Gene Scaling with Proteome** | $N_{\mathrm{dark}}$ vs $N_{\mathrm{prot}}$ | $\rho = +0.94574$ | $2.82 \times 10^{-91}$ | **Near-linear scaling**: dark genes drive 95% of total proteome variation. |
| **Gene Scaling with Proteome** | $N_{\mathrm{RBH}}$ vs $N_{\mathrm{prot}}$ | $\rho = +0.32941$ | $4.68 \times 10^{-6}$ | **Weak scaling**: conserved RBH genes expand far more slowly. |
| **Independent-Count Null** | Permuted $N_{\mathrm{dark}}, N_{\mathrm{RBH}}$ over $N_{\mathrm{prot}}$ | Null mean $\rho = +0.30174$ | $9.999 \times 10^{-5}$ | 95% CI: $[+0.1789, +0.4194]$. Common denominator alone creates positive correlation. Observed negative correlation departs strongly due to bounded sum ($N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$). |

---

### Table 5: Descriptive Accessory Dark Family Sharing Distribution

Mapping of 124,022 contigs across all 237 representatives to 14,645 unique admitted dark families:

| Sharing Tier | Definition | Number of Families | Percentage of Families | Biological Meaning |
|---|---|---|---|---|
| **Singletons** | Present in exactly 1 genome | 3,594 | 24.54% | Lineage-specific genes, rare accessory content, or solitary unclustered calls. |
| **Rare Accessory** | Present in 2 to 11 genomes ($< 5\%$) | 10,312 | 70.41% | Low-frequency accessory pool shared among close relatives. |
| **Intermediate Accessory** | Present in 12 to 47 genomes (5% to 20%) | 727 | 4.96% | Moderately shared lineage-specific families. |
| **Frequent Accessory** | Present in 48 to 118 genomes (20% to 50%) | 12 | 0.08% | Widely shared clade-specific dark families. |
| **Core-like** | Present in $\ge 119$ genomes ($\ge 50\%$) | 0 | 0.00% | No dark family achieves 50% phylum-wide occupancy. |

- **Between-Class Family Richness:**
  * `c__SW-10-69-26` ($n = 185$): median 244.0 distinct dark families per representative (IQR 153.0–363.0).
  * `c__E2` ($n = 51$): median 77.0 distinct dark families per representative (IQR 53.0–115.0).
  * Mann-Whitney test: $U = 1,180.5, p = 1.22 \times 10^{-16}$.
- **Methodological Boundary:** This documents that `SW-10-69-26` possesses higher accessory family richness, but does **not** prove evolutionary expansion vs reduction without outgroup-rooted phylogenetic character mapping.

---

## 6. Informative Negative Outcomes and Scientific Interpretation

1. **Rejection of Reference Proximity Dilution Within Phylum:**  
   The naive expectation that dark fraction tracks sequence similarity to known isolate genomes is **not supported** within the desert *Thermoplasmatota* cohort. In PGLS, the AAI slope is slightly positive ($\beta = +0.00170, p = 3.12 \times 10^{-5}$), meaning that MAGs with higher AAI to reference genomes do not display lower darkness. Reference bias acts as a major ascertainment filter at the root of the phylum (all 237 genomes are phylogenetically distant from known isolates, AAI 55–72%), but within this desert radiation, fine-scale AAI variation does not explain dark fraction.
2. **Quality & Assembly Fragmentation Dominate Residual Variation:**  
   Assembly length ($\beta = +0.11543, p = 1.01 \times 10^{-9}$) and contig fragmentation ($\log_{10}\mathrm{N50}$, $\beta = -0.00900, p = 0.00285$) are strong, persistent predictors of dark fraction. Genomes assembled into smaller contigs have higher dark fractions, partly due to fragmented partial gene calls.
3. **Decoupled Allometric Partitioning:**  
   The negative correlation between dark fraction and shared RBH fraction ($\rho = -0.3942$) is neither a pure mathematical artifact of common-denominator division (which would produce a positive correlation of $+0.3017$), nor evidence of biological displacement. Rather, it reflects **allometric proteome partitioning**: as genome size expands, newly acquired or retained genes are overwhelmingly lineage-specific dark genes ($\rho = 0.95$ with proteome size), which mechanically dilutes the fraction of the proteome composed of slow-expanding conserved RBH genes ($\rho = 0.33$).
4. **No Unqualified Evolutionary Polarity:**  
   While `c__SW-10-69-26` genomes are larger and contain more distinct dark families (median 244 vs 77, $p = 1.22 \times 10^{-16}$), this cannot be asserted as "accessory expansion" in `SW-10-69-26` without ruling out ancestral genome streamlining in `E2`.

---

## 7. Full-Test Recipe and Resource Requirements

If expanded into a formal publication module, the complete confirmatory analysis requires:
1. **Ancestral Proteome Reconstruction:** Fit maximum-likelihood continuous-trait models (Brownian motion vs Ornstein-Uhlenbeck) and discrete gain/loss models (Dollo parsimony / Count) across the full 1,914-taxon GTDB r232 ar53 tree to determine the ancestral state and establish whether `SW-10-69-26` underwent gene gain or `E2` underwent gene reduction.
2. **Reciprocal Codon-Level Gene-Caller Reconciliation:** Benchmark gene calls on raw contigs using Prodigal, Pyrodigal, and Prodigal-gv, mapping exact codon coordinates to isolate genuine complete open reading frames from fragmentation-induced partial calls.
3. **Environmental Covariance:** Regulate for depth, salinity, and geographic sub-basin within the Rub' al-Khali to test whether accessory dark content correlates with ecological niche differentiation.
4. **Compute Resources:** 16 CPUs, 64 GB RAM, walltime ~2.5 hours on IBEX `batch` partition.

---

## 8. Go/No-Go Decision Boundaries

- **GO Condition:** Retain within-*Thermoplasmatota* findings as a rigorous, bounded case study demonstrating:
  1. Strong phylogenetic inertia ($\hat{\lambda} = 0.9158$).
  2. Allometric proteome expansion dominated by lineage-specific dark families.
  3. The absence of fine-scale reference-distance dilution within the phylum.
- **NO-GO Condition:** **Do NOT claim:**
  - That external reference divergence is "refuted" across the tree of life.
  - That the negative ratio correlation is "100% mathematical artifact" (it departs from the positive null due to finite proteome constraints).
  - That `SW-10-69-26` underwent "accessory expansion" as a proven evolutionary polarity.
  - That singleton families are mere "gene-caller fragmentation artifacts".
