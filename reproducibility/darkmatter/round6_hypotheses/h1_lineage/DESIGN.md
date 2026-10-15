# DESIGN (v3): Hypothesis 1 — Within-Thermoplasmatota Dark Proteome Associations

**Version:** v3 (Pre-Execution Design Specification; supersedes `DESIGN_v2_superseded.md`)  
**Date:** 2026-10-13  
**Scope:** Bounded Finite-Cohort, Phylogeny-Conditional Association Study ($n = 237$ Desert *Thermoplasmatota* Species Representatives)

---

## 1. Scientific Objective and Bounded Scope

Within the Rub' al-Khali desert metagenome-assembled genome (MAG) catalogue ($5,438$ refined species representatives), the archaeal phylum *Thermoplasmatota* ($n = 237$ species representatives) exhibits an elevated dark proteome fraction (median $30.96\%$). In GTDB-Tk r232 (`ar53`), these $237$ representatives span three classes:
- `c__SW-10-69-26`: $n = 185$ representatives ($4$ named families across $2$ orders, plus $7$ representatives with unassigned family; $140 / 185$ belong to `f__JAJPHT01`).
- `c__E2`: $n = 51$ representatives ($3$ named families in $1$ order `o__JACPAO01`, plus $2$ representatives with unassigned family).
- `c__Thermoplasmata`: $n = 1$ representative (`f__Methanomassiliicoccaceae`).

### Bounded Research Question
Across the $n = 237$ desert *Thermoplasmatota* species representatives conditioned on continuous branch-length phylogenetic covariance from the GTDB-Tk r232 `ar53` marker tree (`thermo_tree.nwk`) and assembly quality covariates, how does the per-genome dark proteome fraction ($Y_i = N_{\mathrm{dark}, i} / N_{\mathrm{prot, cat}, i}$) vary with:
1. Continuous amino acid identity (`aai`) to the nearest external GTDB r232 or MGnify reference genome?
2. Desert-normalized shared reciprocal-best-hit (RBH) gene coverage (`frac_rbh_desert`)?

### Explicit Epistemic and Methodological Boundaries
1. **Finite-Cohort Observational Scope:** Inferences apply strictly to conditional associations within the observed $237$ desert *Thermoplasmatota* representatives and their placement in the $1,914$-taxon `ar53` phylogeny. Neither positive, negative, nor null slopes establish universal laws across the tree of life.
2. **Two-Sided Exploratory Testing:** All primary and sensitivity regression tests are two-sided ($H_0: \beta_j = 0$ vs $H_1: \beta_j \ne 0$). A negative, null, or positive conditional slope is equally informative for falsifying or bounding simple reference-proximity dilution within this cohort. No post-hoc one-sided confirmatory hypothesis is retrofitted after inspecting parameter estimates.
3. **External Reference Identity (Not Assumed Isolates):** The comparator database constructed in `v2_thermo_relatives.py` comprises $218$ GTDB r232 species representatives and $29$ distinct non-duplicate MGnify genomes in `c__SW-10-69-26` and `c__E2` (after excluding $54$ MGnify assemblies that duplicate a GTDB r232 representative). Across the $237$ desert representatives, the nearest comparator ($\ge 100$ RBHs) is a GTDB r232 representative for $229$ genomes and an MGnify genome for $8$ genomes ($32$ unique comparator accessions). These comparators are **external GTDB r232 and MGnify reference genomes (predominantly uncultured environmental MAGs)**, not cultivated isolates: $145 / 237$ nearest pairs carry MGnify environmental biome metadata ($142$ `Environmental:Terrestrial:Soil`, $3$ `Environmental:Aquatic:Marine:Sediment`) and $92 / 237$ have GTDB-only provenance without cultivation metadata.
4. **Empirical AAI Domain:** Across all $237$ representatives, nearest-reference AAI spans **$48.678\%$ to $80.041\%$** (median **$64.774\%$**, IQR **$58.763\%\text{–}66.392\%$**). Across the $236$ `c__SW-10-69-26` and `c__E2` representatives (the two classes represented in the comparator pool), AAI spans **$50.389\%$ to $80.041\%$** (median **$64.781\%$**; `c__SW-10-69-26` spans $50.389\%\text{–}80.041\%$, `c__E2` spans $58.763\%\text{–}76.964\%$); the single `c__Thermoplasmata` representative has a between-class AAI of **$48.678\%$**.
5. **No Macroevolutionary Polarity or Causal Gene-Caller Claims:** Descriptive differences in genome size or dark family counts between `c__SW-10-69-26` and `c__E2` cannot distinguish evolutionary gene gain from gene reduction without outgroup-rooted ancestral state reconstruction, and covariates measuring contig fragmentation (`log10(N50)`, CheckM2 `frac_complete`) measure conditional quality associations rather than proving that unclustered genes are caller artifacts.

---

## 2. Primary Identifiable Estimand and PGLS Model Specification

### Model Formulation
We specify a single primary Phylogenetic Generalized Least Squares (PGLS) linear model across the $n = 237$ species representatives:

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim \mathcal{N}\left(\mathbf{0}, \sigma^2 \mathbf{V}(\lambda)\right)$$

where:
- $\mathbf{y} \in \mathbb{R}^{237}$ is the catalogue dark proteome fraction $Y_i = N_{\mathrm{dark}, i} / N_{\mathrm{prot, cat}, i}$.
- $\mathbf{V}(\lambda) = \lambda \mathbf{C} + (1 - \lambda)\mathbf{I}_{237}$ for Pagel's $\lambda \in [0, 1]$, where $\mathbf{C} \in \mathbb{R}^{237 \times 237}$ is the phylogenetic correlation matrix derived from the $1,914$-taxon `ar53` FastTree (`thermo_tree.nwk`).
- $\mathbf{X} \in \mathbb{R}^{237 \times 7}$ (full column rank $p = 7$) comprises an intercept and $6$ pre-specified covariates:
  1. `completeness`: CheckM2 completeness ($\%$).
  2. `contamination`: CheckM2 contamination ($\%$).
  3. `log_len`: $\log_{10}(\mathrm{assembly\ length\ in\ bp})$.
  4. `log_n50`: $\log_{10}(\mathrm{contig\ N50\ in\ bp})$.
  5. `aai`: RBH amino acid identity ($\%$) to the nearest external GTDB r232 / MGnify reference genome.
  6. `frac_rbh_desert`: shared RBH protein count divided by the desert MAG's single-genome Pyrodigal protein count ($N_{\mathrm{RBH}, i} / N_{\mathrm{prot, aai}, i}$).

### Primary Estimands and Two-Sided Hypotheses
1. **Conditional External Proximity Slope ($\beta_{\mathrm{AAI}}$):** Expected change in dark fraction per $1$ percentage-point increase in nearest-reference AAI, holding `frac_rbh_desert`, assembly quality, and phylogenetic covariance constant ($H_0: \beta_{\mathrm{AAI}} = 0$ vs $H_1: \beta_{\mathrm{AAI}} \ne 0$, two-sided).
2. **Conditional Shared Coverage Slope ($\beta_{\mathrm{frac\_rbh\_desert}}$):** Expected change in dark fraction per unit increase in desert-normalized RBH coverage, holding `aai`, assembly quality, and phylogenetic covariance constant ($H_0: \beta_{\mathrm{frac\_rbh\_desert}} = 0$ vs $H_1: \beta_{\mathrm{frac\_rbh\_desert}} \ne 0$, two-sided).

### Separation of ML Profile Likelihood and Unbiased Fixed-Covariance Inference
For any candidate $\lambda \in [0, 1]$, let $\mathbf{V}(\lambda) = \mathbf{L}\mathbf{L}^T$ be the lower-triangular Cholesky factorization, with whitened 1D response $\tilde{\mathbf{y}} = \mathbf{L}^{-1}\mathbf{y} \in \mathbb{R}^n$ and whitened design matrix $\tilde{\mathbf{X}} = \mathbf{L}^{-1}\mathbf{X} \in \mathbb{R}^{n \times p}$:
1. **GLS Point Estimate and Residual Sum of Squares:**
   $$\hat{\boldsymbol{\beta}}(\lambda) = (\tilde{\mathbf{X}}^T \tilde{\mathbf{X}})^{-1} \tilde{\mathbf{X}}^T \tilde{\mathbf{y}}, \quad \mathrm{RSS}(\lambda) = \|\tilde{\mathbf{y}} - \tilde{\mathbf{X}}\hat{\boldsymbol{\beta}}(\lambda)\|_2^2$$
2. **Concentrated ML Variance and Profile Log-Likelihood (for $\lambda$ estimation only):**
   $$\hat{\sigma}^2_{\mathrm{ML}}(\lambda) = \frac{\mathrm{RSS}(\lambda)}{n}, \quad \ell_{\mathrm{profile}}(\lambda) = -\frac{n}{2}\ln(2\pi) - \frac{n}{2}\ln\hat{\sigma}^2_{\mathrm{ML}}(\lambda) - \frac{1}{2}\ln|\mathbf{V}(\lambda)| - \frac{n}{2}$$
   Pagel's $\hat{\lambda}_{\mathrm{ML}}$ maximizes $\ell_{\mathrm{profile}}(\lambda)$ over $\lambda \in [0, 1]$. The $95\%$ profile-likelihood confidence interval for $\lambda$ is:
   $$\mathrm{CI}_{95\%}(\lambda) = \left\{\lambda \in [0, 1] : 2\left[\ell_{\mathrm{profile}}(\hat{\lambda}_{\mathrm{ML}}) - \ell_{\mathrm{profile}}(\lambda)\right] \le \chi^2_{1, 0.95} = 3.8414588\right\}$$
3. **Unbiased Fixed-Covariance Variance, Standard Errors, and Student-$t$ Intervals:**
   Conditional on a fixed $\lambda$ (including $\lambda = \hat{\lambda}_{\mathrm{ML}}$), unbiased variance and covariance of $\hat{\boldsymbol{\beta}}(\lambda)$ are:
   $$\hat{\sigma}^2_{\mathrm{unbiased}}(\lambda) = \frac{\mathrm{RSS}(\lambda)}{n - p}, \quad \widehat{\mathrm{Cov}}(\hat{\boldsymbol{\beta}}(\lambda)) = \hat{\sigma}^2_{\mathrm{unbiased}}(\lambda) (\tilde{\mathbf{X}}^T \tilde{\mathbf{X}})^{-1}$$
   Standard errors are $\mathrm{SE}(\hat{\beta}_j) = \sqrt{[\widehat{\mathrm{Cov}}(\hat{\boldsymbol{\beta}})]_{jj}}$, two-sided Wald $t$-statistics are $t_j = \hat{\beta}_j / \mathrm{SE}(\hat{\beta}_j)$ evaluated against a Student-$t$ distribution with $df = n - p = 230$, and $95\%$ confidence intervals conditional on $\lambda$ are:
   $$\hat{\beta}_j \pm t_{0.975, n - p}\,\mathrm{SE}(\hat{\beta}_j)$$
   Because fixed-$\hat{\lambda}_{\mathrm{ML}}$ intervals do not analytically propagate uncertainty in $\hat{\lambda}_{\mathrm{ML}}$, we explicitly label them **conditional on $\hat{\lambda}_{\mathrm{ML}}$** and report full parameter, SE, CI, and $p$-value sensitivity across both $\lambda \in \{0.0, 0.25, 0.50, 0.75, 1.0\}$ and the $95\%$ profile-likelihood bounds $[\lambda_{\mathrm{CI,low}}, \lambda_{\mathrm{CI,high}}]$.
4. **Correct 1D Whitened GLS $R^2(\lambda)$:**
   Let $\tilde{\mathbf{1}} = \mathbf{L}^{-1}\mathbf{1}_n \in \mathbb{R}^n$ (strictly 1D vector of length $n$) and $\hat{\mu}_{\mathrm{GLS}}(\lambda) = (\tilde{\mathbf{1}}^T \tilde{\mathbf{y}}) / (\tilde{\mathbf{1}}^T \tilde{\mathbf{1}})$. The whitened intercept-only total sum of squares is $\mathrm{TSS}_{\mathrm{GLS}}(\lambda) = \sum_{i=1}^n (\tilde{y}_i - \tilde{1}_i \hat{\mu}_{\mathrm{GLS}}(\lambda))^2$, and:
   $$R^2_{\mathrm{GLS}}(\lambda) = 1 - \frac{\mathrm{RSS}(\lambda)}{\mathrm{TSS}_{\mathrm{GLS}}(\lambda)}$$
   Both `y_tilde` and `ones_tilde` are asserted to be 1D arrays of shape `(n,)` so no $(n, n)$ broadcasting can occur. Note that $R^2_{\mathrm{GLS}}(\lambda)$ measures goodness of fit within the fixed $\mathbf{V}(\lambda)^{-1/2}$ metric at that $\lambda$; across different $\lambda$ values, model fit is compared via $\ell_{\mathrm{profile}}(\lambda)$ and information criteria (where `fit_pgls` tabulates fixed-$\lambda$ AIC/BIC with $k = p + 1$ parameters counting $\boldsymbol{\beta}$ and $\sigma^2$, and full-ML comparison at profiled $\hat{\lambda}_{\mathrm{ML}}$ counts $k = p + 2$).
5. **Independent Trusted Implementation Cross-Check:**
   Every fixed-$\lambda$ GLS fit ($\lambda = 0.0$, $\lambda = 0.5$, and $\lambda = \hat{\lambda}_{\mathrm{ML}}$) is verified programmatically against `statsmodels.regression.linear_model.GLS(y, X, sigma=V)` (asserting maximum absolute difference $< 10^{-10}$ in $\hat{\boldsymbol{\beta}}$, $\mathrm{SE}_{\mathrm{unbiased}}$, $t$, $p$, $\ell$, $\hat{\sigma}^2_{\mathrm{unbiased}}$, and $R^2_{\mathrm{GLS}}$).

---

## 3. Phylogenetic Covariance Assumptions and Sensitivities

1. **Tree Provenance and Rooting Validation:**
   `thermo_tree.nwk` is an approximate-ML FastTree (`LG+Gamma`) built from the GTDB-Tk r232 `ar53` alignment across $1,914$ taxa: $237$ desert *Thermoplasmatota* representatives, $1,647$ GTDB r232 *Thermoplasmatota* species representatives, and $30$ archaeal outgroup representatives (`p__Halobacteriota` and `p__Thermoproteota`, listed in `thermo_tree_taxa.tsv`). We verify programmatically that midpoint rooting of the $1,914$-tip tree induces the exact same crown split on the pruned $237$-tip desert tree (`c__E2` vs `c__SW-10-69-26 + c__Thermoplasmata`) and the exact same cophenetic/shared-branch VCV structure as explicit outgroup rooting on the $30$ archaeal outgroups.
2. **Correlation Normalization vs Unnormalized Brownian VCV:**
   Because the FastTree branch lengths are in substitutions per site and non-ultrametric (root-to-tip distances on the pruned $237$-tip tree span $0.2587$ to $0.5283$), the primary matrix $\mathbf{C}$ normalizes raw shared branch lengths $V_{ij}$ to unit diagonal: $C_{ij} = V_{ij} / \sqrt{V_{ii} V_{jj}}$. To verify that unit-diagonal normalization does not drive the results, we run an explicit sensitivity PGLS using the **unnormalized Brownian VCV matrix** scaled by its mean diagonal ($\mathbf{V}_{\mathrm{raw}} / \overline{\mathrm{diag}(\mathbf{V}_{\mathrm{raw}})}$).
3. **Pre-Specified Sensitivity Cohorts and Specifications:**
   - **S1 (Primary Cohort, $n = 237$):** All $237$ desert representatives with correlation matrix $\mathbf{C}$ and `frac_rbh_desert = aai_n_rbh / n_prot_aai`.
   - **S2 (Unnormalized Brownian VCV, $n = 237$):** Uses mean-diagonal-scaled $\mathbf{V}_{\mathrm{raw}}$ instead of correlation matrix $\mathbf{C}$.
   - **S3 (Catalogue-Denominator RBH Fraction, $n = 237$):** Replaces `frac_rbh_of_desert_aai` with `frac_rbh_of_desert_cat = aai_n_rbh / n_proteins_cat` so both response and RBH covariate share the exact catalogue denominator.
   - **S4 (`c__SW-10-69-26` Only, $n = 185$):** Restricts to the dominant class to eliminate between-class leverage.
   - **S5 (Dominant Family `f__JAJPHT01` Omitted, $n = 97$):** Omits the $140$ `f__JAJPHT01` genomes to verify persistence across the remaining clades.
   - **S6 (`c__SW-10-69-26` + `c__E2` Same-Class Comparator Cohort, $n = 236$):** Omits the single `c__Thermoplasmata` genome whose nearest AAI is between-class.
   - **S7 (CheckM2 Caller Complete-Gene Proxy `frac_complete`, $n = 232$):** Adds `frac_complete` (`partial=00` fraction from `/ibex/scratch/projects/c2014/rob/metagenomics/05_refine/_checkm2_refined/protein_files/{genome}.faa`). Exactly $5$ MAGs (`48PRr1__SemiBin_112`, `48PRr1__SemiBin_141`, `48PRr1__SemiBin_510`, `48PRr1__SemiBin_539`, `48PRr1__SemiBin_849`) lack `.faa` files in `_checkm2_refined/protein_files/` and are documented as explicit cohort attrition ($237 \to 232$).

---

## 4. Caller Basis, Non-Disjointness of $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$, and Algebraic Coupling Audit

1. **Distinct Gene-Caller Bases and Non-Exclusivity:**
   - `n_dark` and `n_proteins_cat` (`v2_reps.tsv`) are computed on the metagenomic catalogue gene calls mapped to $90\%$ identity / $80\%$ coverage (`all_cluster90.tsv`) cluster representatives lacking `UniRef50 2026_01` and functional annotations (`dark_rep_ids.txt.gz`).
   - `aai_n_rbh` and `n_prot_aai` (`v2_thermo_relatives_aai.tsv.gz`) are computed on independent single-genome `Pyrodigal 3.7.1 (meta=False)` predictions aligned via `DIAMOND 2.1.16 blastp --very-sensitive` ($e \le 10^{-5}$, identity $\ge 30\%$, query and subject coverage $\ge 70\%$) against the $218$ GTDB r232 and $29$ MGnify comparator genomes.
   - Consequently, $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ come from **different gene callers** (`n_prot_aai - n_proteins_cat` ranges from $-13$ to $+93$ proteins per genome, median $+18$, Pearson $r = 0.999752$) and **non-mutually-exclusive criteria**: a protein whose cluster representative lacks a hit in `UniRef50 2026_01` can still form a reciprocal best hit with an unannotated protein in a recently assembled environmental MAG in GTDB r232 or MGnify. Because `v2_thermo_relatives_aai.py` retained only pair-level RBH counts (`v2_thermo_relatives_aai.tsv.gz`), gene-level intersection is unobserved, and we make **no claim** that $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ form a disjoint partition ($N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$).
2. **Matched-Denominator Algebraic Permutation Illustration:**
   - To illustrate the purely algebraic effect of dividing two independent positive count vectors by a shared variable denominator $D_i$, we permute $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ independently across genomes ($10,000$ iterations, seed $42$, plus-one empirical tail $(k+1)/(N+1)$) and divide **both** permuted numerators by the **same** denominator $D_i$, comparing against the observed correlation of $(N_{\mathrm{dark}, i}/D_i, N_{\mathrm{RBH}, i}/D_i)$ computed with that **exact same denominator** $D_i$ (evaluated both for $D_i = N_{\mathrm{prot, cat}, i}$ and for native denominators).
   - This permutation is strictly an **algebraic illustration** showing that a shared variable denominator induces a positive correlation ($\bar{\rho} \approx +0.302$) between independent numerators, **not a biological null model**. The observed negative correlation between dark fraction and RBH fraction ($\rho = -0.3918$ under shared $N_{\mathrm{prot, cat}}$; $\rho = -0.3942$ under native denominators) coexists with a positive raw-count correlation ($\rho(N_{\mathrm{dark}}, N_{\mathrm{RBH}}) = +0.2419$) and reflects the descriptive association that $N_{\mathrm{dark}}$ scales more strongly with total protein count ($\rho = +0.9457$) than $N_{\mathrm{RBH}}$ ($\rho = +0.3294$).

---

## 5. Descriptive Accessory Dark Family Sharing & Future Work Outside Scope

- Using `contig_to_bin_refined.tsv.gz` and `member_fam_final.tsv`, we summarize the prevalence distribution of the $14,645$ admitted dark families across the $237$ representatives (singletons in $1$ genome, rare $2\text{–}11$, intermediate $12\text{–}47$, frequent $48\text{–}118$, core $\ge 119$) and compare per-genome family counts between `c__SW-10-69-26` and `c__E2` descriptively.
- Outgroup-rooted ancestral state reconstruction (e.g., Dollo parsimony or birth-death gain/loss models across the $1,914$-taxon tree) and raw-contig multi-caller re-prediction are **optional future extensions outside the scope of this bounded study**.
