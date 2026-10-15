# DOSSIER: H1 — Phylogeny-Conditional Proteome Partitioning and Sequence Divergence in Desert *Thermoplasmatota* (v3 — Round 6 Repaired Implementation)

## 1. Executive Summary & Bounded Scientific Claim

Across the **237 desert *Thermoplasmatota* species-cluster representatives** recovered in the Empty Quarter metagenomic assembly catalogue (`185` `c__SW-10-69-26`, `51` `c__E2`, and `1` `c__Thermoplasmata`), we executed a pre-registered, two-sided exploratory phylogenetic generalized least squares (PGLS) and denominator-reconciled association study (`DESIGN.md` v3; IBEX Slurm Job `53704703`). This v3 implementation resolves all five critical issues (`C1–C5`) and optional audit items (`O1–O2`) identified in the Round 6 independent review (`reviews/darkmatter-2026-10-13/round6/h1_review1/review.txt`) and root audit (`reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`, `reviews/darkmatter-2026-10-13/round6/root_checks/h1_v2_gls_audit.json`). Superseded v1 and v2 scripts, tables, and documentation (`DESIGN_v2_superseded.md`, `DOSSIER_v2_superseded.md`, `INPUTS_v2_superseded.tsv`) are preserved intact for auditability.

### Bounded Finite-Cohort Finding
Within this finite cohort of 237 desert *Thermoplasmatota* species representatives, after accounting for shared phylogenetic ancestry on the 1,914-tip GTDB-Tk `ar53` archaeal reference phylogeny (`thermo_tree.nwk`, profiled Pagel's $\hat{\lambda}_{\mathrm{ML}} = 0.915784$, 95% profile-likelihood CI $[0.846303, 0.957823]$, likelihood-ratio test vs $\lambda = 0$: $2\Delta\ell = 229.5585$, $p = 7.4411 \times 10^{-52}$) and controlling for assembly completeness, contamination, $\log_{10}(\text{genome length})$, and $\log_{10}(\text{N50})$, genomic dark proteome fraction (`dark_fraction` $= N_{\mathrm{dark}} / N_{\mathrm{proteins,cat}}$) exhibits two distinct, statistically significant two-sided partial associations (`results/v3_pgls_model_summary.tsv`, `results/v3_pgls_ancestry_summary.txt`):
1. **Conserved gene-content overlap (`frac_rbh_desert` $= N_{\mathrm{RBH}} / N_{\mathrm{prot,aai}}$):** Strong negative partial slope ($\hat{\beta}_{\mathrm{RBH}} = -0.137484$, unbiased $\mathrm{SE}_{\mathrm{unb}} = 0.019456$, $t_{230} = -7.0664$, two-sided $p = 1.879206 \times 10^{-11}$, conditional 95% Student-$t$ CI $[-0.175819, -0.099149]$).
2. **Amino-acid identity (`aai`) across conserved bidirectional best hits to the nearest external GTDB r232 or MGnify reference genome:** Positive partial slope ($\hat{\beta}_{\mathrm{AAI}} = +0.001704$ per percentage-point AAI, unbiased $\mathrm{SE}_{\mathrm{unb}} = 0.000407$, $t_{230} = +4.1847$, two-sided $p = 4.065008 \times 10^{-5}$, conditional 95% Student-$t$ CI $[+0.000902, +0.002506]$).

The repaired 1D whitened fixed-covariance GLS coefficient of determination at $\hat{\lambda}_{\mathrm{ML}} = 0.915784$ is **$R^2_{\mathrm{GLS}} = 0.469051$** ($\mathrm{RSS} = 0.243540$, $\mathrm{TSS}_{\mathrm{GLS}} = 0.458687$), and for non-phylogenetic OLS ($\lambda = 0$) is **$R^2 = 0.495786$**, matching `statsmodels.GLS` to machine precision ($\max|\Delta| < 5.2 \times 10^{-12}$; `results/v3_pgls_statsmodels_verification.tsv`). Across all seven pre-registered sensitivity models—including unnormalized Brownian VCV (`S2`), shared catalogue-denominator RBH (`S3`), restriction to `c__SW-10-69-26` (`S4`, $n=185$), omission of the dominant family `f__JAJPHT01` (`S5`, $n=97$), exclusion of the single between-class `c__Thermoplasmata` representative (`S6`, $n=236$), and inclusion of the CheckM2 `frac_complete` caller-fragmentation proxy (`S7`, $n=232$)—and across the entire 95% profile-likelihood confidence interval for $\lambda$ ($[0.846303, 0.957823]$), both partial associations retain the same sign and remain statistically significant (`results/v3_pgls_sensitivity.tsv`, `results/v3_pgls_lambda_profile.tsv`).

### Explicit Scope Boundaries
- **Finite-cohort observational scope:** All inferences are strictly bounded to the 237 sampled desert *Thermoplasmatota* species-cluster representatives and their nearest external GTDB r232 (`n=229`) and MGnify (`n=8`) reference genomes (which are predominantly uncultured environmental MAGs; `142` annotated as `Environmental:Terrestrial:Soil`, `3` as `Environmental:Aquatic:Marine:Sediment`, and `92` GTDB-only genomes lacking habitat strings; `results/v3_comparator_metadata_summary.tsv`).
- **No macroevolutionary or caller-causal extrapolation:** Because this study is cross-sectional and observational, it does not identify evolutionary gain/loss rates, ancestral gene-family polarities, selection regimes, or causal gene-caller mechanics. Future ancestral gene-tree/species-tree reconciliation (e.g., ALE/Generax) is noted solely as optional future work outside this completed finite-cohort experiment.

---

## 2. Summary of Round 6 Critical Repairs (`C1–C5`) and Audit Verifications

| Issue ID | Round 6 Review / Root Audit Defect in v2 | Executed v3 Repair & Empirical Verification | Evidence Files |
| :--- | :--- | :--- | :--- |
| **C1** | **Invalid GLS $R^2$ (`0.9896`–`0.9979`):** `02_pgls_ancestry_models.py` subtracted a `(1,1)` 2D mean from a 1D whitened vector `Ly`, broadcasting `Ly_centered` into a `(n,n)` matrix and inflating `tss` by a factor of $n$. | Flattened all response and predictor arrays to 1D (`y.reshape(-1)`) and asserted `Ly_centered.ndim == 1` before computing `tss_gls`. Verified against `statsmodels.GLS` at $\lambda \in \{0.0, 0.5, 0.915784, 0.915800\}$ (`max_abs_diff_r2 <= 4.44e-15`). Primary PGLS $R^2_{\mathrm{GLS}} = 0.469051$; OLS $R^2 = 0.495786$. | `scripts/02_v3_pgls_ancestry_models.py` (lines 71–150), `results/v3_pgls_model_summary.tsv`, `results/v3_pgls_statsmodels_verification.tsv` |
| **C2** | **ML vs unbiased residual variance, normal quantiles, truncated $p$-values, and $\lambda$ uncertainty:** v2 used $\hat{\sigma}^2_{\mathrm{ML}} = \mathrm{RSS}/n$, normal multiplier `1.96`, truncated `0.000000` $p$-values, and omitted profile-likelihood CI for $\lambda$. | Separated $\hat{\sigma}^2_{\mathrm{ML}} = \mathrm{RSS}/n$ (for log-likelihood profiling) from $\hat{\sigma}^2_{\mathrm{unb}} = \mathrm{RSS}/(n-p)$ (for $\mathrm{SE}_{\mathrm{unb}}$, $t = \hat{\beta}/\mathrm{SE}_{\mathrm{unb}}$, two-sided `scipy.stats.t.sf(|t|, df=n-p)*2` in scientific notation, and exact $t_{0.975, n-p}$ 95% CIs). Added 95% profile-likelihood CI for $\lambda$ (`[0.846303, 0.957823]` in S1) and evaluated fixed-$\lambda$ conditional estimates across the $\lambda$ profile grid. | `scripts/02_v3_pgls_ancestry_models.py` (lines 90–215), `results/v3_pgls_model_summary.tsv`, `results/v3_pgls_lambda_profile.tsv`, `results/v3_pgls_sensitivity.tsv` |
| **C3** | **Post-hoc one-sided directional framing:** v2 framed $\beta_{\mathrm{AAI}} > 0$ as a confirmatory expectation after observing the positive sign in Round 5 (`r = +0.292`). | Rewrote `DESIGN.md` prior to v3 computation to pre-register two-sided exploratory hypotheses ($H_0: \beta_j = 0$ vs $H_1: \beta_j \ne 0$) for both `aai` and `frac_rbh_desert`. | `DESIGN.md` (Section 1), `scripts/02_v3_pgls_ancestry_models.py` |
| **C4** | **Denominator mismatch in permutation (`n_proteins_cat` vs `n_prot_aai`) and false disjoint-set claim:** v2 compared mixed-caller observed ratios against a common-denominator permutation and claimed $N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$. | Computed observed and permuted correlations on **identical denominators** (`n_proteins_cat` for both, and `n_prot_aai` for both). Removed all disjoint-partition claims; documented why $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ come from separate callers and non-exclusive criteria; framed permutation strictly as an algebraic illustration and differential count scaling as descriptive association. | `scripts/03_v3_coupling_and_comparator_audit.py` (lines 79–196), `results/v3_coupling_and_comparator_audit.tsv`, `results/v3_coupling_and_comparator_summary.txt` |
| **C5** | **Comparator mischaracterization ("non-desert isolates") and misquoted AAI range (`50.9%–80.0%`):** Comparators are GTDB r232 (`229`) and MGnify (`8`) genomes (`142` soil MAGs, `3` marine sediment MAGs, `92` unannotated GTDB-only); empirical min AAI is `48.678%` (`50.389%` in `SW-10-69-26 + E2`). | Replaced all "non-desert isolate" language with exact source-verified terminology ("external GTDB r232 and MGnify reference genomes, predominantly uncultured environmental MAGs") and reported exact empirical AAI ranges per stratum (`48.678%–80.041%` across all `237` reps; `50.389%–80.041%` across `236` same-class reps). | `scripts/03_v3_coupling_and_comparator_audit.py` (lines 48–77), `results/v3_comparator_metadata_summary.tsv` |
| **O1 / O2 / Audit** | **Rooting validation, `DYNAMIC` manifest placeholders, and Mann-Whitney $U$ typo:** v2 relied on unverified midpoint rooting, used `DYNAMIC` in `INPUTS.tsv`, and misreported Mann-Whitney $U = 1,180.5$ instead of $8,292.5$ ($U_{\min} = 1,142.5$). | Proved outgroup rooting on the 30 `Halobacteriota`/`Thermoproteota` outgroup tips yields the exact same 237-tip VCV as midpoint rooting (`max_abs_diff = 0.000e+00`); enumerated all 247 raw input entries with SHA-256 checksums (`INPUTS.tsv`); reported exact `scipy.stats.mannwhitneyu` $U = 8292.5$ ($p = 1.224554 \times 10^{-16}$). | `scripts/01_v3_input_manifest_and_phylo_audit.py`, `INPUTS.tsv`, `results/v3_phylo_vcv_validation.json`, `results/v3_closure_verification.json` |

---

## 3. Phylogenetic Covariance Validation and Outgroup Rooting Equivalence

Script `scripts/01_v3_input_manifest_and_phylo_audit.py` (executed with DendroPy `5.0.8` in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python`) audited the 1,914-tip Newick tree (`thermo_tree.nwk`) and taxon metadata (`thermo_tree_taxa.tsv`), producing `results/v3_phylo_vcv_validation.json` and `results/v3_phylo_vcv_validation.txt`:

- **Tip uniqueness and cohort mapping:** `thermo_tree.nwk` contains `1,914` unique leaf labels (`1,648` `GTDB_r232`, `237` `Desert_MAG`, `29` `MGnify`) with zero duplicates. All `237 / 237` (`100%`) desert *Thermoplasmatota* species-cluster representatives map 1-to-1 to tree leaves with zero cohort attrition.
- **Rooting equivalence (midpoint vs explicit outgroup rooting):** `thermo_tree_taxa.tsv` includes `30` explicit archaeal outgroup references (`17` `Halobacteriota` and `13` `Thermoproteota`) alongside `1,884` `Thermoplasmatota` ingroup tips. Because the 237 desert MAGs reside entirely within the `Thermoplasmatota` ingroup clade and the longest path in `thermo_tree.nwk` spans the basal split between the 30 outgroup taxa and `Thermoplasmatota`, rerooting at the MRCA of the 30 `Halobacteriota`/`Thermoproteota` outgroup tips versus midpoint rooting produces **mathematically identical** 237-tip patristic variance-covariance matrices ($\max|\mathbf{V}_{\mathrm{midpoint}} - \mathbf{V}_{\mathrm{outgroup}}| = 0.000000 \times 10^{0}$) and diagonal-normalized correlation matrices ($\max|\mathbf{C}_{\mathrm{midpoint}} - \mathbf{C}_{\mathrm{outgroup}}| = 0.000000 \times 10^{0}$).
- **Symmetry, positive definiteness, and rank:** Both the unnormalized patristic VCV matrix $\mathbf{V}$ (`root-to-tip` diagonal range $[0.258680, 0.528319]$, median $0.408705$, mean $0.407105$, SD $0.043892$ substitutions/site) and the diagonal-normalized correlation matrix $C_{ij} = V_{ij} / \sqrt{V_{ii} V_{jj}}$ (`off-diagonal` range $[0.000000, 0.995025]$) are exactly symmetric ($\max|\mathbf{M} - \mathbf{M}^\top| = 0.0$), full rank ($\mathrm{rank} = 237$), and strictly positive definite ($\lambda_{\min}(\mathbf{V}) = 2.350676 \times 10^{-3} > 0$, $\lambda_{\max}(\mathbf{V}) = 40.950633$; $\lambda_{\min}(\mathbf{C}) = 4.944177 \times 10^{-3} > 0$, $\lambda_{\max}(\mathbf{C}) = 101.468139$, condition number $20522.75$).

---

## 4. Primary Ancestry-Aware PGLS vs Non-Phylogenetic OLS Results (`C1`, `C2`, `C3`)

Script `scripts/02_v3_pgls_ancestry_models.py` (executed in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python` with `numpy 2.2.6`, `scipy 1.15.2`, `pandas 2.3.3`, `statsmodels 0.14.6`) fits the pre-registered linear model across all $n = 237$ desert *Thermoplasmatota* species representatives ($p = 7$ parameters including intercept, residual $df = n - p = 230$, exact two-sided Student-$t$ critical multiplier $t_{0.975, 230} = 1.970332$):

$$\text{dark\_fraction}_i = \beta_0 + \beta_{\mathrm{AAI}}\,\text{aai}_i + \beta_{\mathrm{RBH}}\,\text{frac\_rbh\_desert}_i + \beta_1\,\text{completeness}_i + \beta_2\,\text{contamination}_i + \beta_3\,\log_{10}(\text{genome\_length}_i) + \beta_4\,\log_{10}(\text{n50}_i) + \varepsilon_i$$

where $\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \sigma^2 \boldsymbol{\Sigma}(\lambda))$ and $\boldsymbol{\Sigma}(\lambda) = \lambda \mathbf{C} + (1 - \lambda)\mathbf{I}_n$.

### 4.1 Model-Level Fit Statistics (`S1_Primary_All_237` vs OLS)

| Metric | Primary PGLS ($\hat{\lambda}_{\mathrm{ML}} = 0.915784$) | Non-Phylogenetic OLS ($\lambda = 0.000000$) |
| :--- | :--- | :--- |
| Sample size ($n$), parameters ($p$), residual $df$ ($n-p$) | $n = 237,\; p = 7,\; df = 230$ | $n = 237,\; p = 7,\; df = 230$ |
| Pagel's $\lambda$ MLE (95% profile-likelihood CI) | $0.915784\;\;([0.846303,\; 0.957823])$ | $0.000000$ (fixed) |
| Gaussian log-likelihood ($\ell$) | $+640.932977$ | $+526.153746$ |
| Likelihood-ratio test vs $\lambda = 0$ ($2\Delta\ell$, $p$-value) | $229.5585\;\;(p = 7.4411 \times 10^{-52})$ | Reference ($\lambda = 0$) |
| Fixed-$\lambda$ AIC / BIC ($k = p + 1 = 8$; full profiled-ML $k = p + 2 = 9$ for PGLS: $-1263.865955\;/\;-1232.653413$) | $-1265.865955\;/\;-1238.121473$ | $-1036.307491\;/\;-1008.563010$ |
| Whitened residual sum of squares ($\mathrm{RSS}$) | $0.243540$ | $0.163665$ |
| Whitened total sum of squares ($\mathrm{TSS}_{\mathrm{GLS}}$) | $0.458687$ | $0.324595$ |
| **Fixed-covariance whitened GLS $R^2_{\mathrm{GLS}}$** | **$0.469051$** | **$0.495786$** |
| Unbiased residual variance $\hat{\sigma}^2_{\mathrm{unb}} = \mathrm{RSS}/(n-p)$ | $0.00105887$ | $0.00071159$ |
| Maximum-likelihood variance $\hat{\sigma}^2_{\mathrm{ML}} = \mathrm{RSS}/n$ | $0.00102759$ | $0.00069057$ |

### 4.2 Parameter Estimates, Unbiased Standard Errors, and Conditional 95% Student-$t$ CIs (`results/v3_pgls_model_summary.tsv`)

> **Important conditioning note (`C2`):** The standard errors ($\mathrm{SE}_{\mathrm{unb}}$), $t$-statistics, $p$-values, and 95% Student-$t$ confidence intervals in the table below are **conditional on $\lambda = \hat{\lambda}_{\mathrm{ML}} = 0.915784$ treated as fixed** (using unbiased residual scale $\hat{\sigma}^2_{\mathrm{unb}} = \mathrm{RSS}/(n-p)$ and $df = 230$). They do not propagate uncertainty in estimating $\lambda$; Section 4.3 explicitly evaluates how coefficient estimates and conditional intervals vary across the 95% profile-likelihood confidence interval for $\lambda$.

| Predictor Variable | PGLS $\hat{\beta}$ | $\mathrm{SE}_{\mathrm{unb}}$ ($df=230$) | $\mathrm{SE}_{\mathrm{ML}}$ ($/n$) | $t_{\mathrm{unb}}$ | Two-Sided $p_{\mathrm{unb}}$ | Conditional 95% Student-$t$ CI ($\hat{\lambda}_{\mathrm{ML}}$) | OLS $\hat{\beta}$ ($\lambda=0$) | OLS $\mathrm{SE}_{\mathrm{unb}}$ | OLS $t_{\mathrm{unb}}$ | OLS Two-Sided $p_{\mathrm{unb}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `Intercept` | $-0.453207$ | $0.111867$ | $0.110203$ | $-4.0513$ | $6.967615 \times 10^{-5}$ | $[-0.673623,\; -0.232792]$ | $-1.252281$ | $0.131908$ | $-9.4936$ | $3.003069 \times 10^{-18}$ |
| `completeness` | $-0.000335$ | $0.000126$ | $0.000124$ | $-2.6586$ | $8.398044 \times 10^{-3}$ | $[-0.000583,\; -0.000087]$ | $-0.000506$ | $0.000199$ | $-2.5367$ | $1.185263 \times 10^{-2}$ |
| `contamination` | $-0.000320$ | $0.000491$ | $0.000484$ | $-0.6507$ | $5.159032 \times 10^{-1}$ | $[-0.001288,\; +0.000648]$ | $-0.000106$ | $0.000888$ | $-0.1196$ | $9.048866 \times 10^{-1}$ |
| `log_len` ($\log_{10}\text{bp}$) | $+0.115427$ | $0.018391$ | $0.018117$ | $+6.2764$ | $1.706357 \times 10^{-9}$ | $[+0.079191,\; +0.151662]$ | $+0.249850$ | $0.020731$ | $+12.0522$ | $2.983976 \times 10^{-26}$ |
| `log_n50` ($\log_{10}\text{bp}$) | $-0.009005$ | $0.003031$ | $0.002986$ | $-2.9709$ | $3.284365 \times 10^{-3}$ | $[-0.014976,\; -0.003033]$ | $-0.022638$ | $0.004832$ | $-4.6854$ | $4.782238 \times 10^{-6}$ |
| **`aai` (%)** | **$+0.001704$** | **$0.000407$** | **$0.000401$** | **$+4.1847$** | **$4.065008 \times 10^{-5}$** | **$[+0.000902,\; +0.002506]$** | **$+0.001931$** | **$0.000383$** | **$+5.0417$** | **$9.349350 \times 10^{-7}$** |
| **`frac_rbh_desert`** | **$-0.137484$** | **$0.019456$** | **$0.019167$** | **$-7.0664$** | **$1.879206 \times 10^{-11}$** | **$[-0.175819,\; -0.099149]$** | **$-0.105413$** | **$0.024346$** | **$-4.3298$** | **$2.228562 \times 10^{-5}$** |

- **External verification against `statsmodels.GLS` (`results/v3_pgls_statsmodels_verification.tsv`):** Across $\lambda \in \{0.0, 0.5, 0.91578426, 0.91580000\}$, our Cholesky GLS solver agrees with `statsmodels.api.GLS(y, X, sigma=Sigma)` to within $\max|\Delta\hat{\beta}| = 5.16 \times 10^{-12}$, $\max|\Delta\mathrm{SE}_{\mathrm{unb}}| = 3.15 \times 10^{-13}$, $\max|\Delta R^2_{\mathrm{GLS}}| = 4.44 \times 10^{-15}$, and $\max|\Delta\hat{\sigma}^2_{\mathrm{unb}}| < 10^{-15}$. At rounded $\lambda = 0.9158$, $R^2_{\mathrm{GLS}} = 0.46905148$, matching `reviews/darkmatter-2026-10-13/round6/root_checks/h1_v2_gls_audit.json` exactly.

### 4.3 Profile-Likelihood Uncertainty in $\lambda$ and Parameter Stability Across $\lambda \in [0, 1]$ (`results/v3_pgls_lambda_profile.tsv`)

To transparently report how uncertainty in $\lambda$ affects inference without pretending that fixed-$\lambda$ Student-$t$ intervals propagate $\lambda$ variance, we evaluated the GLS model across fixed $\lambda$ values spanning $[0, 1]$ and at the exact endpoints of the 95% profile-likelihood confidence set $\{\lambda : 2[\ell(\hat{\lambda}_{\mathrm{ML}}) - \ell(\lambda)] \le \chi^2_{1, 0.95} = 3.841459\}$:

| $\lambda$ | Grid Role | $\ell(\lambda)$ | $R^2_{\mathrm{GLS}}$ | $\hat{\sigma}^2_{\mathrm{unb}}$ | $\hat{\beta}_{\mathrm{AAI}}$ ($\mathrm{SE}_{\mathrm{unb}}$) | Two-Sided $p_{\mathrm{AAI}}$ | Conditional 95% CI ($\beta_{\mathrm{AAI}}$) | $\hat{\beta}_{\mathrm{RBH}}$ ($\mathrm{SE}_{\mathrm{unb}}$) | Two-Sided $p_{\mathrm{RBH}}$ | Conditional 95% CI ($\beta_{\mathrm{RBH}}$) |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `0.000000` | OLS ($\lambda=0$) | $526.1537$ | $0.495786$ | $0.00071159$ | $+0.001931$ ($0.000383$) | $9.3494 \times 10^{-7}$ | $[+0.001177, +0.002686]$ | $-0.105413$ ($0.024346$) | $2.2286 \times 10^{-5}$ | $[-0.153382, -0.057443]$ |
| `0.250000` | Fixed $\lambda=0.25$ | $595.6966$ | $0.431436$ | $0.00046460$ | $+0.001978$ ($0.000350$) | $4.6901 \times 10^{-8}$ | $[+0.001288, +0.002668]$ | $-0.131750$ ($0.022043$) | $8.5988 \times 10^{-9}$ | $[-0.175182, -0.088317]$ |
| `0.500000` | Fixed $\lambda=0.50$ | $617.1670$ | $0.443609$ | $0.00050467$ | $+0.001927$ ($0.000355$) | $1.4725 \times 10^{-7}$ | $[+0.001227, +0.002628]$ | $-0.137759$ ($0.021083$) | $4.0702 \times 10^{-10}$ | $[-0.179300, -0.096218]$ |
| `0.750000` | Fixed $\lambda=0.75$ | $633.8459$ | $0.458605$ | $0.00066878$ | $+0.001854$ ($0.000371$) | $1.1449 \times 10^{-6}$ | $[+0.001123, +0.002585]$ | $-0.140358$ ($0.020155$) | $3.4362 \times 10^{-11}$ | $[-0.180070, -0.100646]$ |
| **`0.846303`** | **95% $\lambda$ CI Lower** | **$639.0122$** | **$0.465032$** | **$0.00082751$** | **$+0.001790$ ($0.000386$)** | **$5.8489 \times 10^{-6}$** | **$[+0.001030, +0.002550]$** | **$-0.139789$ ($0.019767$)** | **$1.8202 \times 10^{-11}$** | **$[-0.178736, -0.100841]$** |
| **`0.915784`** | **Profile MLE $\hat{\lambda}_{\mathrm{ML}}$** | **$640.9330$** | **$0.469051$** | **$0.00105887$** | **$+0.001704$ ($0.000407$)** | **$4.0650 \times 10^{-5}$** | **$[+0.000902, +0.002506]$** | **$-0.137484$ ($0.019456$)** | **$1.8792 \times 10^{-11}$** | **$[-0.175819, -0.099149]$** |
| **`0.957823`** | **95% $\lambda$ CI Upper** | **$639.0122$** | **$0.471087$** | **$0.00134948$** | **$+0.001613$ ($0.000434$)** | **$2.5117 \times 10^{-4}$** | **$[+0.000758, +0.002467]$** | **$-0.133697$ ($0.019258$)** | **$3.8983 \times 10^{-11}$** | **$[-0.171642, -0.095753]$** |
| `1.000000` | Strict Brownian ($\lambda=1$) | $608.6726$ | $0.503660$ | $0.00253978$ | $+0.001484$ ($0.000539$) | $6.3878 \times 10^{-3}$ | $[+0.000422, +0.002546]$ | $-0.123559$ ($0.020164$) | $3.8332 \times 10^{-9}$ | $[-0.163289, -0.083830]$ |

Across the entire 95% profile-likelihood interval $\lambda \in [0.846303, 0.957823]$, the union of conditional 95% Student-$t$ intervals is $[+0.000758, +0.002550]$ for $\beta_{\mathrm{AAI}}$ and $[-0.178736, -0.095753]$ for $\beta_{\mathrm{RBH}}$, bounded well away from zero throughout.

---

## 5. Pre-Registered Sensitivity Analyses (`S1–S7`, `results/v3_pgls_sensitivity.tsv`)

| Model ID & Specification | $n$ ($df$) | $\hat{\lambda}_{\mathrm{ML}}$ (95% Profile CI) | $R^2_{\mathrm{GLS}}$ | $\hat{\sigma}^2_{\mathrm{unb}}$ | $\hat{\beta}_{\mathrm{AAI}}$ ($\mathrm{SE}_{\mathrm{unb}}$) | Two-Sided $p_{\mathrm{AAI}}$ | Conditional 95% CI ($\beta_{\mathrm{AAI}}$) | $\hat{\beta}_{\mathrm{RBH}}$ ($\mathrm{SE}_{\mathrm{unb}}$) | Two-Sided $p_{\mathrm{RBH}}$ | Conditional 95% CI ($\beta_{\mathrm{RBH}}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`S1_Primary_All_237`** (Normalized VCV, native RBH denom) | $237$ ($230$) | $0.915784$ $[0.846303, 0.957823]$ | $0.469051$ | $0.00105887$ | $+0.001704$ ($0.000407$) | $4.065008 \times 10^{-5}$ | $[+0.000902, +0.002506]$ | $-0.137484$ ($0.019456$) | $1.879206 \times 10^{-11}$ | $[-0.175819, -0.099149]$ |
| **`S2_Unnormalized_Brownian_VCV_237`** (Raw patristic VCV $\mathbf{V}$) | $237$ ($230$) | $0.916132$ $[0.846418, 0.958344]$ | $0.475264$ | $0.00106074$ | $+0.001744$ ($0.000407$) | $2.694748 \times 10^{-5}$ | $[+0.000942, +0.002546]$ | $-0.140505$ ($0.019533$) | $8.831619 \times 10^{-12}$ | $[-0.178990, -0.102019]$ |
| **`S3_Catalogue_Denominator_RBH_237`** ($N_{\mathrm{RBH}}/N_{\mathrm{proteins,cat}}$) | $237$ ($230$) | $0.915415$ $[0.845649, 0.957628]$ | $0.468174$ | $0.00105880$ | $+0.001693$ ($0.000407$) | $4.499514 \times 10^{-5}$ | $[+0.000891, +0.002495]$ | $-0.136591$ ($0.019416$) | $2.265112 \times 10^{-11}$ | $[-0.174848, -0.098334]$ |
| **`S4_Restricted_c__SW-10-69-26_185`** (Within `c__SW-10-69-26` only) | $185$ ($178$) | $0.900804$ $[0.800104, 0.960065]$ | $0.500906$ | $0.00097447$ | $+0.001921$ ($0.000491$) | $1.285632 \times 10^{-4}$ | $[+0.000953, +0.002889]$ | $-0.180441$ ($0.022624$) | $1.776834 \times 10^{-13}$ | $[-0.225087, -0.135795]$ |
| **`S5_Dominant_Family_JAJPHT01_Omitted_97`** (Excluding 140 `f__JAJPHT01` reps) | $97$ ($90$) | $0.901145$ $[0.765888, 0.963229]$ | $0.444105$ | $0.00103139$ | $+0.001494$ ($0.000513$) | $4.526286 \times 10^{-3}$ | $[+0.000475, +0.002514]$ | $-0.102379$ ($0.030547$) | $1.176333 \times 10^{-3}$ | $[-0.163065, -0.041692]$ |
| **`S6_Same_Class_Comparators_SW_plus_E2_236`** (Excluding 1 `c__Thermoplasmata` rep) | $236$ ($229$) | $0.893376$ $[0.806220, 0.946533]$ | $0.475040$ | $0.00091782$ | $+0.001568$ ($0.000392$) | $8.387793 \times 10^{-5}$ | $[+0.000797, +0.002340]$ | $-0.145879$ ($0.019191$) | $7.485088 \times 10^{-13}$ | $[-0.183693, -0.108064]$ |
| **`S7_With_Caller_Complete_Proxy_232`** (Adding CheckM2 `frac_complete`, $p=8$) | $232$ ($224$) | $0.913489$ $[0.839117, 0.958295]$ | $0.483091$ | $0.00104168$ | $+0.001550$ ($0.000452$) | $7.245799 \times 10^{-4}$ | $[+0.000659, +0.002441]$ | $-0.104756$ ($0.026891$) | $1.293680 \times 10^{-4}$ | $[-0.157748, -0.051764]$ |

---

## 6. External Comparator Provenance, Empirical AAI Ranges, and Gene-Caller Reconciliation (`C4`, `C5`)

### 6.1 Exact Comparator Provenance and Empirical AAI Distributions (`results/v3_comparator_metadata_summary.tsv`)
- **Comparator pool and nearest-neighbor provenance:** In the upstream pairwise AAI pipeline (`calculate_aai.py` -> `thermo_aai_summary.tsv`), the 237 desert *Thermoplasmatota* species representatives were compared against `247` external reference genomes belonging to classes `c__SW-10-69-26` and `c__E2` (`218` `GTDB_r232_rep` and `29` non-duplicate `MGnify` MAGs). Across the 237 desert representatives, the nearest external comparator (`nearest_ref` maximizing `aai_pct`) is a `GTDB_r232_rep` genome for `229 / 237` (`96.6%`) representatives and an `MGnify` MAG for `8 / 237` (`3.4%`) representatives (`32` unique external accessions).
- **Habitat metadata of nearest comparators:** Cross-referencing `thermo_tree_taxa.tsv` shows that `142 / 237` (`59.9%`) nearest comparators carry `Environmental:Terrestrial:Soil` biome metadata, `3 / 237` (`1.3%`) carry `Environmental:Aquatic:Marine:Sediment` metadata, and `92 / 237` (`38.8%`) are GTDB r232 reference genomes without MGnify habitat strings. None are verified pure-culture laboratory isolates.
- **Exact empirical AAI distributions by taxonomic stratum:**

| Taxonomic Stratum | $n$ | Nearest `GTDB_r232_rep` / `MGnify` (Unique) | Biome: Soil / Marine Sed. / Unannot. | Min AAI (%) | Q25 AAI (%) | Median AAI (%) | Mean AAI (%) | Q75 AAI (%) | Max AAI (%) | Median `frac_rbh_desert` (`aai` / `cat` denom) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`All_Thermoplasmatota_237`** | $237$ | $229\;/\;8$ ($32$) | $142\;/\;3\;/\;92$ | **$48.678$** | $58.763$ | $64.774$ | $63.727$ | $66.392$ | **$80.041$** | $0.409236\;/\;0.413209$ |
| **`Same_Class_Comparators_SW_plus_E2_236`** | $236$ | $228\;/\;8$ ($32$) | $141\;/\;3\;/\;92$ | **$50.389$** | $58.806$ | $64.781$ | $63.791$ | $66.396$ | **$80.041$** | $0.410209\;/\;0.413318$ |
| **`Class_c__SW-10-69-26_185`** | $185$ | $185\;/\;0$ ($15$) | $121\;/\;0\;/\;64$ | **$50.389$** | $57.656$ | $64.794$ | $63.334$ | $66.172$ | **$80.041$** | $0.416028\;/\;0.418299$ |
| **`Class_c__E2_51`** | $51$ | $43\;/\;8$ ($17$) | $20\;/\;3\;/\;28$ | **$58.763$** | $61.609$ | $64.185$ | $65.450$ | $67.731$ | **$76.964$** | $0.356624\;/\;0.357761$ |
| **`Class_c__Thermoplasmata_1`** | $1$ | $1\;/\;0$ ($1$) | $1\;/\;0\;/\;0$ | **$48.678$** | $48.678$ | $48.678$ | $48.678$ | $48.678$ | **$48.678$** | $0.276248\;/\;0.279252$ |

### 6.2 Gene-Caller Reconciliation and Non-Disjointness of $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ (`C4`)
- **Dual gene-caller provenance:** `n_proteins` (`n_proteins_cat`, median $3302.0$, range $[1481, 5388]$) is the protein count from the unified catalogue gene call (`v2_reps.tsv`), whereas `n_prot_aai` (median $3312.0$, range $[1485, 5469]$) is the protein count from fresh single-genome `pyrodigal.GeneFinder(meta=False)` execution in `v2_thermo_relatives_aai.py`. Across all 237 genomes, $N_{\mathrm{prot,aai}} - N_{\mathrm{proteins,cat}}$ has mean $+19.52$, median $+18.0$, range $[-13, +93]$, and Pearson $r = 0.999752$. Across the `232 / 237` genomes with archived CheckM2 `.faa` files (`5` `48PRr1__SemiBin_*` representatives—`112`, `141`, `510`, `539`, `849`—lack `.faa` files in `_checkm2_refined/protein_files/`), `n_prot_aai` matches the exact count of FASTA headers (`max |n_prot_aai - n_proteins_checkm2| = 0`).
- **Why disjointness of $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ is neither established nor assumed:** $N_{\mathrm{dark}}$ counts catalogue genes whose 90%-identity / 80%-coverage MMseqs2 cluster representative lacks any hit in UniRef50 (release `2026_01`) or functional databases (eggNOG, Pfam, KOfam, Bakta). Conversely, $N_{\mathrm{RBH}}$ (`aai_n_rbh`) counts single-genome Pyrodigal calls forming a DIAMOND reciprocal best hit ($\ge 30\%$ identity, $\ge 70\%$ query coverage, $e \le 10^{-5}$) against the nearest external GTDB r232 or MGnify reference genome ($145 / 237$ nearest pairs carry MGnify environmental biome metadata and $92 / 237$ are GTDB-only without cultivation metadata). Because the external comparator set includes environmental MAGs containing uncharacterized proteins, a desert MAG gene can be unannotated in UniRef50/Pfam/KOfam (counted in $N_{\mathrm{dark}}$) and simultaneously share an uncharacterized homolog in a related GTDB/MGnify reference genome (counted in $N_{\mathrm{RBH}}$). Thus disjointness of $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ ($N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$) is neither established nor assumed.

### 6.3 Matched-Denominator Algebraic Permutation Illustration and Descriptive Count Scaling (`results/v3_coupling_and_comparator_audit.tsv`)

To illustrate the purely mathematical behavior of dividing two positive count variables by an identical varying denominator $D$, we compared the observed Spearman rank correlation $\rho(N_{\mathrm{dark}}/D, N_{\mathrm{RBH}}/D)$ against `10,000` independent random permutations of the numerators ($\pi_1(N_{\mathrm{dark}})/D$ vs $\pi_2(N_{\mathrm{RBH}})/D$, `seed = 42`) using **strictly matched denominators** (`D = n_proteins_cat` for both, and `D = n_prot_aai` for both):

| Cohort Subset | Comparison / Denominator Basis | Observed or Null Metric | Spearman $\rho$ | Two-Sided / Plus-One $p$-Value | 95% Permutation Interval / Interpretation |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **`c__SW-10-69-26` ($n=185$)** | Observed ratio: native denominators (`n_dark/n_proteins_cat` vs `aai_n_rbh/n_prot_aai`) | Observed $\rho$ | $-0.394176$ | $2.833078 \times 10^{-8}$ | Primary PGLS variable definitions |
| **`c__SW-10-69-26` ($n=185$)** | **Observed ratio: identical catalogue denominator (`n_proteins_cat`)** | **Observed $\rho$** | **$-0.395303$** | **$2.565575 \times 10^{-8}$** | **Exact denominator match to catalogue permutation** |
| **`c__SW-10-69-26` ($n=185$)** | **Algebraic permutation: identical catalogue denominator (`n_proteins_cat`)** | **Null mean $\rho$** | **$+0.301737$** | **$9.999000 \times 10^{-5}$** | **$[+0.1789, +0.4194]$ (algebraic illustration only)** |
| **`c__SW-10-69-26` ($n=185$)** | **Observed ratio: identical Pyrodigal denominator (`n_prot_aai`)** | **Observed $\rho$** | **$-0.381418$** | **$8.481197 \times 10^{-8}$** | **Exact denominator match to Pyrodigal permutation** |
| **`c__SW-10-69-26` ($n=185$)** | **Algebraic permutation: identical Pyrodigal denominator (`n_prot_aai`)** | **Null mean $\rho$** | **$+0.302959$** | **$9.999000 \times 10^{-5}$** | **$[+0.1807, +0.4199]$ (algebraic illustration only)** |
| **`c__SW-10-69-26` ($n=185$)** | Historical distorted $\min(N_A, N_B)$ denominator (`aai_frac_rbh`) | Observed $\rho$ | $-0.301017$ | $3.134038 \times 10^{-5}$ | Uses external comparator $N_B$ in $74.7\%$ of pairs |
| **`c__SW-10-69-26` ($n=185$)** | Raw unnormalized counts (`n_dark` vs `aai_n_rbh`) | Observed $\rho$ | $+0.241889$ | $9.094245 \times 10^{-4}$ | Positive raw count association across genomes |
| **`c__SW-10-69-26` ($n=185$)** | Descriptive count scaling: `n_dark` vs `n_proteins_cat` | Observed $\rho$ | $+0.945740$ | $2.819623 \times 10^{-91}$ | Descriptive scaling with catalogue proteome size |
| **`c__SW-10-69-26` ($n=185$)** | Descriptive count scaling: `aai_n_rbh` vs `n_proteins_cat` | Observed $\rho$ | $+0.329415$ | $4.680864 \times 10^{-6}$ | Descriptive scaling with catalogue proteome size |
| **`c__SW-10-69-26` ($n=185$)** | Descriptive count scaling: `aai_n_rbh` vs `n_prot_aai` | Observed $\rho$ | $+0.322727$ | $7.454653 \times 10^{-6}$ | Descriptive scaling with Pyrodigal proteome size |
| **`All Thermoplasmatota` ($n=237$)** | Observed ratio: native denominators (`n_dark/n_proteins_cat` vs `aai_n_rbh/n_prot_aai`) | Observed $\rho$ | $-0.246356$ | $1.272235 \times 10^{-4}$ | Primary PGLS variable definitions |
| **`All Thermoplasmatota` ($n=237$)** | **Observed ratio: identical catalogue denominator (`n_proteins_cat`)** | **Observed $\rho$** | **$-0.248353$** | **$1.116208 \times 10^{-4}$** | **Exact denominator match to catalogue permutation** |
| **`All Thermoplasmatota` ($n=237$)** | **Algebraic permutation: identical catalogue denominator (`n_proteins_cat`)** | **Null mean $\rho$** | **$+0.305434$** | **$9.999000 \times 10^{-5}$** | **$[+0.1966, +0.4114]$ (algebraic illustration only)** |
| **`All Thermoplasmatota` ($n=237$)** | **Observed ratio: identical Pyrodigal denominator (`n_prot_aai`)** | **Observed $\rho$** | **$-0.234877$** | **$2.643639 \times 10^{-4}$** | **Exact denominator match to Pyrodigal permutation** |
| **`All Thermoplasmatota` ($n=237$)** | **Algebraic permutation: identical Pyrodigal denominator (`n_prot_aai`)** | **Null mean $\rho$** | **$+0.305432$** | **$9.999000 \times 10^{-5}$** | **$[+0.1956, +0.4118]$ (algebraic illustration only)** |
| **`All Thermoplasmatota` ($n=237$)** | Historical distorted $\min(N_A, N_B)$ denominator (`aai_frac_rbh`) | Observed $\rho$ | $-0.214637$ | $8.820538 \times 10^{-4}$ | Uses external comparator $N_B$ in $74.7\%$ of pairs |
| **`All Thermoplasmatota` ($n=237$)** | Raw unnormalized counts (`n_dark` vs `aai_n_rbh`) | Observed $\rho$ | $+0.495214$ | $4.554587 \times 10^{-16}$ | Positive raw count association across genomes |
| **`All Thermoplasmatota` ($n=237$)** | Descriptive count scaling: `n_dark` vs `n_proteins_cat` | Observed $\rho$ | $+0.962319$ | $6.723093 \times 10^{-135}$ | Descriptive scaling with catalogue proteome size |
| **`All Thermoplasmatota` ($n=237$)** | Descriptive count scaling: `aai_n_rbh` vs `n_proteins_cat` | Observed $\rho$ | $+0.577748$ | $1.677587 \times 10^{-22}$ | Descriptive scaling with catalogue proteome size |
| **`All Thermoplasmatota` ($n=237$)** | Descriptive count scaling: `aai_n_rbh` vs `n_prot_aai` | Observed $\rho$ | $+0.573625$ | $3.890989 \times 10^{-22}$ | Descriptive scaling with Pyrodigal proteome size |

> **Methodological interpretation (`C4`):** Permuting $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ independently across genomes breaks their empirical covariance with proteome size and is presented **solely as an algebraic illustration** showing that dividing two independent positive counts by an identical genome-size denominator $D$ induces a **positive** correlation ($\rho \approx +0.302$ in `c__SW-10-69-26`, $\rho \approx +0.305$ across all $237$ reps), not a negative one. Likewise, the stronger positive rank correlation of $N_{\mathrm{dark}}$ with total proteome size ($\rho = +0.945740$ in `c__SW-10-69-26`) relative to $N_{\mathrm{RBH}}$ ($\rho = +0.329415$) is a **descriptive empirical summary** of how conserved versus unannotated gene counts co-vary with proteome size in this cohort, not a mechanistic or causal proof.

### 6.4 Descriptive Accessory Dark Family Occupancy (`O2`)
Across the `14,645` admitted dark protein clusters (`results/reconciled_dark_family_taxonomy_summary.tsv`), `c__SW-10-69-26` species representatives ($n=185$) carry a median of **`244.0` distinct dark gene families per genome** (IQR `153.0–363.0`, mean `263.49`), compared with **`77.0` distinct dark gene families per genome** (IQR `53.0–115.0`, mean `87.84`) in `c__E2` ($n=51$). The two-sided Mann-Whitney $U$ statistic returned by `scipy.stats.mannwhitneyu` is **$U = 8292.5$** (complementary smaller statistic $U_{\min} = 185 \times 51 - 8292.5 = 1142.5$, two-sided $p = 1.224554 \times 10^{-16}$; correcting the v2 transcription error `1,180.5`). This cross-sectional difference in dark family repertoire breadth is reported strictly as a descriptive summary without inferring ancestral expansion or contraction polarity.

---

## 7. Execution Provenance, Slurm Job Record, and File Inventory

### 7.1 IBEX Slurm Execution Record
- **Batch script:** `scripts/run_h1_v3_pipeline.sbatch`
- **Slurm Job ID:** `53704703` (Partition: `debug`, Node: `gpu510-32`, CPUs: `8`, Memory: `32G`, Walltime elapsed: `00:00:42` / `42.65s`, State: `COMPLETED`, Exit Code: `0:0`)
- **Software environments:**
  - Step 1 (`01_v3_input_manifest_and_phylo_audit.py`): `/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python` (`Python 3.12.12`, `DendroPy 5.0.8`, `numpy 2.4.2`, `scipy 1.17.0`, `pandas 3.0.1`)
  - Steps 2–3 (`02_v3_pgls_ancestry_models.py`, `03_v3_coupling_and_comparator_audit.py`): `/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python` (`Python 3.10.19`, `numpy 2.2.6`, `scipy 1.15.2`, `pandas 2.3.3`, `statsmodels 0.14.6`)
- **Logs:** `logs/h1_v3_53704703.out`, `logs/h1_v3_53704703.err`

### 7.2 Complete v3 File Inventory (`reproducibility/darkmatter/round6_hypotheses/h1_lineage/`)
- **Pre-registered Design & Dossier:**
  - `DESIGN.md` — Pre-registered v3 experimental design (written prior to v3 computation; v2 preserved at `DESIGN_v2_superseded.md`).
  - `DOSSIER.md` — This v3 scientific dossier (v2 preserved at `DOSSIER_v2_superseded.md`).
  - `response.txt` — Point-by-point mapping of Round 6 review items `C1–C5` and `O1–O2` to executed v3 scripts and output tables.
  - `INPUTS.tsv` — Enumerated SHA-256 checksum manifest of all `247` raw input entries (`10` core tables/trees + `232` readable CheckM2 `.faa` files + `5` explicitly documented missing `.faa` entries; v2 preserved at `INPUTS_v2_superseded.tsv`).
- **v3 Analysis Scripts (`scripts/`):**
  - `scripts/01_v3_input_manifest_and_phylo_audit.py`
  - `scripts/02_v3_pgls_ancestry_models.py`
  - `scripts/03_v3_coupling_and_comparator_audit.py`
  - `scripts/run_h1_v3_pipeline.sbatch`
- **v3 Output Tables and Reports (`results/`):**
  - `results/v3_input_file_manifest.tsv` — Identical copy of top-level `INPUTS.tsv` (`247` rows).
  - `results/v3_phylo_vcv_validation.json` & `results/v3_phylo_vcv_validation.txt` — Tree tip uniqueness (`1,914`), cohort mapping (`237/237`), outgroup vs midpoint rooting equivalence (`max_abs_diff = 0.0`), symmetry, and eigenvalue audit.
  - `results/v3_pgls_model_summary.tsv` — Primary PGLS (`lambda_ML = 0.915784`) and OLS (`lambda = 0.0`) parameter table with both `SE_unb` and `SE_ML`, two-sided Student-$t$ $p$-values, and conditional 95% Student-$t$ CIs.
  - `results/v3_pgls_lambda_profile.tsv` — Fixed-$\lambda$ GLS evaluations across $\lambda \in \{0.0, 0.25, 0.50, 0.75, 0.846303, 0.915784, 0.957823, 1.0\}$.
  - `results/v3_pgls_sensitivity.tsv` — Full summary of pre-registered models `S1` through `S7`.
  - `results/v3_pgls_statsmodels_verification.tsv` — Independent verification of `beta`, `SE_unb`, `t`, `p`, `R^2_GLS`, `sigma^2_unb`, and `loglik` against `statsmodels.GLS`.
  - `results/v3_pgls_ancestry_summary.txt` — Formatted statistical report for `S1–S7`.
  - `results/v3_comparator_metadata_summary.tsv` — Exact source, biome, and empirical AAI quantiles across taxonomic strata.
  - `results/v3_coupling_and_comparator_audit.tsv` & `results/v3_coupling_and_comparator_summary.txt` — Matched-denominator permutation illustrations, descriptive count scaling, and Mann-Whitney $U$ dark family audit.
  - `results/v3_closure_verification.json` — Programmatic assertion log verifying all mathematical, statistical, and cohort-closure invariants.
