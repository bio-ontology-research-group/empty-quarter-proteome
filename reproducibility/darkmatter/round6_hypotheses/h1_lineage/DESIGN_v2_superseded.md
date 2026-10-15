# DESIGN: Hypothesis 1 — What Explains Within-Thermoplasmatota Darkness?

## 1. Scientific Objective and Context

In the Rub' al-Khali desert metagenome-assembled genome (MAG) catalogue, the archaeal phylum *Thermoplasmatota* exhibits an exceptionally high dark proteome fraction (median 30.96%, ~4.5-fold enriched in dark families relative to its species count). In GTDB-Tk r232, the 237 desert *Thermoplasmatota* species representatives are assigned predominantly to two provisional classes:
- `c__SW-10-69-26` ($n = 185$ representatives, spanning 4 families across 2 orders)
- `c__E2` ($n = 51$ representatives, spanning 3 families in 1 order)
- `c__Thermoplasmata` ($n = 1$ representative, family `f__Methanomassiliicoccaceae`)

In Round 5, an exploratory review (Argon) proposed three major empirical assertions regarding within-phylum darkness:
1. External AAI has virtually zero correlation with dark fraction across the 237 genomes (Spearman $\rho = 0.012$, $p = 0.85$), which was interpreted as refuting external reference divergence as an explanatory factor.
2. `c__SW-10-69-26` has ~1,000 more proteins per genome and a 3.7 percentage-point higher dark fraction than `c__E2`, which was claimed to establish an "accessory proteome expansion" in `c__SW-10-69-26`.
3. Within `c__SW-10-69-26`, dark fraction correlates negatively with external reciprocal best hit (RBH) fraction ($\rho = -0.301$, $p = 3.1 \times 10^{-5}$), which was interpreted as biological core dilution.

Following Root's critical audit (`ROOT_AUDIT_NOTES.md`), these exploratory assertions require rigorous restructuring:
- **Low Effective Lineage Replication & Pseudoreplication:** An asymptotic cluster-sandwich adjustment on 7 families or 3 orders is mathematically invalid due to extreme sample imbalance (`f__JAJPHT01` contains 140 / 237 genomes; `c__E2` contains only 1 order). Proper ancestry-aware analysis requires continuous branch-length phylogenetic covariance (PGLS) based on the ar53 FastTree, not arbitrary clade clustering.
- **Ratio Coupling vs Biological Dilution:** The negative ratio correlation between $N_{\mathrm{dark}} / N_{\mathrm{prot}}$ and $N_{\mathrm{RBH}} / N_{\mathrm{prot}}$ cannot be interpreted as biological core dilution without accounting for the finite proteome constraint ($N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$) and scaling rates. Conversely, an independent-numerator permutation with a common denominator produces a *positive* null correlation (+0.3017); calling the observed negative correlation "100% mathematical artifact" was a sign error in the initial draft that must be retracted.
- **Distorted RBH Denominator:** In `v2_thermo_relatives_aai.py`, `aai_frac_rbh` was defined as $N_{\mathrm{RBH}} / \min(N_a, N_b)$. Because desert MAGs are larger than nearest external comparators in 74.7% of pairs, the denominator was the *comparator's* proteome size, not the desert genome's.
- **Restriction of Range:** All 237 representatives are phylogenetically distant from known isolate genomes (median AAI 64.7%, near the genus/family boundary). Flat correlation in a truncated range cannot rule out reference bias across the broader tree of life.
- **Macroevolutionary Scope Limits:** Extra genes in larger genomes represent diverse accessory orthogroups and solitary families, but cannot establish evolutionary polarity (gene gain in `SW-10-69-26` vs gene reduction in `E2`) without outgroup-rooted ancestral state reconstruction.

---

## 2. Primary Identifiable Association Estimand and Model

### Mathematical Formulation
To eliminate two-class pseudoreplication, asymptotic cluster distortion, and unmodeled phylogenetic autocorrelation, we define **ONE primary identifiable association model**: Phylogenetic Generalized Least Squares (PGLS):

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim \mathcal{N}\left(\mathbf{0}, \sigma^2 \mathbf{V}(\lambda)\right)$$

where:
- $\mathbf{y} \in \mathbb{R}^{237}$ is the observed dark fraction ($Y_i = N_{\mathrm{dark}, i} / N_{\mathrm{prot}, i}$).
- $\mathbf{V}(\lambda) = \lambda \mathbf{C} + (1 - \lambda)\mathbf{I}_{237}$ is the phylogenetic variance-covariance structure parameterized by Pagel's $\lambda \in [0, 1]$.
- $\mathbf{C} \in \mathbb{R}^{237 \times 237}$ is the phylogenetic correlation matrix derived from the cophenetic branch lengths of the GTDB-Tk r232 ar53 FastTree (`thermo_tree.nwk`) pruned to the 237 desert representatives.
- $\mathbf{X} \in \mathbb{R}^{237 \times 7}$ contains the intercept and 6 a priori defined covariates:
  1. `aai`: continuous amino acid identity to the nearest external GTDB r232 reference genome (percentage units).
  2. `frac_rbh_desert`: shared RBH coverage normalized by the desert MAG protein count ($N_{\mathrm{RBH}} / N_{\mathrm{prot, desert}}$).
  3. `completeness`: CheckM2 completeness percentage.
  4. `contamination`: CheckM2 contamination percentage.
  5. `log_len`: $\log_{10}(\mathrm{assembly\ length\ in\ bp})$.
  6. `log_n50`: $\log_{10}(\mathrm{contig\ N50\ in\ bp})$.

### Primary Estimands
The primary estimands are the conditional partial regression slopes:
1. **External Divergence Slope:** $\beta_{\mathrm{AAI}}$ (change in dark fraction per 1 percentage-point increase in nearest-relative AAI, holding shared gene coverage, assembly length, fragmentation, completeness, and phylogenetic covariance constant).
2. **Shared Coverage Slope:** $\beta_{\mathrm{frac\_rbh\_desert}}$ (change in dark fraction per unit increase in shared RBH coverage, conditional on AAI, quality, and phylogenetic covariance).

### Model Fitting and $\lambda$ Profiling
The model is fit via Cholesky decomposition of $\mathbf{V}(\lambda) = \mathbf{L}\mathbf{L}^T$, transforming the system to $\tilde{\mathbf{y}} = \mathbf{L}^{-1}\mathbf{y}$ and $\tilde{\mathbf{X}} = \mathbf{L}^{-1}\mathbf{X}$. Pagel's $\lambda$ is estimated by profile maximum likelihood across $\lambda \in [0, 1]$. Model-conditional 95% confidence intervals and two-sided Wald $t$-tests ($df = n - p = 230$) are reported.

---

## 3. Populations and Sensitivity Analyses (A Priori Defined)

To guard against model misspecification, extreme clade leverage, and unmeasured caller artifacts, we pre-specify one primary population and four sensitivity checks:

1. **Primary Cohort ($n = 237$):**
   All 237 desert *Thermoplasmatota* species representatives. Full sample retention ($100\%$).
2. **Restricted Cohort (`c__SW-10-69-26` Only, $n = 185$):**
   Omits the 51 `c__E2` and 1 `c__Thermoplasmata` representatives to test whether associations are driven by between-class leverage or hold strictly within the dominant class.
3. **Dominant Family Omission (`f__JAJPHT01` Omitted, $n = 97$):**
   Family `f__JAJPHT01` contains 140 of the 237 genomes (all within `SW-10-69-26`). Omission evaluates whether findings depend on this single dense lineage or persist across the remaining 6 families.
4. **Phylogenetic Signal Grid ($\lambda \in \{0.00, 0.25, 0.50, 0.75, 1.00\}$ vs $\hat{\lambda}_{\mathrm{ML}}$):**
   Evaluates slope stability from ordinary least squares ($\lambda = 0$, independence) to strict Brownian motion ($\lambda = 1$).
5. **Caller Complete-Gene Proxy Sensitivity ($n = 232$):**
   Adds `frac_complete` (fraction of genes with `partial=00` from CheckM2 single-genome Pyrodigal calls) as an explicit fragmentation/caller quality proxy. Sample loss of 5 MAGs ($n = 232$) is honestly reported; kept strictly as a sensitivity model to avoid conflating gene-caller sets in the primary specification.

---

## 4. Reconciled Ratio Coupling and Denominator Audit

### RBH Denominator Reconciliation
We audit three distinct denominators for reciprocal best hit fraction:
1. Argon published: $N_{\mathrm{RBH}} / \min(N_a, N_b)$ (distorted by external comparator size).
2. Desert-normalized: $N_{\mathrm{RBH}} / N_{\mathrm{desert}}$ (true fraction of desert genome proteome with external RBH).
3. Comparator-normalized: $N_{\mathrm{RBH}} / N_{\mathrm{nearest}}$ (fraction of external relative proteome shared).

### Mathematical Coupling vs Biological Mechanism
To rigorously evaluate the correlation between $N_{\mathrm{dark}} / N_{\mathrm{prot}}$ and $N_{\mathrm{RBH}} / N_{\mathrm{prot}}$:
- **Algebraic Null Illustration:** Permute $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ independently across observed $N_{\mathrm{prot}}$ values (10,000 iterations) with $(k+1)/(N+1)$ plus-one empirical tail correction. This demonstrates the pure effect of dividing independent numerators by a common denominator (which induces a *positive* correlation).
- **Mutual Exclusivity Constraint:** Acknowledge that $N_{\mathrm{dark}}$ and $N_{\mathrm{RBH}}$ are mutually exclusive subsets of the proteome ($N_{\mathrm{dark}} + N_{\mathrm{RBH}} \le N_{\mathrm{prot}}$). When total CDS count expands, if dark genes scale faster ($\rho = 0.95$) than RBH genes ($\rho = 0.33$), the ratio correlation becomes negative even though raw counts are positively correlated.
- **Explicit Scope:** Raw count positivity does not imply the ratio relation is purely an artifact; rather, the negative ratio reflects differential allometric scaling within a bounded sum, not independent competitive dilution.

---

## 5. Descriptive Accessory Dark Family Protocol

To inspect the genetic basis of genome size variation without asserting macroevolutionary polarities:
1. Map 124,022 contigs of the 237 representatives to 14,645 unique admitted dark families using `contig_to_bin_refined.tsv.gz` and `member_fam_final.tsv`.
2. Stratify families into genome-sharing tiers:
   - Singletons: present in exactly 1 genome.
   - Rare accessory: 2 to 11 genomes ($< 5\%$).
   - Intermediate accessory: 12 to 47 genomes (5% to 20%).
   - Frequent accessory: 48 to 118 genomes (20% to 50%).
   - Core-like: $\ge 119$ genomes ($\ge 50\%$).
3. Compare family distributions between `SW-10-69-26` and `E2` descriptively.
4. **Scope Limitation:** Singleton families are families observed in 1 genome in this dataset; they include rare lineage-specific genes, mobile elements, and unclustered homologs alongside potential fragments. Differences in family counts describe accessory richness in the sampled cohort, but do not establish whether `SW-10-69-26` expanded or `E2` underwent evolutionary reduction.

---

## 6. Execution and Verification Architecture

- **Phylogenetic Extraction:** DendroPy 5.0.8 in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python` prunes `thermo_tree.nwk` to the 237 tips, computes cophenetic distance, root distance, and VCV correlation matrix $\mathbf{C}$.
- **PGLS Modeling:** Python in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python` performs Cholesky transformation, profile likelihood maximization, and sensitivity fitting.
- **IBEX Slurm Job:** Script `run_h1_v2_pipeline.py` orchestrated via `run_h1_v2_pipeline.sbatch` on partition `debug`.
- **Reproducibility Guarantee:** All scripts, seeds, commands, job IDs, input hashes, and output tables are preserved.
