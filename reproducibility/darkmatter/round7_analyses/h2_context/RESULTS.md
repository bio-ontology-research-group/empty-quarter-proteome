# RESULTS: Full-Catalogue Flanking Gene-Neighbourhood Recurrence Analysis (Hypothesis 2, v3 Protocol)

**Analysis Target:** Full Census of all 1,461 Multi-Species ($\ge 2$ dRep 95% ANI Species Clusters) Confidence-Qualified Monomer Structural Dark Families ($570 \, K$, $891 \, Q \setminus K$)  
**Target Venue:** *Nature Communications* (Supporting Supplementary / Article Findings)  
**Execution Environment:** IBEX Slurm Cluster (`cn603-16-l`, Job ID `53738443`, 16 CPUs, 64 GB RAM, runtime 5 min 21 s, exit code `0:0`)  
**Working Directory:**  
- Remote: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/`  
- Local: `reproducibility/darkmatter/round7_analyses/h2_context/`  
**Input Checksums:** 10,775 verified input files pinned in `results/executed_inputs_manifest.tsv` and `INPUTS.tsv` (`sha256: 56fbc2b5e7b3c045e5a4b6491350ad78976d2b42ced43b9dfff3e2585d35408c`)  
**Output Checksums:** All 13 output tables pinned in `results/outputs_sha256.txt` (excluding self)

---

## 1. Executive Summary & Core Findings

We executed the full-census evaluation of flanking gene-neighbourhood KEGG KO recurrence across the entire catalogue of $N = 1{,}461$ multi-species structural families admitted at Gate 3 ($570 \, K$ low-query-TM families with $\text{qTM} < 0.5$ vs. PDB and AFDB; $891 \, Q \setminus K$ fold-matched sequence-dark families with $\text{qTM} \ge 0.5$), scaling from the 300-family pilot to the complete finite catalogue:

1. **Part (a) Robust Within-Arm Syntenic Enrichment Over Same-Contig Background ($p = 0.0002$ in both arms):**  
   Relative to a same-contig eligible pseudo-focal position reference ($B = 5{,}000$ Monte Carlo draws on each of the $1{,}380$ Tier 1 eligible loci across $7{,}473$ host contigs), both structural dark arms exhibit strong and highly significant within-arm flanking KO recurrence:
   - **Arm $K$ ($n = 173$ Tier 1 families):** Observed recurrence rate $\bar{C}_K^{\text{obs}} = 40.46\%$ (`70 / 173`) vs. same-contig expected rate $\bar{\mu}_K^{\text{pos}} = 8.71\%$ (excess $+31.75$ percentage points, $p = 0.0002$; mean continuous recurrence $\bar{R}_K^{\text{obs}} = 0.3785$ vs. expected $0.0614$, excess $+0.3171, p = 0.0002$).
   - **Arm $Q \setminus K$ ($n = 180$ Tier 1 families):** Observed recurrence rate $\bar{C}_{Q\setminus K}^{\text{obs}} = 32.22\%$ (`58 / 180`) vs. same-contig expected rate $\bar{\mu}_{Q\setminus K}^{\text{pos}} = 6.56\%$ (excess $+25.67$ percentage points, $p = 0.0002$; mean continuous recurrence $\bar{R}_{Q\setminus K}^{\text{obs}} = 0.2940$ vs. expected $0.0520$, excess $+0.2420, p = 0.0002$).
   - Within the 353 evaluable families, recurrence exceeds the chosen same-contig positional reference; this result does not establish a general property of all dark families.

2. **Part (b) Adjusted Null Contrast Between Low-Query-TM ($K$) and Fold-Matched ($Q \setminus K$) Families:**  
   Before covariate adjustment, Arm $K$ displays an apparent $+8.24$ percentage-point difference in raw recurrence ($40.46\%$ vs. $32.22\%$, Fisher exact $p = 0.1214$). However, Tier 1 $K$ families exhibit observational covariate differences, including longer focal proteins (mean $320.8$ aa vs. $233.1$ aa, $\text{SMD} = +0.6439$) and higher flanking KO annotation opportunity ($0.3839$ vs. $0.3554$, $\text{SMD} = +0.1288$), which correlate with baseline detection opportunity across families. These observational associations represent baseline covariate imbalances rather than demonstrated causal mechanisms or proven confounding:
   - **Primary Stratified Label Permutation ($B = 10{,}000$, 16/16 common-support strata, 100% coverage):** The Mantel-Haenszel weighted excess contrast is:
     $$\Delta_{\text{MH}}(E_{\text{fam}}) = +0.035358 \quad (\text{Null SD} = 0.042346, \; p = 0.4069)$$
     and on raw recurrence:
     $$\Delta_{\text{MH}}(C_{\text{fam}}) = +0.053850 \quad (\text{Null SD} = 0.046674, \; p = 0.2598)$$
   - **1:1 Strict Exact-Phylum + Tight-Length Matching ($96$ pairs):** Enforcing 100% exact GTDB-Tk R232 phylum group matching and tight length calipers yields an exact zero raw difference ($31.25\%$ vs. $31.25\%$, sign-flip $p = 1.0000$) and null excess contrast ($\Delta_{\text{strict}}(E_{\text{fam}}) = -0.005440, p = 0.9226$).
   - **Full-Rank Multivariable Regressions ($n = 353$, all converged):** OLS/HC3 on background-adjusted excess $E_{\text{fam}}$ indicates $\hat{\beta}_{\text{is\_k}} = +0.029446 \pm 0.042950$ ($p = 0.4930$), while MLE logistic regression on raw recurrence yields an odds ratio of $\text{OR} = 1.3233$ ($p = 0.3445$, 95% CI $[0.740, 2.366]$).
   - **Conclusion:** Once conditioned on species breadth, coarse host taxonomy, and flanking KO annotation opportunity, low-query-TM ($K$) and fold-matched ($Q \setminus K$) dark families do not show a statistically established difference in background-adjusted excess recurrence. The hypothesis of differential neighbourhood recurrence between low-query-TM and fold-matched sequence-dark families is not supported (consistent with the adjusted null). Importantly, this failure to reject the null holds conditionally under the identifying assumption of within-stratum exchangeability; label permutation evaluates conditional exchangeability given observed strata but does not prove mutual independence across families or genomes.

---

## 2. Complete Attrition Funnel and Admitted Cohort Denominators

Every stage of data reduction from the upstream confidence-qualified monomer structural universe ($N = 2{,}838$) to the evaluable tiers was tracked without silent drops:

| Attrition Gate / Cohort Stage | Total Families | Arm $K$ ($\text{qTM} < 0.5$) | Arm $Q \setminus K$ ($\text{qTM} \ge 0.5$) | Excluded Denominator Breakdown |
|---|---:|---:|---:|---|
| **Gate 0:** Confidence-qualified monomer structural families (`mean_plddt >= 0.70`, `length >= 100` aa) | **2,838** | 1,031 | 1,807 | Baseline confidence-qualified structural dark universe |
| **Gate 1:** $\ge 1$ 90%-sequence-cluster rep on a refined MAG contig (`contig_to_bin_refined.tsv.gz`) | **2,426** | 917 | 1,509 | 412 families with only unbinned 90%-rep contigs (114 $K$, 298 $Q\setminus K$) |
| **Gate 2:** $\ge 1$ 90%-rep member in a valid dRep 95% ANI species cluster (`Cdb.csv`) | **2,378** | 898 | 1,480 | 48 families only in unclustered refined bins (19 $K$, 29 $Q\setminus K$) |
| **Gate 3 (Admitted Full Census):** $\ge 2$ distinct dRep 95% ANI species clusters (`all_evaluable_families.tsv`) | **1,461** | **570** | **891** | 917 single-species families (`n_species == 1`: 328 $K$, 589 $Q\setminus K$) |
| **Tier 1 (Primary Gate):** $\ge 2$ distinct species with complete 4-flank exact-coordinate Bakta-mapped windows | **353** | **173** (30.35%) | **180** (20.20%) | **1,108 Ineligible / Missing (`c_fam = -1`: 397 $K$, 711 $Q\setminus K$)** |
| **Tier 2 (3'-Stop Sensitivity Gate):** $\ge 2$ distinct species with 4 exact or unique same-strand 3'-stop flanks | **637** | **260** (45.61%) | **377** (42.31%) | **824 Ineligible / Missing (`c_fam = -1`: 310 $K$, 514 $Q\setminus K$)** |
| **Tier 3 (Partial-Window Sensitivity Gate):** $\ge 2$ distinct species with $\ge 1$ exact-mapped flank | **1,320** | **521** (91.40%) | **799** (89.67%) | **141 Ineligible / Missing (`c_fam = -1`: 49 $K$, 92 $Q\setminus K$)** |

**Denominators and Locus Accounting:**
- Across the $1{,}461$ Gate 3 families, `scripts/01_select_cohort.py` emitted $7{,}604$ candidate member rows ($4{,}791$ in valid dRep species clusters, $224$ in unclustered refined bins, $2{,}589$ unbinned contigs).
- Outcome-independent primary locus selection selected exactly one primary locus per `(family, species)` pair, yielding **$4{,}542$ designated primary loci** and **$3{,}062$ non-primary candidate rows** ($4{,}542 + 3{,}062 = 7{,}604$).
- Among primary loci, **$1{,}380$** satisfied the Tier 1 primary gate (`eligible_locus_primary == 1`), residing on **$869$ unique contigs** across **$617$ refined bins** and **$568$ dRep species clusters**.
- Family-level admission requires $S_F^{\text{prim}} \ge 2$, admitting **$353$ Tier 1 evaluable families** ($173 \, K, 180 \, Q \setminus K$). All remaining $1{,}108$ families are coded `-1` (`MISSING`), never $0$.

---

## 3. Cryptographic Contig Sequence SHA-256 + Bakta Log Verification

Across the $7{,}604$ candidate member rows, $5{,}015$ reside on $2{,}126$ unique refined MAG bins:
- **Verified by Exact Nucleotide SHA-256 + Bakta `.log` `orig-id`:** **`2,041 / 2,126` (`96.00%`)** refined bins passed exact uppercase forward-strand nucleotide SHA-256 identity between `output_bins/{bin}.fa.gz` and `bakta/{bin}/{bin}.fna` plus explicit `qc: revised sequence: id=contig_N, orig-id=<raw_contig>` matching in `{bin}.log`.
- **Handling of Duplicate Contig Lengths:** Even though **`1,561 / 2,041` (`76.48%`)** verified bins contain multiple contigs with identical nucleotide lengths, sequence SHA-256 combined with `.log` parsing achieved **100% unique 1:1 mapping** with **zero hash mismatches** and **zero reverse-complement mismatches**.
- **Explicit Cataloguing of Missing-Annotation Bins:** Exactly **`85 / 2,126` (`4.00%`)** refined bins lacked Bakta/eggNOG annotation files on IBEX. All 85 bins are explicitly recorded in `results/bin_contig_verification.tsv` and `results/pilot_join_metrics.json` (`bins_missing_annotations`). All $108$ primary-locus flanks on these bins were classified as `source_file_unavailable` and excluded from callable complete windows.

---

## 4. Flank Join Audit and Seven-Category Mutual Exclusivity

Every single flanking position ($\{-2, -1, +1, +2\}$) across all $7{,}604$ candidate rows (**$30{,}416$ flanks total**) was classified into seven mutually exclusive categories:

| Flank Audit Category | All Rows ($7{,}604$ rows, $30{,}416$ flanks) | Designated Primary Loci ($4{,}542$ rows, $18{,}168$ flanks) | Non-Primary Rows ($3{,}062$ rows, $12{,}248$ flanks) | Tier 1 Eligible Primary Loci ($1{,}380$ rows, $5{,}520$ flanks) | Tier 1 Ineligible Primary Loci ($3{,}162$ rows, $12{,}648$ flanks) |
|---|---:|---:|---:|---:|---:|
| `mapped_KO` (exact `(contig, start, stop, strand)` + KO) | 4,319 (14.20%) | 4,106 (22.60%) | 213 (1.74%) | **2,060 (37.32%)** | 2,046 (16.18%) |
| `mapped_noKO` (exact `(contig, start, stop, strand)`, no KO) | 7,581 (24.92%) | 7,241 (39.86%) | 340 (2.78%) | **3,460 (62.68%)** | 3,781 (29.89%) |
| `stop3prime_matched_KO` (unique 3'-stop+strand+frame+overlap + KO) | 789 (2.59%) | 750 (4.13%) | 39 (0.32%) | 0 (0.00%) | 750 (5.93%) |
| `stop3prime_matched_noKO` (unique 3'-stop+strand+frame+overlap, no KO) | 1,507 (4.95%) | 1,451 (7.99%) | 56 (0.46%) | 0 (0.00%) | 1,451 (11.47%) |
| `coordinate_mismatch_unmapped` (no exact or unique 3'-stop match) | 1,288 (4.23%) | 1,202 (6.62%) | 86 (0.70%) | 0 (0.00%) | 1,202 (9.50%) |
| `source_file_unavailable` (unbinned or missing Bakta bin) | 6,409 (21.07%) | 108 (0.59%) | 6,301 (51.45%) | 0 (0.00%) | 108 (0.85%) |
| `contig_edge_truncated` (position falls outside $[0, N_{\mathrm{genes}}-1]$) | 8,523 (28.02%) | 3,310 (18.22%) | 5,213 (42.56%) | 0 (0.00%) | 3,310 (26.17%) |
| **Total Flanks Evaluated** | **30,416 (100.0%)** | **18,168 (100.0%)** | **12,248 (100.0%)** | **5,520 (100.0%)** | **12,648 (100.0%)** |

**Cross-Caller 3'-Stop Concordance:** Among the $3{,}403$ non-exact binned primary flanks ($750 + 1{,}451 + 1{,}202$), **$2{,}201$ ($64.68\%$)** share the exact 3'-stop coordinate, strand, reading frame (`start % 3`), and positive nucleotide overlap with a unique Bakta CDS on the same verified contig, rescuing $750$ KO-bearing flanks in the Tier 2 sensitivity cohort ($637$ evaluable families: $260 \, K, 377 \, Q \setminus K$).

---

## 5. Part (a): Within-Arm Recurrence Relative to Same-Contig Positional Reference

For each Tier 1 eligible locus $\ell$, pseudo-focal positions were sampled uniformly from gene indices $j \in \mathcal{I}_c^{\mathrm{prim}}$ on the same host contig whose four flanks $(j-2, j-1, j+1, j+2)$ all satisfy the identical Tier 1 complete-window gate:
- Across all $1{,}380$ Tier 1 primary loci, **$1{,}340 / 1{,}380$ ($97.10\%$)** reside on contigs with $M_c^{\mathrm{prim}} \ge 2$ eligible positions ($1{,}169$ have $M_c^{\mathrm{prim}} \ge 5$; median $18.0$, mean $42.68$, $\max = 795$). Only $40$ loci ($2.90\%$) have $M_c^{\mathrm{prim}} = 1$.
- Under $B = 5{,}000$ Monte Carlo draws per locus, both arms show strong, non-random within-arm syntenic enrichment relative to random positions on the same host contigs:

| Structural Arm | Evaluable Families ($N$) | Observed Rate $\bar{C}_{\mathrm{arm}}^{\mathrm{obs}}$ | Same-Contig Positional Expected $\bar{\mu}_{\mathrm{arm}}^{\mathrm{pos}}$ (Null SD) | Within-Arm Excess $\bar{E}_{\mathrm{arm}} = \bar{C}^{\mathrm{obs}} - \bar{\mu}^{\mathrm{pos}}$ | $p_{\mathrm{pos}}$ ($C_{\mathrm{fam}}$) | Observed Mean $\bar{R}_{\mathrm{arm}}^{\mathrm{obs}}$ | Positional Expected $\bar{\rho}_{\mathrm{arm}}^{\mathrm{pos}}$ | Continuous Excess $\bar{E}_{\mathrm{arm}}^{(R)}$ ($p_{\mathrm{pos}}$) | Multi-Position Subset ($M_c \ge 2$ in $\ge 2$ sp) Excess |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Low-qTM ($K$)** | 173 | **40.46%** (70/173) | 8.71% (1.66%) | **+31.75 pp** | **0.0002** | 0.3785 | 0.0614 | +0.3171 (**0.0002**) | +32.26 pp ($n = 168$) |
| **Fold-Matched ($Q \setminus K$)** | 180 | **32.22%** (58/180) | 6.56% (1.50%) | **+25.67 pp** | **0.0002** | 0.2940 | 0.0520 | +0.2420 (**0.0002**) | +25.85 pp ($n = 168$) |

---

## 6. Part (b): Between-Arm Contrast ($K$ vs. $Q \setminus K$) Across Three Analysis Frameworks

### 6.1 Primary Estimator: Covariate-Stratified Label Permutation ($B = 10{,}000$)

Evaluable families were stratified jointly by `species_breadth_bin` (`sp2` vs. `sp3plus`), verified GTDB-Tk R232 `phylum_group` (`Actinomycetota`, `Thermoplasmatota`, `other_bacteria_archaea`, `multi`), and `ko_opportunity_bin` (`<= median` vs. `> median`):
- **100.0% Common-Support Coverage:** All $16 / 16$ joint strata have common support ($n_{s, K} \ge 1$ and $n_{s, Q\setminus K} \ge 1$), covering **$100.0\%$ of families** in both Tier 1 ($173/173 \, K, 180/180 \, Q\setminus K$) and Tier 2 ($260/260 \, K, 377/377 \, Q\setminus K$).
- Under the identifying assumption of conditional exchangeability within joint strata, arm labels were permuted strictly within strata ($B = 10{,}000$ iterations). Importantly, this conditional exchangeability is a modeling assumption for stratified comparison, not an empirical proof that families or host genomes are mutually independent; unobserved phylogenetic or genome-sharing dependencies are not eliminated by permutation:

| Metric / Estimator | Tier 1 Primary Gate ($n = 353$ families) | Tier 2 3'-Stop Sensitivity ($n = 637$ families) |
|---|---|---|
| **Common Support Strata** | **16 / 16 strata** ($173/173 \, K$, $180/180 \, Q\setminus K$, 100% coverage) | **16 / 16 strata** ($260/260 \, K$, $377/377 \, Q\setminus K$, 100% coverage) |
| **MH-Weighted Excess Contrast $\Delta_{\mathrm{MH}}(E_{\mathrm{fam}})$** | **$+0.035358$** (Null SD $= 0.042346$, **$p = 0.4069$**) | **$+0.018391$** (Null SD $= 0.031223$, **$p = 0.5618$**) |
| **Pooled Common-Support Excess Contrast $\Delta_{\mathrm{pool}}(E_{\mathrm{fam}})$** | $+0.060839$ (Null SD $= 0.040917$, $p = 0.2244$) | $+0.037724$ (Null SD $= 0.029783$, $p = 0.2942$) |
| **MH-Weighted Raw Contrast $\Delta_{\mathrm{MH}}(C_{\mathrm{fam}})$** | **$+0.053850$** (Null SD $= 0.046674$, **$p = 0.2598$**) | **$+0.043859$** (Null SD $= 0.034337$, **$p = 0.2323$**) |
| **Pooled Common-Support Raw Contrast $\Delta_{\mathrm{pool}}(C_{\mathrm{fam}})$** | $+0.082402$ (Null SD $= 0.045100$, $p = 0.1589$) | $+0.068568$ (Null SD $= 0.032753$, $p = 0.1178$) |
| **MH-Weighted Continuous Excess $\Delta_{\mathrm{MH}}(E_{\mathrm{fam}}^{(R)})$** | **$+0.047337$** ($p = 0.2347$) | **$+0.032909$** ($p = 0.2514$) |

In both tiers, once stratified by species breadth, coarse host lineage, and flanking KO annotation opportunity, the between-arm excess recurrence difference is not statistically established ($p = 0.4069$ in Tier 1; $p = 0.5618$ in Tier 2).

---

### 6.2 Diagnostic Estimator: 1:1 Covariate Caliper Matching

We evaluated matching without replacement under both soft and strict caliper regimes to assess sensitivity to observational covariate distributions:
- **Soft-Caliper Matching ($148$ pairs):** Balances species breadth ($\text{SMD} = +0.0661$) and KO opportunity ($\text{SMD} = +0.0297$), reducing observational length imbalance from pre-match $\text{SMD} = +0.6439$ to post-match $\text{SMD} = +0.3622$ ($281.6$ aa vs. $243.8$ aa), with a $62.2\%$ exact phylum match rate. The matched between-arm contrast is null: $\Delta_{\mathrm{match}}(E_{\mathrm{fam}}) = +0.020781$ (sign-flip $p = 0.6628$) and $\Delta_{\mathrm{match}}(C_{\mathrm{fam}}) = +0.040541$ ($p = 0.5243$).
- **Strict Exact-Phylum + Tight-Length Sensitivity ($96$ pairs):** Enforces 100% exact GTDB-Tk R232 `phylum_group` matching ($96/96$ pairs same phylum group) and tight length calipers ($|\Delta z_{\mathrm{len}}| \le 0.50$, reducing post-match observational length $\text{SMD}$ to $+0.1242$ and KO opportunity $\text{SMD}$ to $+0.0013$). The matched contrast is identically zero on raw recurrence ($31.25\%$ vs. $31.25\%$, $\Delta = 0.0000, p = 1.0000$) and slightly negative on excess: $\Delta_{\mathrm{strict}}(E_{\mathrm{fam}}) = -0.005440$ ($p = 0.9226$).

---

### 6.3 Multivariable Regression Models ($n = 353$, Rank $= 9/9$, All Converged)

Multivariable regression models simultaneously evaluating the arm indicator (`is_k`) alongside observational covariates (species breadth, KO opportunity, focal protein length, mean pLDDT, and host phylum group) confirmed the stratified permutation results:

| Covariate / Parameter | OLS / HC3 on Excess $E_{\mathrm{fam}}$ ($\hat{\beta} \pm \mathrm{SE}_{\mathrm{HC3}}, p$) | OLS / HC3 on Raw $C_{\mathrm{fam}}$ ($\hat{\beta} \pm \mathrm{SE}_{\mathrm{HC3}}, p$) | MLE Logistic on Raw $C_{\mathrm{fam}}$ ($\hat{\beta} \pm \mathrm{SE}$, OR, $p$) | L2-Penalized Ridge Logistic ($\lambda = 0.5$) ($\hat{\beta} \pm \mathrm{SE}$, OR, $p$) |
|---|---|---|---|---|
| **Arm Indicator (`is_k`: $K$ vs. $Q \setminus K$)** | **$+0.029446 \pm 0.042950$ ($p = 0.4930$)** | **$+0.049661 \pm 0.046674$ ($p = 0.2873$)** | **$+0.280149 \pm 0.296352$ ($\text{OR} = 1.32, p = 0.3445$)** | **$+0.294897 \pm 0.278099$ ($\text{OR} = 1.34, p = 0.2890$)** |
| **Species Breadth ($S_F^{\mathrm{prim}}$)** | $-0.006722 \pm 0.020778$ ($p = 0.7463$) | $+0.066643 \pm 0.021244$ ($p = 0.0017$) | $+0.419208 \pm 0.141145$ ($\text{OR} = 1.52, p = 0.0030$) | $+0.396617 \pm 0.134989$ ($\text{OR} = 1.49, p = 0.0033$) |
| **KO Annotation Opportunity ($A_{\mathrm{flank}}^{\mathrm{prim}}$)** | **$+0.786956 \pm 0.080334$ ($p < 10^{-15}$)** | **$+0.942835 \pm 0.082188$ ($p < 10^{-15}$)** | **$+5.647450 \pm 0.744404$ ($p < 10^{-13}$)** | **$+4.503901 \pm 0.603709$ ($p < 10^{-13}$)** |
| **Sequence Length (per 100 aa)** | $-0.000142 \pm 0.013702$ ($p = 0.9917$) | $-0.003815 \pm 0.014596$ ($p = 0.7938$) | $-0.045775 \pm 0.123615$ ($\text{OR} = 0.96, p = 0.7112$) | $-0.057709 \pm 0.117950$ ($\text{OR} = 0.94, p = 0.6247$) |
| **Mean pLDDT (scaled per 0.10)** | $+0.037694 \pm 0.037404$ ($p = 0.3136$) | $+0.045500 \pm 0.040196$ ($p = 0.2577$) | $+0.236272 \pm 0.244049$ ($\text{OR} = 1.27, p = 0.3330$) | $+0.224332 \pm 0.231742$ ($\text{OR} = 1.25, p = 0.3330$) |
| **Lineage: `Actinomycetota` vs. `Thermoplasmatota`** | $+0.068105 \pm 0.065047$ ($p = 0.2951$) | $+0.086399 \pm 0.070133$ ($p = 0.2180$) | $+0.441274 \pm 0.387514$ ($\text{OR} = 1.55, p = 0.2548$) | $+0.480195 \pm 0.354689$ ($\text{OR} = 1.62, p = 0.1758$) |
| **Lineage: `multi` vs. `Thermoplasmatota`** | $-0.160060 \pm 0.056174$ ($p = 0.0044$) | $-0.198166 \pm 0.060164$ ($p = 0.0010$) | $-1.329295 \pm 0.402361$ ($\text{OR} = 0.26, p = 0.0010$) | $-1.146257 \pm 0.362903$ ($\text{OR} = 0.32, p = 0.0016$) |
| **Lineage: `other_bacteria_archaea` vs. `Thermoplasmatota`** | $+0.014012 \pm 0.065090$ ($p = 0.8296$) | $+0.017581 \pm 0.072484$ ($p = 0.8084$) | $+0.118737 \pm 0.385751$ ($\text{OR} = 1.13, p = 0.7582$) | $+0.124190 \pm 0.354055$ ($\text{OR} = 1.13, p = 0.7258$) |
| **Model Diagnostics** | `converged=True`, $R^2 = 0.2418$ | `converged=True`, $R^2 = 0.3044$ | `converged=True`, pseudo-$R^2 = 0.2708$ | `converged=True`, all $\mathrm{SE} \le 0.60$ |

---

## 7. Cross-Family Dependence Audit & Descriptive Species Resampling

1. **Locus Disjointness:** All $877$ Tier 1 eligible primary loci across the $353$ families reside on $869$ unique `(bin, contig)` pairs. Only $8$ contigs are shared by $\ge 2$ families, and there are **zero overlapping $\pm 2$-flank windows** across structural families.
2. **Host Taxon Overlap:** The $877$ loci span $617$ unique refined bins ($145$ shared by $\ge 2$ families) and $568$ unique dRep species clusters ($153$ shared by $\ge 2$ families).
3. **Descriptive Species-Cluster Bootstrap:** Resampling the $568$ dRep species clusters with replacement ($B = 2{,}000$) and weighting each family's pre-computed fixed outcomes by its mean species multiplicity yields:
   $$\mathrm{SE}_{\mathrm{boot}}(\Delta C_{\mathrm{fam}}) = 0.0353, \qquad \mathrm{SE}_{\mathrm{boot}}(\Delta E_{\mathrm{fam}}) = 0.0313$$
   As mandated, this bootstrap is reported strictly as a descriptive fixed-outcome sensitivity check; it does not recompute multi-species recurrence from resampled species pools or resolve cross-family phylogenetic dependence.
4. **Disjoint-Species Subset Sensitivity:** In the sub-cohort of $n = 112$ families with zero shared species across the entire Tier 1 cohort ($51 \, K, 61 \, Q \setminus K$), the unadjusted excess contrast is $\Delta(E_{\mathrm{fam}}) = +0.0634$ ($+0.2619$ in $K$ vs. $+0.1984$ in $Q \setminus K$), a descriptive sensitivity that does not establish an adjusted between-arm effect.

---

## 8. Flanking KO Annotations and Mobilome Keyword Screening Bounds

Across all $128$ Tier 1 conserved families ($70 \, K$, $58 \, Q \setminus K$):
- **Screening of Top Recurrent Flank (`best_ko_mobilome_flag`):** All $128$ selected top recurrent KO annotations were screened against the 11-keyword mobilome list (`transposase`, `integrase`, `phage`, `plasmid`, `conjugat`, `insertion`, `recombinase`, `mobilis`, `mobiliz`, `isla`, `prophage`). Exactly **`0 / 128` (`0.0%`)** matched this keyword list (`best_ko_mobilome_flag == 0`).
- **Strict Scope Bounds:** As noted in the review caveats, this keyword screen applies strictly to the top recurrent KO itself:
  - It does **not** prove the absence of mobile genetic elements, prophages, or plasmids across the genomic windows.
  - Across all $5{,}520$ Tier 1 flanks, $62.68\%$ ($3{,}460 / 5{,}520$) have no assigned KO (`mapped_noKO`) and many encode hypothetical proteins or uncharacterized domains.
  - Top recurrent KOs include diverse housekeeping, transport, and uncharacterized metabolic functions (e.g., $K03642$ `peptidase Do`, $K01142$ `exodeoxyribonuclease III`, $K03116$ `single-strand DNA-binding protein`, $K06888$ `DUF255`, and $K02221$ `YggT`).
  - No functional mechanism, operon co-transcription, or adaptive claims are made.

---

## 9. Programmatic Closure Verification

All 11 independent closure assertions in `scripts/03_analyse_conservation.py` passed (`true`):

1. `c1_all_tier1_eligible_loci_have_4_exact_mapped_flanks`: `true` ($5{,}520 / 5{,}520$ flanks in `{"mapped_KO", "mapped_noKO"}`).
2. `c1_all_ineligible_families_coded_minus_one`: `true` ($1{,}108$ Tier 1, $824$ Tier 2, $141$ Tier 3 families coded `-1`).
3. `c1_all_r_fam_within_0_1_bounds`: `true` ($0.0 \le R_{\mathrm{fam}} \le 1.0$ across all $1{,}461$ families in all 3 tiers).
4. `c1_numerator_denominator_species_closure`: `true` (`best_ko_species_count <= n_species_complete_window`).
5. `c2_zero_unassigned_r232_phylum_in_cohort`: `true` ($0$ unassigned phyla across all $1{,}461$ families).
6. `c2_all_binned_loci_have_positive_contig_len_bp`: `true` ($100\%$ of binned candidate rows have `contig_len_bp > 0`).
7. `c3_all_positional_null_draws_satisfy_primary_gate`: `true` (all $1{,}380$ Tier 1 loci have $M_c^{\mathrm{prim}} \ge 1$).
8. `c3_between_arm_stratified_common_support_verified`: `true` ($16 / 16$ common-support strata with $100\%$ coverage).
9. `c3_ols_hc3_and_logit_models_converged_full_rank_finite_se`: `true` (rank $9/9$, all models converged, finite SEs).
10. `c4_flank_category_totals_close_exactly`: `true` ($18{,}168 + 12{,}248 = 30{,}416$; $5{,}520 + 12{,}648 = 18{,}168$).
11. `c5_cohort_attrition_totals_close_exactly`: `true` ($353 + 1{,}108 = 1{,}461$; $637 + 824 = 1{,}461$; $1{,}320 + 141 = 1{,}461$).

---

## 10. Complete Artifact Manifest and Checksum Pinned Ledgers

All 13 primary output files generated by the full census are pinned with exact byte sizes and SHA-256 digests in `results/outputs_sha256.txt` (excluding self):

| Filename | File Size | SHA-256 Checksum |
|---|---:|---|
| `results/all_evaluable_families.tsv` | 339,429 B | `a5f6749765d76c2df3603394ac0c7539c45145868749eaec974740f2ba497fab` |
| `results/bin_contig_verification.tsv` | 1,933,978 B | `3c5c53f3f68cad8cfbf28baf68f323754c85574ee482fe29fc19d6f8301845fb` |
| `results/executed_inputs_manifest.tsv` | 3,113,314 B | `56fbc2b5e7b3c045e5a4b6491350ad78976d2b42ced43b9dfff3e2585d35408c` |
| `results/pilot_cohort.tsv` | 171,061 B | `c5abfa0ab93cc3bbff5c4484f0084ebf16dae0763a01d1012018420cda7f5363` |
| `results/pilot_common_support_strata.tsv` | 1,254 B | `bf6eef88b7a0adfb8e21c2653eb424a9272c59b9147493eca0e83e05d6a8dbd3` |
| `results/pilot_conservation_summary.json` | 58,870 B | `d002c26fa0bb66e51acb75fbbeed464a0b88a699a5b7a415e316843b3d9e26cb` |
| `results/pilot_contig_genes.tsv` | 31,719,940 B | `dc3e1e76a4472a84cf261793d5e66c4a0651d1bafb794c660f152b84b503b45d` |
| `results/pilot_family_conservation.tsv` | 321,446 B | `ac95c81d279151819b94a2631644c391844c2330225ab2f3fa5f395cd0b8a6c5` |
| `results/pilot_join_metrics.json` | 5,628 B | `0750690b369f99e0c1094f96af8fa06daccf51d7a3f8302569f46909ec3345fb` |
| `results/pilot_matched_pairs.tsv` | 20,995 B | `f399bf686255fdde97b571000a1d50fdf6d029a269edd725eec8d48a13f273a0` |
| `results/pilot_members_candidates.tsv` | 1,736,676 B | `992bf14bbcd4ecde8eee24242531bfa6607d269bd173997ca964d9c552239ddb` |
| `results/pilot_members.tsv` | 2,133,836 B | `16734c95dde15e901485b85f1af9f4c18acb4a28241f898581f782f03549f3f6` |
| `results/pilot_neighbourhoods.tsv` | 5,182,017 B | `6d2980b68005cebeaa57a3beebac233e43d1a577fb49aa695f441d8c8b22eca6` |
| `INPUTS.tsv` (Top-level working manifest) | 3,113,314 B | `56fbc2b5e7b3c045e5a4b6491350ad78976d2b42ced43b9dfff3e2585d35408c` |

Execution logs:
- `logs/full_53738443.out` (4,404 B, sha256: `66d0ef06143c720aeb401ba0333d690a2cf913346e4c27807d9b50b57e795c6c`)
- `logs/full_53738443.err` (0 B, clean empty standard error)

## Argon review bounds and artifact naming

The primary complete-window requirement excludes 1,108/1,461 families (75.84%): 173/570 K families (30.35%) and180/891 QnotK families (20.20%) remain. Tier2 retains260/570 K (45.61%) and377/891 QnotK (42.31%). These differential eligibility rates constrain generalization beyond the eligible families. Tier3 (1,320 families) is an unadjusted partial-window sensitivity with unequal flank completeness and species breadth; no Tier3 positional reference or adjusted contrast was calculated. Its raw p-value must not be used as evidence for a structural-group effect.

Cross-family sharing remains:877 primary loci among353 evaluable families span869bin-contig pairs (8shared),617bins (145shared) and568species (153shared). The species-disjoint sensitivity contains112families (51K/61QnotK), with mean excess recurrence .261855/.198430; it does not resolve all phylogenetic dependence. Phylum covariates use four coarse groups: Actinomycetota, Thermoplasmatota, cross-phylum families, and other bacteria/archaea; the last pools14phyla. The85bins without annotation files remain explicit missing inputs. The keyword screen cannot establish absence of mobile elements or assign function to focal dark proteins.

The legacy `pilot_*` filenames under `results/` contain the full1,461-family census in this round; they are retained for code/provenance compatibility. `scripts/run_full.sbatch` is the full workflow. Its duplicate `run_pilot.sbatch` is retained only as a historical alias and also invokes the full mode; use run_full.sbatch.
