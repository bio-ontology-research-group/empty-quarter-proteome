# DESIGN.md — Hypothesis 4 (V3): Within-Model Recovery of Confident Structural Regions Under Symmetric Reference-Conditional Cropped Search and Horvitz–Thompson Design Weighting

**Analysis Directory**: `reproducibility/darkmatter/round6_hypotheses/h4_structure`  
**Remote Workbench**: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure`  
**Date**: 2026-10-13  
**Status**: Revised Specification (V3 — Repaired Stratified Pilot with Horvitz–Thompson Finite-Population Estimators + Full $N = 5,851$ Region-Eligibility Census, Addressing `h4_review1` C1–C5 and `ROOT_AUDIT_NOTES.md`)  
**Superseded Versions**: `scripts/v1_superseded/`, `results/v1_superseded/`, `scripts/v2_superseded/`, `results/v2_superseded/`  

---

## 1. Scientific Context and Explicit Scope Narrowing

### 1.1 Context in Baseline Manuscript
In the Empty Quarter MAGs dark-proteome analysis (§1.3, lines 307–351), Foldseek structural searching (`foldseek/10-941cd33`, `easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 50`) of $2,838$ confidence-qualified representative candidate models ($\ge 100\text{ aa}$, mean $\text{pLDDT} \ge 0.70$; one representative per fold cluster) and $3,013$ confidence-qualified blinded PDB-derived positive controls (`KF...`) against PDB (`pdb_20260730/foldseek_pdb`) and AFDB (`afdb_uniprot_v6_uniprot_2025_03`) identified:
- **$1,031$ dual no-hit candidate models ($K = 1,031$, $36.33\%$)** under query-normalised TM score ($qTM_{pdb} < 0.5$ and $qTM_{afdb} < 0.5$).
- Only **$95$ candidate models ($3.35\%$)** remained dual no-hits under alignment-normalised TM score ($alnTM_{pdb} < 0.5$ and $alnTM_{afdb} < 0.5$).
- Positive controls showed an unfiltered dual no-hit rate of $3.35\%$ ($101 / 3,013$) under $qTM < 0.5$ and $2.26\%$ ($68 / 3,013$) under $alnTM < 0.5$.

### 1.2 Scope Narrowing from Original `round5/EXTENSIONS.md` and V1/V2
1. **Confident Structural Regions vs Validated Domains**: Contiguous residue segments passing outcome-blind per-residue confidence gates ($\ge 50\text{ aa}$, 15-aa rolling-window smoothed $\text{pLDDT} \ge 0.70$, bridging $\le 5\text{ aa}$ dips with gap mean $\text{pLDDT} \ge 0.55$) are operational **confident structural regions**, not validated evolutionary or functional domains.
2. **Model Uncertainty vs Intrinsic Disorder**: Low pLDDT reflects algorithmic prediction uncertainty and low predicted local distance difference test scores, **not established biological intrinsic disorder or flexible linkers**.
3. **Reference-Conditional Divergence vs Novel Folds**: Persistent low structural similarity ($qTM_{crop} < 0.5$) measures **reference-conditional divergence under heuristic Foldseek search against PDB**, not proof of a novel fold.
4. **What Cropped Search Measures**: Physical cropping and de novo Foldseek search of a confident region of length $L_{crop} \le L_{wc}$ tests whether an isolated confident segment achieves $qTM_{crop} \ge 0.5$ against PDB when normalised by $L_{crop}$ and aligned without flanking residues. Recovery reflects the combined effect of segment-length normalisation ($L_{crop}$ vs $L_{wc}$) and region-isolated 3Di/TMalign alignment, without asserting a single biophysical mechanism (such as inter-domain reorientation).
5. **Repairs from V2 (`h4_review1` C1–C5)**:
   - **C1**: Replaces post-hoc filtering of `--max-seqs 50` output (which saturated at 50 hits for close PDB homologs) with same-run `--max-seqs 1000` searches for both whole chains and cropped regions, explicitly separating `POST_EXCLUSION_EMPTY` from `PERSISTENT_LOW_SCORE` and tracking hit-list saturation (`total_hits_raw >= 1000`), while noting that a finite prefilter cap (`1000`) mitigates but does not eliminate homology-saturation censoring at stringent sequence-identity cutoffs (`< 20%`).
   - **C2**: Separates `NO_CONFIDENT_REGION` (`n_confident_regions == 0`, never searched as a crop) from `RAW_ZERO_RETRIEVAL` (`n_confident_regions >= 1` with 0 returned hits), and retains the exact full-cohort ($N = 5,851$) physical region-cropping census (`results/v3_full5851_cropped_regions_manifest.tsv`).
   - **C3**: Executes whole-chain and cropped-region PDB searches in the **same V3 run** with identical parameters and output fields (`fident`, `cigar`), applies Tiers 0–4 symmetrically to both whole-chain and cropped evaluations, explicitly partitions PDB-consistent vs AFDB-only strata (`K_local_pdb` vs `K_local_afdb_only`, `Q_not_K_pdb` vs `Q_not_K_afdb_only`), and tracks both match recovery (for whole-chain unmatched models) and match retention/loss (for whole-chain matched models).
   - **C4**: Documents the exact empirical 16-CPU scaling benchmark (Slurm job `53704591`) showing why full-cohort `--max-seqs 1000` TM-align across $5,851$ whole chains + $7,773$ cropped regions ($\sim 1.35 \times 10^7$ 3D superpositions, $\sim 183\text{ min}$) exceeds the $\le 16\text{ CPU}, \le 64\text{ GB}, \le 1\text{ h}$ single-job constraint, and implements explicit **Horvitz–Thompson / Hájek inverse-probability design-weighted finite-population estimators** (with finite-population-corrected standard errors and 95% CIs across the 7 sampling design strata, $\sum w_i = 5,851.00$) alongside the exact $N = 5,851$ region-eligibility census.
   - **C5**: Computes exact aligned residue pairs (`M` operations parsed from Foldseek `cigar` strings, verified equal to `qspan + tspan - alnlen` on 100% of rows) distinct from query coordinate span (`qend - qstart + 1`), uses the parent protein model as the primary unit of inference across all tiers, repairs the 3,013-control paired self-hit exclusion test ($b = 0, c = 42$, exact McNemar $p = 4.55 \times 10^{-13}$), and reports sequence identity strictly over retrieved hits.

---

## 2. Bounded Exploratory Hypothesis and Finite-Cohort Estimands

### 2.1 Bounded Hypothesis
Among confidence-qualified sequence-dark representative protein models ($\ge 100\text{ aa}$, mean $\text{pLDDT} \ge 0.70$; one representative per fold cluster) with low whole-chain PDB structural similarity ($qTM_{wc} < 0.5$), physical cropping and de novo Foldseek alignment of outcome-independent confident structural regions ($\ge 50\text{ aa}$, 15-aa rolling-window smoothed $\text{pLDDT} \ge 0.70$, bridging $\le 5\text{ aa}$ dips with mean $\text{pLDDT} \ge 0.55$) recovers reference-conditional structural matches ($qTM_{crop} \ge 0.5$) at a higher model-level rate in candidates with whole-chain local hits ($K_{local\_pdb}$, and overall $K_{local\_all}$) than in candidates lacking local hits ($K_{dual\_none}$), when evaluated symmetrically across sequence-identity exclusion tiers against internal fold-matched dark comparators ($Q \setminus K$) and PDB positive controls.

### 2.2 Primary Statistical Unit and Reference-Consistent Strata
- **Primary Statistical Unit**: The representative protein model $M$ (1 model per fold cluster; $N = 2,838$ dark candidates and $N = 3,013$ controls; total finite cohort $N = 5,851$). All multi-region outcomes are aggregated to the parent model $M$ before cross-stratum inference.
- **Exhaustive Reference-Consistent Stratification ($N = 5,851$ Population / $n = 300$ Stratified Pilot)**:
  Because baseline screen novelty $K$ was defined across two databases (PDB and AFDB) whereas cropped searches test recovery against PDB (`pdb_20260730/foldseek_pdb`), the finite cohort is partitioned into six mutually exclusive sub-strata (plus combined `K_local_all` and `Q_not_K_all`):
  1. **`K_local_pdb` ($N = 849$, pilot $n = 90$)**: Baseline $qTM_{pdb} < 0.5, qTM_{afdb} < 0.5$, and **$alnTM_{pdb} \ge 0.5$** (primary reference-matched local-hit stratum).
  2. **`K_local_afdb_only` ($N = 87$, pilot $n = 10$)**: Baseline $qTM_{pdb} < 0.5, qTM_{afdb} < 0.5, alnTM_{pdb} < 0.5$, and **$alnTM_{afdb} \ge 0.5$** (local hit present in AFDB only).
     - *Combined `K_local_all`*: $N = 849 + 87 = 936$ (pilot $n = 100$).
  3. **`K_dual_none` ($N = 95$, pilot $n = 50$)**: Baseline $qTM_{pdb} < 0.5, qTM_{afdb} < 0.5, alnTM_{pdb} < 0.5$, and $alnTM_{afdb} < 0.5$.
  4. **`Q_not_K_pdb` ($N = 1,107$, pilot $n = 55$)**: Baseline **$qTM_{pdb} \ge 0.5$** (reference-matched internal dark comparator in PDB).
  5. **`Q_not_K_afdb_only` ($N = 700$, pilot $n = 20$)**: Baseline $qTM_{pdb} < 0.5$ and **$qTM_{afdb} \ge 0.5$** (whole-chain matched in AFDB only).
     - *Combined `Q_not_K_all`*: $N = 1,107 + 700 = 1,807$ (pilot $n = 75$).
  6. **`Controls` ($N = 3,013$, pilot $n = 75$)**: Confidence-qualified blinded PDB-derived positive controls (`KF...`).

### 2.3 Explicit Stratified Sampling Design and Horvitz–Thompson Weights
The deterministic stratified pilot ($n = 300$, `PILOT_SEED = 20261013`) partitions the $N = 5,851$ finite population into $H = 7$ mutually exclusive design strata $h \in \{1, \dots, 7\}$ sampled without replacement:

| Design Stratum $h$ | Definition | Population $N_h$ | Sample $n_h$ | Inclusion Prob. $\pi_h = n_h / N_h$ | H–T Weight $w_h = N_h / n_h$ | FPC Factor $1 - n_h / N_h$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `K_local_len_100_200` | `K_local_all`, $100 \le L < 200\text{ aa}$ | 128 | 25 | 0.19531250 | 5.120000 | 0.80468750 |
| `K_local_len_200_300` | `K_local_all`, $200 \le L < 300\text{ aa}$ | 435 | 25 | 0.05747126 | 17.400000 | 0.94252874 |
| `K_local_len_300_400` | `K_local_all`, $300 \le L < 400\text{ aa}$ | 231 | 25 | 0.10822511 | 9.240000 | 0.89177489 |
| `K_local_len_400_10000` | `K_local_all`, $L \ge 400\text{ aa}$ | 142 | 25 | 0.17605634 | 5.680000 | 0.82394366 |
| `K_dual_none_SRS` | `K_dual_none`, simple random sample | 95 | 50 | 0.52631579 | 1.900000 | 0.47368421 |
| `Q_not_K_SRS` | `Q_not_K_all`, simple random sample | 1,807 | 75 | 0.04150526 | 24.093333 | 0.95849474 |
| `Controls_SRS` | `Controls`, simple random sample | 3,013 | 75 | 0.02489213 | 40.173333 | 0.97510787 |
| **Total** | **All 7 design strata** | **5,851** | **300** | — | **$\sum w_i = 5,851.00$** | — |

For any sub-population domain $D$ (with domain indicator $d_i \in \{0, 1\}$) and binary outcome $y_i \in \{0, 1\}$, the design-consistent Horvitz–Thompson / Hájek ratio estimator is:
$$\hat{p}_{D, HT} = \frac{\sum_{i \in s} w_i d_i y_i}{\sum_{i \in s} w_i d_i}, \qquad \widehat{\mathrm{Var}}(\hat{p}_{D, HT}) = \frac{1}{\left(\sum_{i \in s} w_i d_i\right)^2} \sum_{h=1}^7 N_h^2 \left(1 - \frac{n_h}{N_h}\right) \frac{s_{z,h}^2}{n_h}$$
where $z_i = d_i (y_i - \hat{p}_{D, HT})$ is the Taylor-linearized domain residual and $s_{z,h}^2$ is the sample variance of $z_i$ within design stratum $h$. When $D$ is a union of design strata (such as `K_local_all`, `K_dual_none`, `Q_not_K_all`, or `Controls`), $\sum_{i \in s} w_i d_i = N_D$ exactly and $\widehat{\mathrm{Var}}(\hat{p}_{D, HT})$ reduces to the exact unbiased stratified SRS variance with finite-population correction.

- **Boundary Taylor/Wald Intervals (`[0.0, 0.0]` / `[100.0, 100.0]`) and Empty-Domain Placeholders (`0.0` / `[0.0, 0.0]`)**:
  - In `ht_domain_estimate()` and `wilson_ci()`, an **empty domain** ($\sum_{i \in s} w_i d_i = 0$ or $n_D = 0$—occurring only for whole-chain-unmatched recovery in `Q_not_K_pdb`, where all $55 / 55$ models are whole-chain matched by definition, and whole-chain-matched retention in `K_dual_none`, where $0 / 50$ models are whole-chain matched) emits `0.0` with `[0.0, 0.0]` solely as an implementation placeholder to keep the JSON schema numeric; it is **not** an estimable zero rate.
  - At **boundary sample outcomes** ($k = 0$ or $k = n_D$ in a non-empty domain, such as $51 / 51 = 100.0\%$ region-eligible retention in `Q_not_K_pdb`, $0 / 300$ `POST_EXCLUSION_EMPTY`, or $0 / n_h$ pilot `NO_CONFIDENT_REGION` counts in `K_local_afdb_only`, `Q_not_K_afdb_only`, and `Controls`), all sample Taylor residuals satisfy $z_i = 0$, yielding `fpc_se_pct = 0.0` and degenerate Taylor/Wald intervals (`[0.0, 0.0]` or `[100.0, 100.0]`). These boundary Taylor/Wald intervals **do not establish finite-population absence or certainty**: exact finite-population `NO_CONFIDENT_REGION` rates are established directly by the $N = 5,851$ census (`v3_full5851_cropped_regions_manifest.tsv`), residual homology-saturation censoring is tracked via `total_hits_raw >= 1000`, and none of these boundary or empty-domain cells are used for the primary cross-stratum hypothesis inference.

### 2.4 Symmetric Exclusion Tiers (Applied Identically to Whole Chains and Crops)
For every query (both whole-chain $M$ and cropped region $r \in \mathcal{R}(M)$), hits returned by Foldseek (`--max-seqs 1000`) are filtered under five nested exclusion tiers:
- **Tier 0 (Unfiltered)**: All returned hits.
- **Tier 1 (Source-PDB Excluded)**: Exclude hits where `target[:4].lower() == pdb4` (`pdb4` is the 4-letter parent PDB code from `current_reference_control_manifest_20260802.tsv` for `Controls`; empty for dark candidates, so Tier 1 $\equiv$ Tier 0 for candidates).
- **Tier 2 (`SeqId < 30%`)**: Tier 1 hits with `fident < 0.30`.
- **Tier 3 (`SeqId < 25%`)**: Tier 1 hits with `fident < 0.25`.
- **Tier 4 (`SeqId < 20%`)**: Tier 1 hits with `fident < 0.20`.

### 2.5 Mutually Exclusive Outcome States and Transition Matrix
At each Tier $k \in \{0, 1, 2, 3, 4\}$:
1. **Whole-Chain Status (`wc_t{k}_status`)**:
   - `RAW_ZERO_RETRIEVAL`: `wc_total_hits_raw == 0`
   - `POST_EXCLUSION_EMPTY`: `wc_total_hits_raw > 0` and `wc_t{k}_remaining_hits == 0`
   - `LOW_SCORE`: `wc_t{k}_remaining_hits > 0` and `wc_t{k}_qtm < 0.5`
   - `MATCHED`: `wc_t{k}_remaining_hits > 0` and `wc_t{k}_qtm >= 0.5`
2. **Cropped-Region Status (`reg_t{k}_status`)**:
   - `RAW_ZERO_RETRIEVAL`, `POST_EXCLUSION_EMPTY`, `LOW_SCORE`, `MATCHED` (defined identically per region $r$).
3. **Model-Level Aggregated Crop Status (`crop_t{k}_status`)**:
   - `NO_CONFIDENT_REGION`: `n_confident_regions == 0`
   - `RAW_ZERO_RETRIEVAL`: `n_confident_regions >= 1` and all regions have `total_hits_raw == 0`
   - `POST_EXCLUSION_EMPTY`: `n_confident_regions >= 1`, $\ge 1$ region has `total_hits_raw > 0`, and all regions have `t{k}_remaining_hits == 0`
   - `PERSISTENT_LOW_SCORE`: `n_confident_regions >= 1`, $\ge 1$ region has `t{k}_remaining_hits > 0`, and $\max_r qTM_{crop, t_k}(r) < 0.5$
   - `MATCHED`: `n_confident_regions >= 1` and $\max_r qTM_{crop, t_k}(r) \ge 0.5$
4. **Symmetric Within-Model Transition (`trans_t{k}`, comparing `wc_t{k}_status` to `crop_t{k}_status` at the SAME Tier $k$)**:
   - If `wc_t{k}_status != "MATCHED"` (whole-chain unmatched at Tier $k$):
     - `WC_UNMATCHED_TO_CROP_RECOVERED` (`crop_t{k}_status == "MATCHED"`)
     - `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE` (`crop_t{k}_status == "PERSISTENT_LOW_SCORE"`)
     - `WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY` (`crop_t{k}_status == "POST_EXCLUSION_EMPTY"`)
     - `WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL` (`crop_t{k}_status == "RAW_ZERO_RETRIEVAL"`)
     - `WC_UNMATCHED_NO_CONFIDENT_REGION` (`crop_t{k}_status == "NO_CONFIDENT_REGION"`)
   - If `wc_t{k}_status == "MATCHED"` (whole-chain matched at Tier $k$):
     - `WC_MATCHED_RETAINED_IN_CROP` (`crop_t{k}_status == "MATCHED"`)
     - `WC_MATCHED_LOST_LOW_SCORE` (`crop_t{k}_status == "PERSISTENT_LOW_SCORE"`)
     - `WC_MATCHED_LOST_POST_EXCLUSION_EMPTY` (`crop_t{k}_status == "POST_EXCLUSION_EMPTY"`)
     - `WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL` (`crop_t{k}_status == "RAW_ZERO_RETRIEVAL"`)
     - `WC_MATCHED_LOST_NO_CONFIDENT_REGION` (`crop_t{k}_status == "NO_CONFIDENT_REGION"`)

### 2.6 Aligned-Pair Coverage vs Query Coordinate-Span Coverage
From Foldseek's 21-column output (`...qlen,tlen,alnlen,qstart,qend,tstart,tend,cigar`):
- **Exact Aligned Residue Pairs (`aligned_pairs`)**: Sum of `M` operation lengths parsed directly from `cigar`, asserted equal to `(qend - qstart + 1) + (tend - tstart + 1) - alnlen` on every hit row.
- **Aligned-Pair Query Coverage (`aligned_pair_qcov`)**: `aligned_pairs / qlen`.
- **Query Span (`qspan`) and Span Coverage (`qspan_cov`)**: `qspan = qend - qstart + 1` and `qspan_cov = qspan / qlen` (which includes query insertion residues aligned to target gaps inside `[qstart, qend]`).

---

## 3. Execution Architecture, Scaling Benchmark, and Verification Closure Checks

### 3.1 Empirical Full-Cohort Scaling Benchmark vs Repaired Stratified Pilot on IBEX
1. **Full-Cohort Cropping and Scaling Benchmark (`scripts/v3_full.sbatch`, Slurm Job `53704591` on `cn506-11-l`, 16 CPUs, 64 GB RAM)**:
   - Completed outcome-independent physical region cropping across all **$N = 5,851$ parent models**, producing **$7,773$ cropped confident region PDBs** across $5,690$ region-eligible models and identifying exactly **$161$ `NO_CONFIDENT_REGION` models** (retained with source and crop SHA-256 hashes in `results/v3_full5851_cropped_regions_manifest.tsv`).
   - In Step 3 (`foldseek easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 1000` on $5,851$ whole chains), `prefilter` passed a mean of $990$ (median $1,000$) target sequences per query ($5.79 \times 10^6$ pairs) in $146.5\text{ s}$ and `structurealign` completed in $32.9\text{ s}$, after which 3D `tmalign` progressed through $9 / 65$ ($13.85\%$) of whole-chain queries in $702\text{ s}$ ($11\text{m } 42\text{s}$). This establishes an empirical 16-CPU walltime of $84.5\text{ min}$ for whole-chain `tmalign` alone and $\sim 183\text{ min}$ ($\sim 3.05\text{ hours}$) for combined whole-chain ($5,851$ queries) + cropped-region ($7,773$ queries) searches at `--max-seqs 1000`—exceeding the $\le 16\text{ CPU}, \le 64\text{ GB}, \le 1\text{ h}$ allocation limit.
2. **Repaired Stratified Pilot Execution (`scripts/v3_pilot.sbatch`, Slurm Job `53706771` on `gpu510-32`, 16 CPUs, 64 GB RAM)**:
   - Executes `scripts/v3_crop_and_search_regions.py --mode pilot --threads 16 --max-seqs 1000` on the stratified pilot cohort ($n = 300$ parent models, $413$ cropped regions), running both whole-chain ($300$ queries) and cropped-region ($413$ queries) searches in the same run against `pdb_20260730/foldseek_pdb`, computing Horvitz–Thompson finite-population estimates with FPC standard errors, and integrating the exact $N = 5,851$ region-eligibility census from `v3_full5851_cropped_regions_manifest.tsv`.

### 3.2 Automated Code-Data Closure Assertions (`scripts/v3_verify_closure.py`)
1. **Full-Population & Pilot Partition Closure**: $849 + 87 + 95 + 1,107 + 700 + 3,013 = 5,851$ models in `v3_full5851_cropped_regions_manifest.tsv` ($7,773$ regions + $161$ zero-region models) and $90 + 10 + 50 + 55 + 20 + 75 = 300$ models in `v3_model_summaries.tsv` ($413$ regions + $11$ zero-region models; $\sum w_i = 5,851.00$).
2. **Category Sum & H–T Decomposition Closure**: For every stratum and Tier $k \in \{0..4\}$, the 5 `crop_t{k}_status` counts and 10 `trans_t{k}` transition counts sum to $n_s$, the 4 region status counts sum to $R_s$, and the 5 Horvitz–Thompson weighted status percentages sum to $100.0\%$.
3. **No-Region Exclusion Invariant**: Every model with `n_confident_regions == 0` has `crop_t{k}_status == "NO_CONFIDENT_REGION"` and is never counted in `RAW_ZERO_RETRIEVAL`.
4. **Tier Monotonicity & Symmetry Invariant**: For every query (whole-chain and region), remaining hit counts and max $qTM$ satisfy $\text{Tier } 0 \ge \text{Tier } 1 \ge \text{Tier } 2 \ge \text{Tier } 3 \ge \text{Tier } 4$, and for dark candidates (`pdb4 == ""`), $\text{Tier } 0 \equiv \text{Tier } 1$.
5. **Aligned-Pair Identity Invariant**: For every parsed hit row, `cigar_M_count == qspan + tspan - alnlen <= qspan <= qlen`.
6. **Input Manifest File Hashes vs Directory Locators and Reference Database Pinning**: In `INPUTS.tsv` (15 rows total), the 11 regular file entries carry exact byte sizes and 64-hex SHA-256 digests, whereas the 4 directory/database entries (`foldseek_pdb`, `afdb`, `query_pdb`, `control_folds`) are historical locator placeholders rather than file hashes (`afdb` is a historical upstream baseline database not searched in V3). The 32 actual component files of the searched PDB database (`pdb_20260730/foldseek_pdb`) are independently SHA-256 pinned in `reviews/darkmatter-2026-10-13/round6/root_checks/pdb_reference_verification.json` (Slurm job `53709618`; `4,530,919,251` bytes across the 32 actual components, excluding the `6,978`-byte self-entry of `foldseek_pdb.sha256` from the historical 33-entry total of `4,530,926,229` bytes), and all $5,851$ queried parent PDB files and $7,773$ cropped region PDB files are individually SHA-256 pinned in `results/v3_full5851_cropped_regions_manifest.tsv`.
