# DOSSIER: Hypothesis 4 (V3) — Within-Model Recovery of Confident Structural Regions Under Symmetric Reference-Conditional Cropped Search and Horvitz–Thompson Design Weighting

**Investigation Directory**: `reproducibility/darkmatter/round6_hypotheses/h4_structure`  
**Remote Workbench**: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure`  
**Execution Context**:
- **V3 Full-Cohort Cropping & 16-CPU Scaling Benchmark**: Slurm Job `53704591` on Ibex (`debug` partition, node `cn506-11-l`, 16 CPU cores, 64 GB RAM; physically cropped all $N = 5,851$ parent models into $7,773$ confident regions across $5,690$ region-eligible models and identified $161$ `NO_CONFIDENT_REGION` models in `results/v3_full5851_cropped_regions_manifest.tsv`, and benchmarked `foldseek easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 1000` on $5,851$ whole chains).
- **V3 Repaired Stratified Pilot Execution (`--max-seqs 1000`, Same-Run Whole-Chain + Cropped)**: Slurm Job `53706771` on Ibex (`debug` partition, node `gpu510-32`, 16 CPU cores, 64 GB RAM; whole-chain search on $n = 300$ models completed in $255.63\text{ s}$ [$115,916$ rows]; cropped-region search on $413$ regions completed in $160.59\text{ s}$ [$129,914$ rows]).
- **Full-Cohort Disconnected-Maxima & Control Leakage Audit**: Slurm Job `53700367` on Ibex (`cn506-11-l`, 16 CPU cores, 32 GB RAM, runtime $26\text{ s}$).  
**Harness Environment**: Python `3.10.20` (`metabolism-env`), Foldseek `10-941cd33` (`941cd33ff0771cd2e3f144e3293e22a2b87e9fda`).  
**Date**: 2026-10-13  
**Status**: Completed V3 Revision (Addressing `h4_review1` C1–C5 and `ROOT_AUDIT_NOTES.md`; V1 and V2 preserved in `scripts/v1_superseded/`, `results/v1_superseded/`, `scripts/v2_superseded/`, `results/v2_superseded/`)  

---

## 1. Executive Summary & Core Empirical Findings

In the Empty Quarter MAGs dark-proteome baseline screen (§1.3, lines 307–351), Foldseek structural searching (`foldseek/10-941cd33`, `easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 50`) of $2,838$ confidence-qualified representative candidate models ($\ge 100\text{ aa}$, mean $\text{pLDDT} \ge 0.70$; one model per fold cluster) and $3,013$ confidence-qualified blinded PDB positive controls (`KF...`) against PDB (`pdb_20260730/foldseek_pdb`) and AFDB (`afdb_uniprot_v6_uniprot_2025_03`) showed a large gap between query-normalised ($qTM < 0.5$: $1,031 / 2,838 = 36.33\%$ dual no-hits) and alignment-normalised ($alnTM < 0.5$: $95 / 2,838 = 3.35\%$ dual no-hits) TM scores.

To test whether confidence-qualified dark models with low whole-chain PDB structural similarity ($qTM_{wc} < 0.5$) contain outcome-independent confident structural regions ($\ge 50\text{ aa}$, 15-aa rolling-window smoothed $\text{pLDDT} \ge 0.70$, bridging $\le 5\text{ aa}$ dips with mean $\text{pLDDT} \ge 0.55$) that recover reference-conditional PDB matches ($qTM_{crop} \ge 0.5$) when physically cropped and searched de novo, we executed the V3 analysis addressing all five critical issues (`C1–C5`) from `h4_review1`:

1. **Disconnected Maxima ($N = 5,851$) & Paired Control Self-Hit Exclusion ($N = 3,013$)**:
   - In the baseline screen, `best_qtm` and `best_alntm` were maximized independently across targets (`05_integrate.py`). Across queries returning $\ge 2$ hits, $3,888 / 5,694$ ($68.28\%$) in PDB and $3,889 / 5,661$ ($68.70\%$) in AFDB—and **$2,313 / 2,816$ ($82.14\%$) of retrieved dark candidates**—had their top $qTM$ and top $alnTM$ hits on **different target structures** (`results/exact_hit_vs_maxima.tsv`).
   - Among the $3,013$ positive controls in the baseline PDB screen (`--max-seqs 50`), $2,844$ ($94.39\%$) retrieved their own 4-letter source PDB entry. Excluding exact 4-letter source PDB IDs increases the whole-chain PDB $qTM < 0.5$ control floor from **$117 / 3,013$ ($3.88\%$)** to **$159 / 3,013$ ($5.28\%$, Wilson 95% CI $[4.53\%, 6.13\%]$)** (`results/control_leakage_audit.tsv`). Because removing self-hits can only reduce or preserve a query's maximum $qTM$, the paired discordant counts are **$b = 0$ (Tier 0 unmatched $\to$ Tier 1 matched) and $c = 42$ (Tier 0 matched $\to$ Tier 1 unmatched)**, giving an exact two-sided McNemar binomial $p = 2 \times 0.5^{42} = 4.55 \times 10^{-13}$ (correcting V2's mislabeled unpaired $\chi^2$ $p = 0.010$).

2. **Full-Population ($N = 5,851$) Region-Eligibility Census & Separation of `NO_CONFIDENT_REGION` from `RAW_ZERO_RETRIEVAL` (`C2`)**:
   - Outcome-independent physical cropping across the entire $N = 5,851$ cohort (`results/v3_full5851_cropped_regions_manifest.tsv`, Slurm job `53704591`) generated **$7,773$ cropped confident region PDBs across $5,690$ region-eligible models** and identified **$161$ models ($2.75\%$) with zero eligible confident regions $\ge 50\text{ aa}$ (`NO_CONFIDENT_REGION`)**: $37 / 849$ ($4.36\%$) in `K_local_pdb`, $3 / 87$ ($3.45\%$) in `K_local_afdb_only` ($40 / 936 = 4.27\%$ in `K_local_all`), $7 / 95$ ($7.37\%$) in `K_dual_none`, $70 / 1,107$ ($6.32\%$) in `Q_not_K_pdb`, $31 / 700$ ($4.43\%$) in `Q_not_K_afdb_only` ($101 / 1,807 = 5.59\%$ in `Q_not_K_all`), and $13 / 3,013$ ($0.43\%$) in `Controls`.
   - In the $n = 300$ stratified pilot (`413` cropped regions from `289` region-eligible models + `11` `NO_CONFIDENT_REGION` models), separating `NO_CONFIDENT_REGION` (`3 K_local_pdb`, `4 K_dual_none`, `4 Q_not_K_pdb`) from Foldseek `RAW_ZERO_RETRIEVAL` (`1 K_local_pdb`, `10 K_dual_none`, `0 Q_not_K`, `1 Controls`) shows that **3 of the 4 unretrieved `K_local_pdb` models and 4 of the 14 unretrieved `K_dual_none` models were never searched as crops** because they lacked any contiguous confident region $\ge 50\text{ aa}$.

3. **Elimination of Control Hit-List Saturation Artifact via `--max-seqs 1000` (`C1`)**:
   - In V2 (`--max-seqs 50`), post-hoc sequence-identity filtering caused $17 / 103$, $22 / 103$, and $27 / 103$ control regions to become `POST_EXCLUSION_EMPTY` at Tiers 2 (`< 30%`), 3 (`< 25%`), and 4 (`< 20%`), artificially depressing unconditional control region match rates from $97.09\%$ to $80.58\%$, $75.73\%$, and $67.96\%$.
   - Re-running Foldseek with **`--max-seqs 1000`** on the identical $103$ control regions (Slurm job `53706771`) **reduced `POST_EXCLUSION_EMPTY` to $0 / 103$ ($0.00\%$) across all five tiers (Tiers 0–4)**, rescuing $2$, $19$, $23$, and $31$ control regions at Tiers 1, 2, 3, and 4 (`results/v3_recovery_summary.json`). Under `--max-seqs 1000`, unconditional control region match rates are **$102 / 103$ ($99.03\%$) at Tiers 1 and 2** and **$101 / 103$ ($98.06\%$) at Tiers 3 and 4** ($100.00\%$ and $99.02\%$ conditional on remaining hits). At the **parent-model level** ($n = 75$ controls), control match rates are **$74 / 75$ ($98.67\%$, H–T 95% CI $[96.09\%, 100.0\%]$) at Tiers 1–2** and **$73 / 75$ ($97.33\%$, H–T 95% CI $[93.71\%, 100.0\%]$) at Tiers 3–4**.

4. **Symmetric Same-Run Whole-Chain vs Cropped Recovery & Horvitz–Thompson Design-Weighted Estimates (`C3`, `C4`)**:
   - Searching both whole chains ($n = 300$) and cropped regions ($413$ regions) in the same V3 run (`--max-seqs 1000`, Tiers 0–4 applied symmetrically) across reference-consistent PDB vs AFDB-only strata establishes:
     - **`K_local_pdb` ($N = 849$, pilot $n = 90$, baseline $qTM_{pdb} < 0.5, alnTM_{pdb} \ge 0.5$)**:
       - Across all $90$ models, **$75 / 90$ ($83.33\%$ unweighted; **$83.72\%$ Horvitz–Thompson design-weighted**, FPC SE $4.28\%$, 95% CI $[75.33\%, 92.11\%]$)** achieve a cropped-region match ($qTM_{crop} \ge 0.5$) across Tiers 0–4 (`75 / 87 = 86.21%` unweighted / **`86.60%` H–T weighted** among region-eligible models).
       - Decomposing the $90$ `K_local_pdb` models against the same-run `--max-seqs 1000` whole-chain search: $9 / 90$ ($10.00\%$) already reach $qTM_{wc} \ge 0.5$ on the whole chain when `--max-seqs` is raised from $50$ to $1000$ ($8$ retained in crops, $1$ `NO_CONFIDENT_REGION`), while among the **$81 / 90$ models that remain whole-chain unmatched ($qTM_{wc} < 0.5$) at `--max-seqs 1000`**, **$67 / 81$ ($82.72\%$ unweighted; **$83.04\%$ H–T weighted**, FPC SE $4.58\%$, 95% CI $[74.06\%, 92.02\%]$; and $67 / 79 = 84.81\%$ unweighted / **$85.55\%$ H–T weighted** among region-eligible models)** recover a structural match ($qTM_{crop} \ge 0.5$) after physical cropping.
     - **`K_local_all` ($N = 936$, pilot $n = 100$, combining `K_local_pdb` [$n=90$] and `K_local_afdb_only` [$n=10$])**:
       - Overall cropped match rate is **$82 / 100$ ($82.00\%$ unweighted; **$82.15\%$ H–T weighted**, FPC SE $4.23\%$, 95% CI $[73.86\%, 90.45\%]$)** (`84.54%` unweighted / `84.65%` H–T weighted among region-eligible models), and same-run whole-chain-unmatched recovery is **$73 / 90$ ($81.11\%$ unweighted; **$81.33\%$ H–T weighted**, FPC SE $4.52\%$, 95% CI $[72.47\%, 90.19\%]$)**.
     - **`K_dual_none` ($N = 95$, pilot $n = 50$, baseline $qTM < 0.5, alnTM < 0.5$ in both DBs)**:
       - All $50 / 50$ models remain whole-chain unmatched at `--max-seqs 1000` ($44$ `LOW_SCORE`, $6$ `RAW_ZERO_RETRIEVAL`). Only **$18 / 50$ ($36.00\%$ unweighted and **$36.00\%$ H–T weighted**, FPC SE $4.72\%$, 95% CI $[26.75\%, 45.25\%]$; $18 / 46 = 39.13\%$ among region-eligible models)** recover a cropped match across Tiers 0–4, while $18 / 50$ ($36.00\%$) exhibit `PERSISTENT_LOW_SCORE`, $10 / 50$ ($20.00\%$) exhibit `RAW_ZERO_RETRIEVAL`, and $4 / 50$ ($8.00\%$) have `NO_CONFIDENT_REGION`.
     - **Model-Level Hypothesis Contrast (`K_local_pdb` vs `K_dual_none`)**:
       - Cropped recovery in `K_local_pdb` exceeds `K_dual_none` by **$+47.72$ percentage points H–T weighted** ($83.72\%$ vs $36.00\%$, FPC $z = 7.49$, two-sided $p = 6.91 \times 10^{-14}$; unweighted $75 / 90$ vs $18 / 50$, Fisher's exact two-sided $p = 2.35 \times 10^{-8}$) and by **$+47.04$ percentage points among models unmatched on the whole chain at `--max-seqs 1000`** ($83.04\%$ vs $36.00\%$ H–T weighted, FPC $z = 7.15$, $p = 8.53 \times 10^{-13}$; unweighted $67 / 81$ vs $18 / 50$, Fisher's exact $p = 7.46 \times 10^{-8}$).
     - **Reference-Consistent Internal Dark Comparator (`Q_not_K_pdb`, $N = 1,107$, pilot $n = 55$)**:
       - All $55 / 55$ models match on the whole chain at `--max-seqs 1000`; **$51 / 51$ ($100.00\%$) of region-eligible models retain their PDB match in cropped regions** (`WC_MATCHED_RETAINED_IN_CROP = 51`, `0` lost to low score), and the only $4 / 55$ ($7.27\%$) models losing their match upon cropping are the $4$ models with `NO_CONFIDENT_REGION`.

5. **Exact CIGAR Aligned-Pair Coverage vs Query Coordinate Span & Retrieved-Only Sequence Identity (`C5`)**:
   - Parsing exact `M` operations from Foldseek `cigar` strings across all $115,916$ whole-chain hit rows and $129,914$ cropped hit rows (`100%` verified equal to `qspan + tspan - alnlen`) shows that **$96.44\%$ ($111,788 / 115,916$) of whole-chain hits and $93.48\%$ ($121,448 / 129,914$) of cropped hits contain query insertions opposite target gaps (`qspan > aligned_pairs`)**. Across retrieved Tier 1 regions, true aligned-pair query coverage (`aligned_pairs / qlen`) is systematically lower than coordinate-span coverage (`qspan / qlen`): **`82.79%` vs `92.75%`** in `K_local_pdb`, **`74.86%` vs `88.77%`** in `K_dual_none`, **`89.53%` vs `96.14%`** in `Q_not_K_pdb`, and **`98.06%` vs `99.55%`** in `Controls`.
   - Computed strictly over retrieved regions (`t1_remaining_hits > 0`), mean top-hit sequence identity (`fident`) is **`6.17%` (max `24.40%`)** in `K_local_pdb` ($n=143$), **`5.18%` (max `21.20%`)** in `K_dual_none` ($n=45$), **`7.52%` (max `20.40%`)** in `Q_not_K_pdb` ($n=58$), and **`69.67%` (max `100.00%`)** in `Controls` ($n=102$). Four candidate regions have top-1 `t1_fident >= 20%` (`2` in `K_local_pdb`, `1` in `K_dual_none`, `1` in `Q_not_K_pdb`), and `417` raw candidate cropped hit rows (`114` whole-chain rows) have `fident >= 20%` (`25` cropped rows `>= 25%`, `3` cropped rows `>= 30%`). Even when those hits are removed at Tiers 2–4, every affected candidate region retains alternative hits with `qTM >= 0.5`, leaving candidate model and region match counts invariant across Tiers 0–4.

---

## 2. Empirical Scaling Benchmark, Stratified Sampling Design, and Full-Cohort Region Census

### 2.1 Empirical 16-CPU Scaling Benchmark (`Job 53704591` vs `Job 53706771`)
We first launched the full $N = 5,851$ cohort pipeline (`scripts/v3_full.sbatch`, Slurm Job `53704591` on `cn506-11-l`, 16 CPUs, 64 GB RAM, `01:00:00` limit; log retained at `logs/h4_v3_53704591.out`):
- **Step 2 (Full-Cohort Cropping, Completed)**: Parsed and cropped all $5,851$ parent models into $7,773$ confident region PDBs ($5,690$ eligible models + $161$ zero-region models) with SHA-256 digests in `results/v3_full5851_cropped_regions_manifest.tsv`.
- **Step 3 (`--max-seqs 1000` Whole-Chain Search Benchmark)**: Because `foldseek easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 1000` passes a mean of $990$ (median $1,000$) prefiltered sequences per query ($5.79 \times 10^6$ pairs vs $20,237$ pairs in V2's `--max-seqs 50` pilot), `prefilter` took $146.5\text{ s}$ and `structurealign` took $32.9\text{ s}$, after which 3D `tmalign` progressed through $9 / 65$ ($13.85\%$) of the $5,851$ whole-chain queries in $702\text{ s}$ ($11\text{m } 42\text{s}$). This corresponds to $84.5\text{ min}$ for whole-chain `tmalign` alone and **$\sim 183\text{ min}$ ($3.05\text{ hours}$) for combined whole-chain ($5,851$ queries) + cropped-region ($7,773$ queries) searches on 16 CPUs**, exceeding the $\le 16\text{ CPU}, \le 64\text{ GB}, \le 1\text{ h}$ single-job constraint.
- **Repaired Stratified Pilot Execution (`Job 53706771`)**: We therefore preserved `results/v3_full5851_cropped_regions_manifest.tsv` (giving exact finite-population region-eligibility counts for all $5,851$ models) and executed the repaired stratified pilot ($n = 300$ models, $413$ cropped regions, `--max-seqs 1000` on both whole chains and cropped regions) in Slurm Job `53706771` (`gpu510-32`, 16 CPUs, 64 GB RAM; `logs/h4_v3_53706771.out`), which completed whole-chain search in $255.63\text{ s}$ and cropped-region search in $160.59\text{ s}$ ($6\text{m } 56\text{s}$ combined search time).

### 2.2 Full-Cohort ($N = 5,851$) Region-Eligibility Census (`v3_full5851_cropped_regions_manifest.tsv`)

| Sub-Stratum / Combined Stratum | Finite Population Models ($N$) | Models With $\ge 1$ Confident Region | Models With $0$ Confident Regions (`NO_CONFIDENT_REGION`) | Census `% NO_CONFIDENT_REGION` | Total Cropped Regions ($R$) | Mean Regions per Model |Pilot Sample ($n$) | Pilot `NO_CONFIDENT_REGION` ($n_0$, Unw % / H–T %) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`K_local_pdb`** ($qTM_{pdb} < 0.5, alnTM_{pdb} \ge 0.5$) | 849 | 812 | 37 | **4.36%** | 1,269 | 1.495 | 90 | 3 (3.33% / 3.33%) |
| **`K_local_afdb_only`** ($alnTM_{pdb} < 0.5, alnTM_{afdb} \ge 0.5$) | 87 | 84 | 3 | **3.45%** | 130 | 1.494 | 10 | 0 (0.00% / 0.00%) |
| *Combined `K_local_all`* | *936* | *896* | *40* | ***4.27%*** | *1,399* | *1.495* | *100* | *3 (3.00% / 2.95%)* |
| **`K_dual_none`** (Dual $qTM < 0.5$ & $alnTM < 0.5$) | 95 | 88 | 7 | **7.37%** | 119 | 1.253 | 50 | 4 (8.00% / 8.00%) |
| **`Q_not_K_pdb`** (Baseline $qTM_{pdb} \ge 0.5$) | 1,107 | 1,037 | 70 | **6.32%** | 1,247 | 1.126 | 55 | 4 (7.27% / 7.27%) |
| **`Q_not_K_afdb_only`** ($qTM_{pdb} < 0.5, qTM_{afdb} \ge 0.5$) | 700 | 669 | 31 | **4.43%** | 883 | 1.261 | 20 | 0 (0.00% / 0.00%) |
| *Combined `Q_not_K_all`* | *1,807* | *1,706* | *101* | ***5.59%*** | *2,130* | *1.179* | *75* | *4 (5.33% / 5.33%)* |
| **`Controls`** (Blinded PDB controls `KF...`) | 3,013 | 3,000 | 13 | **0.43%** | 4,125 | 1.369 | 75 | 0 (0.00% / 0.00%) |
| **Total Finite Cohort** | **5,851** | **5,690** | **161** | **2.75%** | **7,773** | **1.328** | **300** | **11 (3.67% / 2.29% H–T)** |

### 2.3 Explicit Sampling Probabilities and Horvitz–Thompson Design Weights (`C4`)
In `K_local_all` ($N = 936$, pilot $n = 100$), equal allocation ($n_h = 25$) across four fixed length bins creates unequal inclusion probabilities $\pi_h = n_h / N_h$ and Horvitz–Thompson weights $w_h = N_h / n_h$:

| `K_local_all` Length Bin $h$ | Population $N_h$ (Share %) | Pilot Sample $n_h$ | Inclusion Prob. $\pi_h$ | H–T Weight $w_h$ | V3 Cropped `MATCHED` ($qTM_{crop} \ge 0.5$) | V3 `PERSISTENT_LOW_SCORE` | V3 `RAW_ZERO_RETRIEVAL` | V3 `NO_CONFIDENT_REGION` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`[100, 200)` aa** | 128 (13.68%) | 25 | 0.19531250 | 5.12 | 17 / 25 (**68.00%**) | 6 / 25 (24.00%) | 0 / 25 (0.00%) | 2 / 25 (8.00%) |
| **`[200, 300)` aa** | 435 (46.47%) | 25 | 0.05747126 | 17.40 | 20 / 25 (**80.00%**) | 4 / 25 (16.00%) | 0 / 25 (0.00%) | 1 / 25 (4.00%) |
| **`[300, 400)` aa** | 231 (24.68%) | 25 | 0.10822511 | 9.24 | 22 / 25 (**88.00%**) | 2 / 25 (8.00%) | 1 / 25 (4.00%) | 0 / 25 (0.00%) |
| **`[400, 10000)` aa** | 142 (15.17%) | 25 | 0.17605634 | 5.68 | 23 / 25 (**92.00%**) | 2 / 25 (8.00%) | 0 / 25 (0.00%) | 0 / 25 (0.00%) |
| **Unweighted Mean ($n=100$)** | — | 100 | — | — | **82.00%** (Wilson CI $[73.33, 88.30]$) | **14.00%** | **1.00%** | **3.00%** |
| **H–T Weighted ($N=936$)** | **936 (100.0%)** | — | — | $\sum w_i = 936.0$ | **82.15%** (SE $4.23\%$, CI $[73.86, 90.45]$) | **13.91%** (SE $3.85\%$, CI $[6.37, 21.44]$) | **0.99%** (SE $0.93\%$, CI $[0.00, 2.81]$) | **2.95%** (SE $1.93\%$, CI $[0.00, 6.73]$) |

---

## 3. Symmetric Same-Run Whole-Chain and Cropped Evaluation Results (`--max-seqs 1000`)

### 3.1 Model-Level Status Decomposition and Within-Model Transitions (Tier 1: Source-PDB Excluded)
For all candidate strata, Tier 0 $\equiv$ Tier 1 (and in V3 `--max-seqs 1000`, candidate model and region outcomes are completely invariant across Tiers 0, 1, 2, 3, and 4).

| Stratum ($N_{pop}$ / $n_{pilot}$) | V3 Whole-Chain Status at `--max-seqs 1000` (`MATCHED` / `LOW` / `EMPTY` / `RAW_ZERO`) | V3 Cropped Model Status at `--max-seqs 1000` (`MATCHED` / `PERS_LOW` / `EMPTY` / `RAW_ZERO` / `NO_REG`) | Unweighted Model Crop Match % (Wilson 95% CI) | H–T Weighted Model Crop Match % $\pm$ FPC SE (95% CI) | Same-Run `WC_UNMATCHED -> CROP_RECOVERED` ($k / n_{unm}$, Unw % / H–T % [95% CI]) | Region-Eligible Crop Match % (Unw % / H–T % [95% CI]) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`K_local_pdb`** ($849 / 90$) | `9 / 81 / 0 / 0` | `75 / 11 / 0 / 1 / 3` | **83.33%** ($[74.31, 89.63]$) | **83.72% $\pm$ 4.28%** ($[75.33, 92.11]$) | **67 / 81** (**82.72%** / **83.04%** $[74.06, 92.02]$) | **75 / 87** (**86.21%** / **86.60%** $[78.74, 94.46]$) |
| **`K_local_afdb_only`** ($87 / 10$) | `1 / 9 / 0 / 0` | `7 / 3 / 0 / 0 / 0` | **70.00%** ($[39.68, 89.22]$) | **69.75% $\pm$ 15.99%** ($[38.41, 100.0]$) | **6 / 9** (**66.67%** / **68.20%** $[35.35, 100.0]$) | **7 / 10** (**70.00%** / **69.75%** $[38.41, 100.0]$) |
| *Combined `K_local_all`* ($936 / 100$) | `10 / 90 / 0 / 0` | `82 / 14 / 0 / 1 / 3` | **82.00%** ($[73.33, 88.30]$) | **82.15% $\pm$ 4.23%** ($[73.86, 90.45]$) | **73 / 90** (**81.11%** / **81.33%** $[72.47, 90.19]$) | **82 / 97** (**84.54%** / **84.65%** $[76.73, 92.58]$) |
| **`K_dual_none`** ($95 / 50$) | `0 / 44 / 0 / 6` | `18 / 18 / 0 / 10 / 4` | **36.00%** ($[24.14, 49.86]$) | **36.00% $\pm$ 4.72%** ($[26.75, 45.25]$) | **18 / 50** (**36.00%** / **36.00%** $[26.75, 45.25]$) | **18 / 46** (**39.13%** / **39.13%** $[29.33, 48.94]$) |
| **`Q_not_K_pdb`** ($1,107 / 55$) | `55 / 0 / 0 / 0` | `51 / 0 / 0 / 0 / 4` | **92.73%** ($[82.74, 97.14]$) | **92.73% $\pm$ 3.45%** ($[85.96, 99.49]$) | $0 / 0$ (51/55 retained; 4 lost to `NO_REG`) | **51 / 51** (**100.00%** / **100.00%**) |
| **`Q_not_K_afdb_only`** ($700 / 20$) | `9 / 11 / 0 / 0` | `17 / 3 / 0 / 0 / 0` | **85.00%** ($[63.96, 94.76]$) | **85.00% $\pm$ 7.87%** ($[69.58, 100.0]$) | **8 / 11** (**72.73%** / **72.73%** $[46.79, 98.67]$) | **17 / 20** (**85.00%** / **85.00%** $[69.58, 100.0]$) |
| *Combined `Q_not_K_all`* ($1,807 / 75$) | `64 / 11 / 0 / 0` | `68 / 3 / 0 / 0 / 4` | **90.67%** ($[81.97, 95.41]$) | **90.67% $\pm$ 3.31%** ($[84.18, 97.16]$) | **8 / 11** (**72.73%** / **72.73%** $[46.79, 98.67]$) | **68 / 71** (**95.77%** / **95.77%** $[91.16, 100.0]$) |
| **`Controls` (Tier 1)** ($3,013 / 75$) | `73 / 1 / 0 / 1` | `74 / 0 / 0 / 1 / 0` | **98.67%** ($[92.83, 99.76]$) | **98.67% $\pm$ 1.32%** ($[96.09, 100.0]$) | **1 / 2** (50.00%; 73/73 WC matched retained) | **74 / 75** (**98.67%** / **98.67%** $[96.09, 100.0]$) |

> **Methodological Note on Boundary Taylor/Wald Intervals (`[0.0, 0.0]` / `[100.0, 100.0]`) and Empty-Domain Placeholders (`0.0` / `[0.0, 0.0]`)**: In `results/v3_recovery_summary.json` (`ht_domain_estimate` and `wilson_ci`), an **empty domain** ($\sum_{i \in s} w_i d_i = 0$ or $n_D = 0$—occurring only for whole-chain-unmatched recovery in `Q_not_K_pdb`, where all $55 / 55$ models are whole-chain matched by definition, and whole-chain-matched retention in `K_dual_none`, where $0 / 50$ models are whole-chain matched) emits `0.0` with `[0.0, 0.0]` solely as an implementation placeholder, not an estimable zero rate. At **boundary sample outcomes** ($k = 0$ or $k = n_D$ in a non-empty domain, such as $51 / 51 = 100.0\%$ region-eligible retention in `Q_not_K_pdb`, $0 / 300$ `POST_EXCLUSION_EMPTY`, or $0 / n_h$ pilot `NO_CONFIDENT_REGION` counts in `K_local_afdb_only`, `Q_not_K_afdb_only`, and `Controls`), Taylor-linearized Wald variance has sample residuals $z_i = 0$ and yields degenerate `[0.0, 0.0]` or `[100.0, 100.0]` intervals (`fpc_se_pct = 0.0`), which **do not establish finite-population absence or certainty**. Exact finite-population `NO_CONFIDENT_REGION` rates are established directly by the $N = 5,851$ census (`Table 2.2`), residual homology-saturation censoring is tracked via `total_hits_raw >= 1000` (`8 / 103` control regions), and none of these boundary or empty-domain cells are used for the primary cross-stratum inference.

**Exact 10-State Transition Breakdown (`trans_t1`)**:
- **`K_local_pdb` ($n = 90$)**: Among $81$ models with `wc_t1_status == "LOW_SCORE"`, $67$ transition to `WC_UNMATCHED_TO_CROP_RECOVERED`, $11$ to `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE`, $1$ to `WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL`, and $2$ to `WC_UNMATCHED_NO_CONFIDENT_REGION`. Among $9$ models that reached `wc_t1_status == "MATCHED"` when `--max-seqs` was raised from $50$ to $1000$, $8$ are `WC_MATCHED_RETAINED_IN_CROP` and $1$ (`41PRr3::41PRr3_702372_18`, `wc_t1_qtm = 0.5799`) is `WC_MATCHED_LOST_NO_CONFIDENT_REGION`.
- **`K_local_afdb_only` ($n = 10$)**: Although selected via baseline $alnTM_{afdb} \ge 0.5$ with $alnTM_{pdb} < 0.5$ (e.g. `V12Sr1::V12Sr1_267603_7`, which had $0$ hits in the historical PDB screen), $10 / 10$ return whole-chain PDB hits at `--max-seqs 1000` ($1$ `MATCHED` retained in crops; $9$ `LOW_SCORE`, of which $6$ transition to `WC_UNMATCHED_TO_CROP_RECOVERED` [including `V12Sr1::V12Sr1_267603_7` at `crop_t1_max_qtm = 0.8633`] and $3$ to `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE`).
- **`K_dual_none` ($n = 50$)**: All $50$ are whole-chain unmatched at `--max-seqs 1000` ($44$ `LOW_SCORE`, $6$ `RAW_ZERO_RETRIEVAL`); $18$ transition to `WC_UNMATCHED_TO_CROP_RECOVERED`, $18$ to `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE`, $10$ to `WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL`, and $4$ to `WC_UNMATCHED_NO_CONFIDENT_REGION`.
- **`Q_not_K_pdb` ($n = 55$)**: All $55$ are whole-chain `MATCHED` at `--max-seqs 1000`; $51 / 55$ ($92.73\%$) are `WC_MATCHED_RETAINED_IN_CROP` ($100.0\%$ of the $51$ region-eligible models, including `V34Dr2::V34Dr2_233315_1`, whose cropped region `reg1` increased from $qTM = 0.4161$ at `--max-seqs 50` in V2 to $qTM = 0.5089$ at `--max-seqs 1000` in V3), and $4 / 55$ ($7.27\%$) are `WC_MATCHED_LOST_NO_CONFIDENT_REGION`.
- **`Q_not_K_afdb_only` ($n = 20$)**: At `--max-seqs 1000`, $9 / 20$ match on the whole chain in PDB (`WC_MATCHED_RETAINED_IN_CROP = 9`), and among the $11$ whole-chain `LOW_SCORE` models, $8$ transition to `WC_UNMATCHED_TO_CROP_RECOVERED` and $3$ to `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE`.
- **`Controls` ($n = 75$)**: At Tiers 0–2, $73 / 75$ match on the whole chain and all $73 / 73$ ($100.0\%$) are `WC_MATCHED_RETAINED_IN_CROP`; of the $2$ whole-chain unmatched controls, $1$ (`LOW_SCORE` on whole chain) transitions to `WC_UNMATCHED_TO_CROP_RECOVERED` and $1$ (`RAW_ZERO_RETRIEVAL` on whole chain) transitions to `WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL`. At Tiers 3–4 (`SeqId < 25%` and `< 20%`), $1$ additional control becomes whole-chain `LOW_SCORE` and `WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE` (`73 / 75 = 97.33%` cropped model match rate).

### 3.2 Paired Comparison of V2 (`--max-seqs 50`) vs V3 (`--max-seqs 1000`) Across Multi-Tier Sequence-Identity Exclusions (`C1`)

#### A. Positive Control Models ($n = 75$) and Regions ($R = 103$) Across Tiers 0–4
| Exclusion Tier | V3 Control Whole-Chain Models ($n = 75$: `MATCH` / `LOW` / `EMPTY` / `ZERO`) | V3 Control Cropped Models ($n = 75$: `MATCH` / `PERS_LOW` / `EMPTY` / `ZERO` / `NO_REG`) | V3 Control Model Match % (Unw / H–T [95% CI]) | V2 (`--max-seqs 50`) Control Regions ($R = 103$: `MATCH` / `EMPTY`) | V3 (`--max-seqs 1000`) Control Regions ($R = 103$: `MATCH` / `LOW` / `EMPTY` / `ZERO`) | V3 Control Region Match % (Uncond / Cond on Remaining Hits) | Control Regions Rescued by `--max-seqs 1000` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Tier 0 (Unfiltered)** | `73 / 1 / 0 / 1` | `74 / 0 / 0 / 1 / 0` | **98.67%** ($[96.09, 100.0]$) | `100 / 0` (97.09%) | `102 / 0 / 0 / 1` | **99.03% / 100.00%** (102/102) | +2 |
| **Tier 1 (Source-PDB Excl.)** | `73 / 1 / 0 / 1` | `74 / 0 / 0 / 1 / 0` | **98.67%** ($[96.09, 100.0]$) | `100 / 0` (97.09%) | `102 / 0 / 0 / 1` | **99.03% / 100.00%** (102/102) | +2 |
| **Tier 2 (`SeqId < 30%`)** | `73 / 1 / 0 / 1` | `74 / 0 / 0 / 1 / 0` | **98.67%** ($[96.09, 100.0]$) | `83 / 17` (80.58%) | `102 / 0 / 0 / 1` | **99.03% / 100.00%** (102/102) | **+19** (17 empty + 2 low) |
| **Tier 3 (`SeqId < 25%`)** | `72 / 2 / 0 / 1` | `73 / 1 / 0 / 1 / 0` | **97.33%** ($[93.71, 100.0]$) | `78 / 22` (75.73%) | `101 / 1 / 0 / 1` | **98.06% / 99.02%** (101/102) | **+23** (22 empty + 1 low) |
| **Tier 4 (`SeqId < 20%`)** | `72 / 2 / 0 / 1` | `73 / 1 / 0 / 1 / 0` | **97.33%** ($[93.71, 100.0]$) | `70 / 27` (67.96%) | `101 / 1 / 0 / 1` | **98.06% / 99.02%** (101/102) | **+31** (27 empty + 4 low) |

- **Hit-List Saturation Diagnostic (`total_hits_raw >= 1000`)**: At `--max-seqs 1000`, $8 / 75$ ($10.67\%$) control whole-chain queries and $8 / 103$ ($7.77\%$) control cropped-region queries reach the $1,000$-hit cap (e.g., `KF000126_reg1` has `total_hits_raw = 1000`, retaining $36$ hits at Tier 2, $23$ at Tier 3, and $10$ at Tier 4 with `t4_qtm = 0.7040`), whereas $0 / 225$ candidate whole chains and $0 / 310$ candidate cropped regions reach $1,000$ hits. Although `--max-seqs 1000` reduced `POST_EXCLUSION_EMPTY` to $0$ across all $75$ control models and $103$ control regions in Tiers 1–4, $8$ control regions remain saturated at $1,000$ raw hits, meaning a finite prefilter cap mitigates but does not theoretically eliminate homology-saturation censoring in large protein families.

#### B. Region-Level Properties, Exact CIGAR Aligned-Pair Coverage vs Span Coverage, and Retrieved-Only `fident` (`C5`)
| Stratum | Total Regions ($R$) | Mean Length (aa) / Mean pLDDT | Raw Zero Retrieval ($R_0$) / Saturated ($\ge 1000$) | Retrieved Regions at Tier 1 ($R_{ret}$) | Tier 1–4 Region Match Count ($R_{match}$, Uncond % / Cond %) | Mean Aligned-Pair Coverage (`aligned_pairs / qlen`, Retrieved Only) | Mean Query Span Coverage (`qspan / qlen`, Retrieved Only) | Mean / Max Top-Hit `fident` (Retrieved Only) | Regions With Top-1 `t1_fident >= 20% / 25% / 30%` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`K_local_pdb`** | 150 | 118.62 / 0.8523 | 7 (4.67%) / 0 | 143 (95.33%) | **122 / 150** (**81.33% / 85.31%**) | **82.79%** | **92.75%** | **6.17% / 24.40%** | `2 / 0 / 0` |
| **`K_local_afdb_only`** | 15 | 118.27 / 0.8617 | 1 (6.67%) / 0 | 14 (93.33%) | **8 / 15** (**53.33% / 57.14%**) | **75.34%** | **85.40%** | **4.66% / 8.30%** | `0 / 0 / 0` |
| *Combined `K_local_all`* | *165* | *118.59 / 0.8532* | *8 (4.85%) / 0* | *157 (95.15%)* | ***130 / 165*** (***78.79% / 82.80%***) | ***82.12%*** | ***92.09%*** | ***6.04% / 24.40%*** | *`2 / 0 / 0`* |
| **`K_dual_none`** | 60 | 122.63 / 0.8473 | 15 (25.00%) / 0 | 45 (75.00%) | **22 / 60** (**36.67% / 48.89%**) | **74.86%** | **88.77%** | **5.18% / 21.20%** | `1 / 0 / 0` |
| **`Q_not_K_pdb`** | 58 | 132.47 / 0.8505 | 0 (0.00%) / 0 | 58 (100.00%) | **58 / 58** (**100.00% / 100.00%**) | **89.53%** | **96.14%** | **7.52% / 20.40%** | `1 / 0 / 0` |
| **`Q_not_K_afdb_only`** | 27 | 112.19 / 0.8342 | 1 (3.70%) / 0 | 26 (96.30%) | **22 / 27** (**81.48% / 84.62%**) | **83.47%** | **92.02%** | **4.65% / 9.90%** | `0 / 0 / 0` |
| *Combined `Q_not_K_all`* | *85* | *126.02 / 0.8453* | *1 (1.18%) / 0* | *84 (98.82%)* | ***80 / 85*** (***94.12% / 95.24%***) | ***87.65%*** | ***94.87%*** | ***6.63% / 20.40%*** | *`1 / 0 / 0`* |
| **`Controls` (Tier 1 / Tier 3–4)** | 103 | 194.85 / 0.8884 | 1 (0.97%) / 8 | 102 (99.03%) | **102 (99.03% / 100.0%)** T1–T2; **101 (98.06% / 99.02%)** T3–T4 | **98.06%** | **99.55%** | **69.67% / 100.0%** | `91 / 87 / 84` |

---

## 4. Methodological Boundaries and Interpretive Constraints

1. **Confident Structural Regions vs Validated Domains**: Contiguous segments passing $\ge 50\text{ aa}$ and smoothed $\text{pLDDT} \ge 0.70$ are operational *confident structural regions*, not validated evolutionary or functional domains. Low pLDDT measures model prediction uncertainty, **not established biological intrinsic disorder or flexible linkers**.
2. **What Physical Cropping Measures**: Cropping a region of length $L_{crop} \le L_{wc}$ and running de novo Foldseek search tests whether the isolated segment achieves $qTM_{crop} \ge 0.5$ when normalised by $L_{crop}$ and aligned without flanking residues. It reflects the combined algorithmic effect of segment-length normalisation and region-isolated 3Di/TM-align search, without proving a single biophysical mechanism (such as inter-domain linker reorientation).
3. **Reference-Conditional Divergence vs Novel Folds**: Persistent low scores ($qTM_{crop} < 0.5$, observed in $18 / 50 = 36.00\%$ of `K_dual_none` models and $11 / 90 = 12.22\%$ [$11.84\%$ H–T weighted] of `K_local_pdb` models) measure reference-conditional divergence under Foldseek search against `pdb_20260730/foldseek_pdb`, **not proof of novel protein folds**.
4. **Prefilter Censoring and Remote Superfamily Homology**: Increasing `--max-seqs` from $50$ to $1,000$ eliminated post-exclusion empty lists across all $75$ pilot controls ($103$ control regions) at `< 30%`, `< 25%`, and `< 20%` sequence identity, but $8 / 103$ control regions still saturated the $1,000$-hit cap and sequence-identity filters do not remove remote superfamily homologs (e.g. CATH/ECOD/SCOP cousins with `< 20%` sequence identity).
5. **Boundary Taylor/Wald Intervals (`[0,0]` / `[100,100]`) and Empty-Domain Placeholders (`0` / `[0,0]`)**: Taylor/Wald boundary `[0.0, 0.0]` or `[100.0, 100.0]` intervals (`fpc_se_pct = 0.0` when $k = 0$ or $k = n_D$ in a non-empty domain) do not establish finite-population absence or certainty, and empty-domain `0.0` / `[0.0, 0.0]` values (`w_den == 0` or $n_D = 0$, occurring only where a conditioning sub-domain has zero models by definition) are numeric schema implementation placeholders rather than estimable zeros; none of these cells are used for the primary cross-stratum inference.

---

## 5. Verification Receipts, Execution Logs, and Code-Data Closure

### 5.1 Pinned Remote Inputs (`INPUTS.tsv`) and Reference Database Component Verification
| Input Path | Access | Size (Bytes) | SHA-256 Digest |
| :--- | :---: | :---: | :--- |
| `per_rep_structures.tsv` | accessible | 1,436,582 | `9d2aa6f85f64ba1458a86ab60d144a58d35ad55163a1a56d6d177668906e8c5b` |
| `per_control_structures.tsv` | accessible | 250,513 | `33869ed27706b643b10874b1f54584b08f69a91c76744dfb1549ddb923a42996` |
| `current_reference_control_manifest_20260802.tsv` | accessible | 1,502,037 | `a65ba1b05cf3c398b29a976a96c39f89c0544f2fa354a90a792fb6d65f16e8b6` |
| `pfam_domtbl.txt` | accessible | 1,393,688 | `8f18b1750e83aae210b5896a7c107b91ea037cda7c00b8d53abe6b9d2a826c8f` |
| `per_protein_pfam.tsv` | accessible | 660,653 | `32b35974d5f5e012117cc938ab4068923e702b61d30448d47fef1bfbe6f968c2` |
| `pdb_hits.raw.tsv` | accessible | 65,056,711 | `02a49fe473cede0f01d7aa01d4f120bf752008f7c616773b555cc17570a871dd` |
| `afdb_hits.raw.tsv` | accessible | 120,552,103 | `c35f6b76ae092742807ce690a9d2661f187c3c511a536940c21f2b6bd90f97a9` |
| `foldseek_pdb.version` | accessible | 121 | `45bd6274a9fb5174bdf9f4f87e9ba66533006f07eae52581234db78752b99a28` |
| `05c_screen.sh` | accessible | 9,132 | `56b4ceff1f7296bb92397fbc4b398bd6be59b69089e8387185371ed0d162dca5` |
| `05_integrate.py` | accessible | 5,711 | `e1bd1dbc6483bd19f36b4de8c858fa58a7c3a8e0828a7ca3ffec169496ede8ad` |
| `foldseek_version.txt` | accessible | 41 | `7dcba988f3614b317e425ea8a241ea9efdb5e284d556dd47c65a0420b6494786` |

- **Regular-File SHA-256 Hashes vs Directory/Database Locator Entries in `INPUTS.tsv`**:
  - Of the 15 rows in `INPUTS.tsv`, the **11 regular file entries** listed above carry actual byte sizes and 64-hex SHA-256 digests. The remaining **4 entries** (`foldseek_pdb` [`directory_database`], `afdb` [`directory_database`], `query_pdb` [`20905_pdb_symlinks`, `directory`], and `control_folds` [`32_shards`, `directory`]) are **historical directory/database locator entries**, not file hashes.
  - **Searched PDB Reference Database (`pdb_20260730/foldseek_pdb`)**: All **32 actual Foldseek database component files** are independently SHA-256 verified and pinned in `reviews/darkmatter-2026-10-13/round6/root_checks/pdb_reference_verification.json` (Slurm job `53709618`; reference manifest `/ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/structural_calibration/references/pdb_20260730/foldseek_pdb.sha256`, SHA-256 `bd591f60d06786864a3ac0540557af56d7dfafd7bdf1f461b0f07ed4b567eaea`). The 32 actual database components sum to **`4,530,919,251` bytes**; the historical 33-entry manifest total of `4,530,926,229` bytes included the `6,978`-byte self-entry of `foldseek_pdb.sha256` (`4,530,919,251 + 6,978 = 4,530,926,229`), which was the only entry that mismatched itself while all 32 actual database components matched `100%`.
  - **Upstream AFDB Baseline (`afdb_uniprot_v6_uniprot_2025_03/database/afdb`) and Individual Model PDBs**: AFDB is a historical upstream baseline database (used only via the SHA-256-pinned baseline outputs `per_rep_structures.tsv`, `per_control_structures.tsv`, and `afdb_hits.raw.tsv`) and is **not searched in V3**. All $5,851$ queried parent PDB files and $7,773$ cropped region PDB files are individually SHA-256 pinned in `results/v3_full5851_cropped_regions_manifest.tsv` (and the $300$ pilot parent PDBs + $413$ cropped PDBs in `results/v3_cropped_regions_manifest.tsv`).

### 5.2 Generated V3 Output SHA-256 Checksums (`results/v3_results.sha256.tsv`)
| File Name | Location | Size (Bytes) | SHA-256 Digest | Description |
| :--- | :---: | :---: | :--- | :--- |
| `v3_full5851_cropped_regions_manifest.tsv` | Local & IBEX | 1,835,820 | `52cd3071948013edc6ca18d9a6d37d0d326c512d1a23db9f0e525e2fc6656eed` | Full $N = 5,851$ cohort region manifest ($7,773$ cropped regions + $161$ zero-region models; source & crop SHA-256 hashes) |
| `v3_cropped_regions_manifest.tsv` | Local & IBEX | 112,535 | `9b11053182562d2b1bc80bd3e1d1bdc9e2fc8c467b52002c28566ceece359272` | Stratified pilot ($n = 300$) region manifest ($413$ cropped regions + $11$ zero-region models; source & crop SHA-256 hashes) |
| `v3_wc_pdb_hits.raw.tsv` | IBEX (`results/`) | 26,585,879 | `517d3a4858e0f419bb67ca4e78fba957a789f74768995a3b248b1dfd39745191` | Raw V3 whole-chain Foldseek hits ($n = 300$ models, `--max-seqs 1000`, $115,916$ rows, 21 cols with `cigar`) |
| `v3_crop_pdb_hits.raw.tsv` | IBEX (`results/`) | 27,142,286 | `d86f09d371dab34b60dfea3c11f86b44672f980d1c00351f7661a3755e8b87df` | Raw V3 cropped-region Foldseek hits ($413$ regions, `--max-seqs 1000`, $129,914$ rows, 21 cols with `cigar`) |
| `v3_region_evaluations.tsv` | Local & IBEX | 188,351 | `e5e52b2ff24130a50669652b8c52c359b5eab8f8f8c3e941dfa199f0a27773bc` | Per-region evaluation table across Tiers 0–4 ($413$ regions, remaining hit counts, CIGAR aligned pairs vs span) |
| `v3_model_summaries.tsv` | Local & IBEX | 193,577 | `98f57db2416bdc3f00f034d6e0127b47fd6dd22f0e88b2a037ac5537a7ba4cc9` | Per-model summary table across Tiers 0–4 ($300$ models, symmetric whole-chain & cropped statuses, 10-state transitions, H–T weights) |
| `v3_recovery_summary.json` | Local & IBEX | 218,523 | `ff4ee8034ffd1c459a3673c7b7b8a045f5c80fae9e9c7a53584ddd12a88ac30d` | Complete JSON summary (unweighted & Horvitz–Thompson estimates, FPC SEs, full $N = 5,851$ region census, V2 vs V3 comparison) |
| `v3_audited_example_hits.tsv` | Local & IBEX | 45,232 | `ef282a768b2a00c4146065a5a89e30a968618a9119389dfefaa77f3e4f30f344` | $230$ representative raw hit rows (including $59$ top-1 region rows and $53$ top-1 whole-chain rows verified against summary tables) |

### 5.3 Reproducible Commands and Independent Closure Check
```bash
# 1. Execute V3 Stratified Pilot + Full 5,851 Region Census Integration on IBEX
sbatch /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure/scripts/v3_pilot.sbatch

# 2. Run Independent Local Code-Data Closure Verification
uv run python reproducibility/darkmatter/round6_hypotheses/h4_structure/scripts/v3_verify_closure.py
```
All 8 programmatic closure checks in `scripts/v3_verify_closure.py` pass with zero errors.
