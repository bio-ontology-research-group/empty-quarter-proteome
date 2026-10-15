# Evidence Dossier: Hypothesis 3 — Reference-Release Composition Under Matched Protein Amount and Length

**Date:** 2026-10-13  
**Author:** Gemini Flash High (for Robert Hoehndorf)  
**Assigned Local Directory:** `reproducibility/darkmatter/round6_hypotheses/h3_references`  
**Assigned Remote Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references`  
**Status:** Both Critical Dependencies Executed (Top-1 Bias Proven & Paired Public SMAG-vs-GEM Matched Search Completed); Defensible Falsifiable Design Established  

---

## 1. Executive Summary and Problem Statement

In the Empty Quarter metagenomic study, unannotated "dark" protein families exhibit a striking divergence in recovery between two landmark metagenome-assembled genome (MAG) catalogues of comparable total sequence volume:
- **SMAG (Ma et al. 2023, *Nat Commun* 14:7318):** 129,550,340 proteins from 40,039 global soil MAGs.
  - Dark family hit rate: **59.53%** (63,122 / 106,042); **40.47% absence**.
  - Control family hit rate: **98.43%** (104,382 / 106,042); **1.57% absence**.
- **GEM (Nayfach et al. 2021, *Nat Biotechnol* 39:499–509):** 122,978,191 proteins from 52,515 multi-biome MAGs.
  - Dark family hit rate: **35.62%** (37,767 / 106,042); **64.38% absence**.
  - Control family hit rate: **97.03%** (102,892 / 106,042); **2.97% absence**.

*(Note: Reviewer Argon previously inverted these values and fabricated an "8.2-fold" contrast; verified project artifacts confirm a 1.67× difference in hit rates, or +23.91 percentage points).*

The central scientific question of Hypothesis 3 is: **Does the higher dark-family recovery in SMAG reflect genuine biological/ecological representation of soil-associated microbial lineages, or is it an artifact of technical differences in database volume, protein length distributions, internal sequence redundancy, gene-caller conventions, or search heuristics?**

### Root Audit Corrections Incorporated in v2:
1. **No Calibration from Controls:** A control discordance ratio (e.g., 3.57× or 4.14×) does **not** establish how much technical bias explains the dark discordance ratio (10.15× or 12.19×). Controls differ fundamentally in length distribution (62.7% $<100$ aa in control queries vs 47.7% in dark queries), and conservation dynamics do not scale proportionally between universally conserved marker genes and unannotated dark families. All claims that controls "calibrate" or "explain away" specific proportions of technical bias have been removed.
2. **No "Residual Lineage Signal" Assertion:** We do not claim that residual asymmetry after partial adjustment proves phylogenetic lineage effects.
3. **No Association-as-Causation for Short Queries:** While short queries correlate with caller differences, we do not claim short-query effects are pure caller artifacts without direct cross-caller gene-level re-calling on the exact same contigs.
4. **Honest Caller-Conditional Estimand:** Because nucleotide source genomes for all 52,515 GEM MAGs are not staged on IBEX, gene-caller harmonization across the full catalogues is unachievable without massive external downloads. We explicitly define our estimand as **"reference-release composition under matched sequence count and length distribution, conditional on published caller releases"**, which is narrower than a causal biome claim or a fully caller-harmonized biological effect.
5. **Executed Critical Paired-Search Dependency:** Rather than stopping at permission barriers or relying on post-hoc filtering, we ingested public SMAG data from Zenodo (Record 8223844), performed de novo gene calling, constructed an exact size-matched GEM comparison, audited cross-catalogue overlap and redundancy, and executed the live paired DIAMOND blastp search under Slurm.

---

## 2. Statistical Estimand and Primary Unit

- **Primary Unit of Analysis:** The orthologous protein family representative $i \in \{1, \dots, N\}$ ($N = 106,042$ dark families from `FAMILIES_FINAL.txt`, paired with $N = 106,042$ length- and prevalence-matched annotated control families from `CONTROL_FAMILIES.txt`).
- **Paired Binary Outcome:** For each query $i$ evaluated against database arm $A$ (SMAG release subset) and database arm $B$ (GEM matched subset):
  $$(Y_{i,A}, Y_{i,B}) \in \{0, 1\}^2$$
  where $Y_{i,k} = 1$ if query $i$ aligns to $\ge 1$ reference protein with $E \le 10^{-3}$ and $\ge 30\%$ sequence identity under `diamond blastp --sensitive`, and $0$ otherwise.
- **Primary Estimand (Discordance Difference):**
  $$\delta_{\text{disc}} = P(Y_{i,A}=1, Y_{i,B}=0) - P(Y_{i,A}=0, Y_{i,B}=1) = \frac{b - c}{N}$$
  where $b$ is the number of arm-$A$-only hits and $c$ is the number of arm-$B$-only hits.
- **Hypothesis Test:** McNemar's exact paired test for marginal homogeneity:
  $$\chi^2 = \frac{(|b - c| - 1)^2}{b + c} \sim \chi^2_1$$
  with discordant odds ratio $\psi = b / c$.
- **Length Stratification:** Evaluated across five predefined protein length strata ($<100$, $100–149$, $150–249$, $250–399$, and $\ge 400$ aa) to inspect whether discordant recovery is uniform or length-dependent.
- **Convenience Conditioning Disclosure:** For the pilot matched search, the SMAG subset was derived from the first two ordered archive parts (`mag.tar.gz.aa` and `mag.tar.gz.ab`) of Zenodo Record 8223844. We openly condition on this archive-order convenience subset; it proves pipeline, gene calling, and paired re-search feasibility, but does not constitute an unconditioned random sample of all global soils.

---

## 3. Falsifiable Hypotheses

- **Null Hypothesis ($H_0$):**  
  Under equal sequence count, comparable residue volume, and identical search parameters, the marginal probability of recovering an Empty Quarter dark protein family in the SMAG release subset is equal to that in the GEM matched subset:
  $$H_0: P(Y_{\text{SMAG}}=1, Y_{\text{GEM}}=0) = P(Y_{\text{SMAG}}=0, Y_{\text{GEM}}=1) \quad (\psi = 1)$$
- **Alternative Hypothesis ($H_1$):**  
  The recovery rates differ systematically between the two reference releases under matched sequence count:
  $$H_1: P(Y_{\text{SMAG}}=1, Y_{\text{GEM}}=0) \ne P(Y_{\text{SMAG}}=0, Y_{\text{GEM}}=1) \quad (\psi \ne 1)$$

---

## 4. Key Confounds and Controlled Evidence

| Confound | Historical State | Controlled Treatment in H3 | Residual Limitation |
|---|---|---|---|
| **Database Size** | SMAG 129.6M vs GEM 123.0M (1.053×) | Exact sequence matching in pilot (693,227 vs 694,598; 0.20% diff) | Evaluated on pilot slice; full 120M scale requires full download |
| **Top-1 Reporting Bias** | Both searched with `-k 1` | Proved post-hoc filtering causes 87.3% undercount; executed fresh physical re-search | None for pilot; must be maintained in any full extension |
| **Sequence Length** | Uncontrolled | Length-stratified across 5 bins ($<100$ to $\ge 400$ aa) | Biological length variation within bins remains |
| **Internal Redundancy** | Unknown in raw hits | Audited exact 32-mer MD5 sequence hashes: SMAG 0.15% vs GEM 2.26% | Species-level dereplication requires full CheckM2/ANI clustering |
| **Cross-Catalogue Overlap**| Unknown | Audited exact identical sequences: 0.00% overlap | Homologous (non-identical) overlap requires all-vs-all diamond |
| **Gene Caller** | SMAG: Pyrodigal 3.7.1; GEM: Prodigal 2.6.3 | Recalled SMAG with Pyrodigal 3.7.1; GEM conditioned on published release | Estimand is explicitly caller-conditional, not fully harmonized |
| **Sampling Structure** | Global soil vs Multi-biome | Tested on 2,000 deterministic queries (1,000 dark, 1,000 ctrl) | SMAG pilot uses archive-order convenience subset |

---

## 5. Empirical Feasibility Evidence

### 5.1 Proof of the Top-1 Subsampling Fallacy (Job 53700569)

To determine whether matched-subset comparisons could be performed post hoc by filtering the existing top-1 hit tables (`smag_hits.m8` and `gem_hits.m8`), we executed an empirical proof on IBEX (node `cn506-11-l`, 8 CPUs, 145 seconds):
- **Subsample:** 2,500 GEM genomes randomly sampled (seed `20261013`; 5,815,177 proteins, 1.835 Gb residues; ~5% of GEM).
- **Queries:** 2,000 deterministic queries (1,000 dark, 1,000 control; 200 per length bin; seed `20261013`).
- **Search:** Fresh DIAMOND blastp (`--sensitive -e 1e-3 -k 1`) vs post-hoc filtering of full GEM top-1 hit table.

**Results (`results/pilot_top1_bias_summary.tsv`):**
- **Dark Families:** 
  - True empirical hits: **126 / 1,000 (12.60%)**
  - Post-hoc filtered hits: **16 / 1,000 (1.60%)**
  - Masked valid hits: **110**
  - **Post-Hoc Undercount Rate:** $\mathbf{87.30\%}$
- **Control Families:**
  - True empirical hits: **953 / 1,000 (95.30%)**
  - Post-hoc filtered hits: **58 / 1,000 (5.80%)**
  - Masked valid hits: **895**
  - **Post-Hoc Undercount Rate:** $\mathbf{93.91\%}$

**Conclusion:** When a reference database is subsampled, top-1 search tables mask the overwhelming majority of valid homologous alignments because secondary hits with $E \le 10^{-3}$ were suppressed by the `-k 1` cap during the full search. **Post-hoc subsampling is fatally biased and invalid; physical re-search is non-negotiable.**

---

### 5.2 Public Source Ingestion and Verification

To bypass local permission locks on `/ibex/user/velazqal/eqmags_31/smag/` without privilege escalation, we ingested the public SMAG archive from Zenodo Record 8223844:
- API metadata verified: `reviews/darkmatter-2026-10-13/round6/public_source_checks/zenodo_8223844.json`.
- Archive parts downloaded and verified with official MD5 checksums (`results/smag_zenodo_download_manifest.tsv`):
  - `mag.tar.gz.aa`: 104,857,600 bytes (MD5: `5fb8a52b97940d1d1c14aecda57ac5bc`) — Verified.
  - `mag.tar.gz.ab`: 104,857,600 bytes (MD5: `de5f1ecb60b8ceac42b616e2780b2d04`) — Verified.
- Genome extraction: Extracted into project scratch space. Inspected all records; confirmed 212 complete, intact MAGs (0 truncated trailing records).
- De novo gene calling: Executed `Pyrodigal v3.7.1 -p meta --no-stop-codon` across 8 parallel worker processes (`scripts/05_extract_and_genecall_smag.py`).
  - Total SMAG proteins generated: **693,227** across 212 MAGs (202,894,612 residues).
  - Runtime: 203.8 seconds (1.04 MAGs/sec).
  - Manifest written: `results/smag_pilot_genome_manifest.tsv`.

---

### 5.3 Redundancy, Length, and Overlap Audit

Using `scripts/06_match_gem_and_search.py`, we sampled 292 GEM genomes (seed `20261013`) to achieve exact protein count parity with the 212 SMAG MAGs. We then performed an exact sequence-level audit (`results/smag_vs_gem_redundancy_overlap_audit.tsv`):

| Metric | SMAG Pilot Subset | GEM Matched Subset | Contrast / Notes |
|---|---|---|---|
| Genomes Sampled | 212 MAGs | 292 MAGs | Both from published MAG collections |
| Total Proteins | **693,227** | **694,598** | Difference: 1,371 seqs (**0.20% match**) |
| Total Residues | 202,894,612 | 221,718,516 | GEM has +9.28% residue volume |
| Mean Protein Length | 292.68 aa | 319.20 aa | GEM proteins are slightly longer on average |
| Median Protein Length | 252.00 aa | 271.00 aa | Median length difference: +19.0 aa |
| Unique Sequence Hashes | 692,194 | 678,868 | Exact 32-mer MD5 protein hashes |
| **Internal Sequence Redundancy** | **0.15%** (1,033 duplicates) | **2.26%** (15,730 duplicates) | GEM has 15× higher internal sequence redundancy |
| **Cross-Catalogue Shared Sequences** | **0** | **0** | **0.00% exact overlap** between the two subsets |

---

### 5.4 Paired Matched Search Results (Job 53704317)

Both subsets were indexed with `diamond makedb` and searched against the 2,000 deterministic queries (1,000 dark, 1,000 control) under identical settings (`--sensitive -e 1e-3 -k 1 -p 8`). Results are stored in `results/smag_vs_gem_matched_search_results.tsv`:

| Cohort | $N$ Queries | Both Hit | SMAG Only ($b$) | GEM Only ($c$) | Neither (Absent) | SMAG Hit Rate | GEM Hit Rate | McNemar $\chi^2$ | $p$-value | Discordant Ratio ($b/c$) |
|---|---|---|---|---|---|---|---|---|---|---|
| **Dark** | 1,000 | 24 | **132** | **13** | 831 | **15.60%** | **3.70%** | **96.03** | $\mathbf{1.13 \times 10^{-22}}$ | **10.15×** |
| **Control** | 1,000 | 875 | 58 | 14 | 53 | **93.30%** | **88.90%** | **25.68** | $\mathbf{4.03 \times 10^{-7}}$ | **4.14×** |
| **All Combined** | 2,000 | 899 | 190 | 27 | 884 | **54.45%** | **46.30%** | **120.94** | $\mathbf{3.94 \times 10^{-28}}$ | **7.04×** |

#### Length-Stratified Matched Search Breakdown (`results/smag_vs_gem_matched_length_stratified.tsv`):

| Cohort | Length Bin (aa) | $N$ | Both Hit | SMAG Only ($b$) | GEM Only ($c$) | Neither | SMAG Hit Rate (%) | GEM Hit Rate (%) | Discordant Ratio ($b/c$) |
|---|---|---|---|---|---|---|---|---|---|
| **Dark** | $<100$ | 200 | 10 | **27** | 6 | 157 | 18.50% | 8.00% | **4.50×** |
| **Dark** | $100–149$ | 200 | 2 | **24** | 0 | 174 | 13.00% | 1.00% | **$\infty$ (24 vs 0)** |
| **Dark** | $150–249$ | 200 | 6 | **25** | 1 | 168 | 15.50% | 3.50% | **25.00×** |
| **Dark** | $250–399$ | 200 | 3 | **27** | 2 | 168 | 15.00% | 2.50% | **13.50×** |
| **Dark** | $\ge 400$ | 200 | 3 | **29** | 4 | 164 | 16.00% | 3.50% | **7.25×** |
| **Control**| $<100$ | 200 | 136 | 26 | 7 | 31 | 81.00% | 71.50% | **3.71×** |
| **Control**| $100–149$ | 200 | 162 | 20 | 2 | 16 | 91.00% | 82.00% | **10.00×** |
| **Control**| $150–249$ | 200 | 180 | 9 | 5 | 6 | 94.50% | 92.50% | **1.80×** |
| **Control**| $250–399$ | 200 | 197 | 3 | 0 | 0 | 100.00% | 98.50% | **$\infty$ (3 vs 0)** |
| **Control**| $\ge 400$ | 200 | 200 | 0 | 0 | 0 | 100.00% | 100.00% | Equal (0 vs 0) |

---

## 6. Synthesis and Scientific Interpretation

1. **Subsampling Parity Does Not Eliminate Dark Disparity:**  
   Even when the reference database sequence count is matched to within 0.20% (693k vs 694k proteins) and GEM actually possesses slightly more total residues (+9.28%) and longer proteins (+26.5 aa mean), SMAG recovers **4.2 times more total dark hits** (156 vs 37) and exhibits a **10.15× asymmetry** in discordant hits (132 vs 13, McNemar $\chi^2 = 96.03, p = 1.13 \times 10^{-22}$).
2. **Robustness Across Length Strata:**  
   The SMAG recovery advantage is not an artifact restricted to tiny edge-truncated peptides. In full-length proteins ($\ge 400$ aa), SMAG recovers 29 unique families compared to 4 for GEM (7.25× discordant ratio). Across all length bins $\ge 100$ aa, SMAG-only dark hits exceed GEM-only hits by $7.25\times$ to $\infty$.
3. **Control Cohort Sensitivity Behavior:**  
   In annotated controls, recovery reaches 100.00% in both databases for long proteins ($\ge 400$ aa, 200/200 in both), demonstrating that both catalogues have full sensitivity for universally conserved orthologs. Discordance in controls is concentrated in the $<150$ aa range, consistent with shorter alignments being more sensitive to edge heuristics and database representation.
4. **Honest Scope of Inference:**  
   This pilot proves the feasibility of downloading public SMAG parts, verifying integrity, calling genes, and executing fair paired searches without privilege escalation. However, because this pilot utilized an archive-order convenience subset of 212 SMAG MAGs and published GEM Prodigal annotations, it does not support an unconditioned claim of causal soil adaptation or geographic endemism. The result must be reported strictly as **"differential reference composition between published SMAG and GEM releases under matched protein sequence volume"**.

---

## 7. Full-Test Recipe

To scale this verified workflow to a definitive manuscript-grade test across the full 106,042 universe:

### Arm Design (3-Arm Controlled Comparison):
1. **Arm 1: SMAG Soil (18,000 MAGs)**  
   - Download remaining Zenodo 8223844 parts on IBEX (39 GB total).
   - Randomly sample 18,000 MAGs across the 7,725 species clusters.
   - Call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon).
2. **Arm 2: GEM Soil (18,000 MAGs)**  
   - Filter the 52,515 GEM MAGs using Nayfach 2021 metadata (`ecosystem == "Soil"`), yielding ~18,000 MAGs.
   - Extract published protein FASTA.
3. **Arm 3: GEM Non-Soil Multi-Biome (18,000 MAGs)**  
   - Randomly sample 18,000 non-soil MAGs from GEM (marine, aquatic, wastewater, host-associated).
   - Extract published protein FASTA.

### Computational Protocol:
- **Residue Volume Harmonization:** Subsample genomes within each arm so that total residue counts match within $\pm 1.0\%$.
- **Dereplication Audit:** Run `mmseqs easy-cluster` at 95% and 50% sequence identity to measure effective sequence diversity per arm.
- **Search Execution:** Search all 212,084 queries (`query_all.faa`) against each arm using identical parameters:
  `diamond blastp -q query_all.faa -d {arm.dmnd} --sensitive -e 1e-3 -k 1 -p 16`
- **Resource Footprint:** 
  - Zenodo download: ~39 GB disk, ~2 hours transfer.
  - Pyrodigal calling (18k MAGs): ~16 CPUs, ~4 hours walltime.
  - DIAMOND search: 3 runs × ~30 minutes on IBEX `batch` partition. Total compute: ~120 CPU-hours.

---

## 8. Informative Outcomes and Go/No-Go Boundaries

### Informative Positive Outcome:
If Arm 1 (SMAG Soil) and Arm 2 (GEM Soil) yield comparable dark hit rates that both significantly exceed Arm 3 (GEM Non-Soil) ($\Delta_{\text{Soil vs Non-Soil}} > 0, p < 10^{-15}$ across all length strata):
- **Conclusion:** Confirms that soil-biome reference content specifically enhances recovery of Empty Quarter dark families, independent of catalogue release or gene caller.

### Informative Negative Outcome:
If Arm 2 (GEM Soil) and Arm 3 (GEM Non-Soil) yield identical dark hit rates, and the excess recovery is confined entirely to Arm 1 (SMAG Soil vs GEM Soil):
- **Conclusion:** Proves that the historical "SMAG advantage" is **not a biological soil-biome effect**, but an artifact of catalogue-specific properties (e.g. Zhejiang pipeline assembly heuristics, contig filtering thresholds, or Pyrodigal calling density).
- **Manuscript Action:** Refutes any claim of soil-specific phylogenetic conservation. §3.1 must present SMAG purely as an incremental technical reference, striking all speculation regarding soil-specific evolutionary adaptations.

### Go / No-Go Boundaries:
- **GO:**
  - Storage allocation on IBEX project scratch confirms $\ge 80$ GB free space.
  - Complete Zenodo 8223844 download finishes with 100% MD5 verification across all 373 parts.
  - GEM soil MAG metadata join unambiguously identifies $\ge 15,000$ high-quality soil MAGs.
- **NO-GO:**
  - Storage quota exhausted on project scratch before completing Zenodo download.
  - Confounding between biome and phylum within GEM (e.g. if GEM soil lacks major desert Actinobacteriota, rendering the biome comparison a disguised taxonomic mismatch). In that event, fall back to taxonomic stratification within the existing release-matched framework.
