# Experimental Design (v2): Hypothesis 3 — Soil-Reference Composition After Technical Matching

**Version:** 2.0 (Post-Audit Revision)  
**Date:** 2026-10-13  
**Author:** Gemini Flash High (for Robert Hoehndorf)  
**Assigned Local Directory:** `reproducibility/darkmatter/round6_hypotheses/h3_references`  
**Assigned Remote Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references`  

---

## 1. Scientific Context and Scope Reconciliation

### 1.1 Background and Problem Statement
In §3.1 of the manuscript, Empty Quarter dark protein families are queried against reference genome catalogues. The two comparably sized genome catalogues are:
- **SMAG (Ma et al. 2023, *Nat Commun* 14:7318):** 129,550,340 proteins from 40,039 global soil MAGs; 40.47% absence (59.53% hit rate).
- **GEM (Nayfach et al. 2021, *Nat Biotechnol* 39:499–509):** 122,978,191 proteins from 52,515 multi-biome MAGs; 64.38% absence (35.62% hit rate).

Reviewer Argon erroneously inverted absence and hits, inventing a 4.93% GEM hit rate and an "8.2-fold" contrast. In reality, SMAG recovers 25,355 more dark families (+23.91 percentage points in hit rate) than GEM at near-identical catalogue scale (129.6M vs 123.0M sequences, a 1.053× ratio).

### 1.2 Dossier Corrections and Scope Narrowing (per Root Audit)
Following root audit feedback, the experimental scope and inferential claims are strictly bounded:
1. **No Calibrated Technical Bias from Controls:** Control discordance (3.57×) does *not* prove that technical bias explains a proportional share of dark discordance (12.19×). Control families differ drastically from dark families in length distribution (62.7% of controls are $<100$ aa vs 47.7% of dark families) and in conservation regime. We explicitly reject attributing the difference to "residual lineage signal" or asserting short-query effects are pure caller artifacts.
2. **Defensible Estimand (Conditional on Caller):** Because GEM nucleotide assemblies are not available in project scratch (only translated `.faa` files exist in `databases/gem`), GEM cannot be re-called with Pyrodigal without external genome re-fetching. Therefore, the honest, falsifiable estimand is:
   $$\text{"Reference-release composition under matched sequence count and length distribution, conditional on published caller releases"}$$
   This estimand is explicitly narrower than a causal biome effect or a fully harmonized-caller comparison.
3. **Conditioning on Archive-Order Convenience Subset:** For the bounded pilot, we download the first ordered chunks of the public Zenodo SMAG archive (Zenodo record 8223844), extract only 100% complete tar genome members, and openly condition on this archive-order convenience subset. This establishes source access, gene-calling, and search execution feasibility without claiming to represent a random draw from global soils.

---

## 2. Statistical Estimand, Units, and Hypotheses

- **Primary Unit of Analysis:** The orthologous protein family representative ($n = 106,042$ dark families; paired with $n = 106,042$ length/prevalence-matched annotated control families).
- **Deterministic Pilot Query Set:** A stratified random sample of $n = 2,000$ representatives ($1,000$ dark, $1,000$ control; $200$ per length stratum: $<100$, $100–149$, $150–249$, $250–399$, $\geq 400$ aa; seed `20261013`).
- **Primary Estimand:**
  $$\Delta_{\text{matched}} = P(\text{Hit in SMAG}_{\text{sub}} \mid \text{query } i) - P(\text{Hit in GEM}_{\text{sub}} \mid \text{query } i)$$
  measured under identical alignment settings (`diamond blastp --sensitive --evalue 1e-3 -k 1`) on physically matched database slices.
- **Matched Criteria:**
  - Sequence count matched within $\pm 0.5\%$.
  - Total amino acid residue count and length distributions recorded and compared.
  - Query set fixed and paired.
- **Null Hypothesis ($H_0$):**  
  Under technical matching of database sequence count and length distributions, the hit rate for Empty Quarter dark families in the SMAG reference subset does not exceed that in the GEM reference subset ($\Delta_{\text{matched}} \leq 0$).
- **Alternative Hypothesis ($H_1$):**  
  The SMAG reference subset maintains a higher hit rate for dark families than the GEM reference subset ($\Delta_{\text{matched}} > 0$).

---

## 3. Methodological Design and Provenance

### 3.1 Public SMAG Source Retrieval and Verification
- Primary source: Zenodo record 8223844 ([Ma et al. 2023](https://doi.org/10.5281/zenodo.8223844)), CC-BY-4.0.
- Archive parts: `mag.tar.gz.aa`, `mag.tar.gz.ab` (each 104,857,600 bytes).
- Checksums verified against public Zenodo API metadata (`reviews/darkmatter-2026-10-13/round6/public_source_checks/zenodo_8223844.json`):
  - `mag.tar.gz.aa`: MD5 `5fb8a52b97940d1d1c14aecda57ac5bc`
  - `mag.tar.gz.ab`: MD5 `de5f1ecb60b8ceac42b616e2780b2d04`
- Download executed directly on IBEX with checksum verification and retry handling.

### 3.2 Integrity-Protected Genome Extraction
- The stream of concatenated parts is decompressed and parsed using Python's `tarfile` module.
- Strict completeness guard: Tar members are verified against their header byte size (`member.size`). Any member cut off by the end of the archive boundary is discarded with explicit logging. Zero truncated genomes are admitted into the analysis.
- Retained genomes are written to individual `.fa` files and recorded in `results/smag_pilot_genome_manifest.tsv`.

### 3.3 Gene Calling on Fresh SMAG Genomes
- Tool: `Pyrodigal v3.7.1` (path `/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal` / Python module).
- Parameters: `-p meta --no-stop-codon` (identical to the project's Empty Quarter pipeline and original SMAG compilation).
- Header convention: `<MAG_ID>::<gene_index>` matching existing catalogue conventions.
- Compiled into `results/smag_pilot_recalled.faa`.

### 3.4 GEM Matched-Subset Sampling
- Source: 52,515 individual genome `.faa.gz` files in `/ibex/scratch/projects/c2014/rob/metagenomics/databases/gem/faa/faa/`.
- Sampling: Deterministic genome sampling (seed `20261013`) to achieve exact sequence count parity ($\pm 0.5\%$) with the extracted SMAG proteins.
- Compiled into `results/gem_pilot_matched.faa`.

### 3.5 Paired DIAMOND Search and Analysis
- Both FASTA files are indexed via `diamond makedb` (v2.0.15).
- Both databases are queried with the identical deterministic 2,000-query sample (`pilot_query_2000.faa`) using identical flags:
  `diamond blastp -q pilot_query_2000.faa -d {db} --sensitive -e 1e-3 -k 1 -p {threads}`.
- Outputs scored and decomposed into:
  - 4-quadrant contingency tables (Both, SMAG-sub only, GEM-sub only, Neither).
  - Length-stratified hit rates.
  - McNemar paired test and symmetry ratio.

---

## 4. Resource Allocation and Compute Policy

- All compute executed on IBEX under project account `--account=pi-hohndor`.
- Pilot capped at: $\leq 8$ CPUs, $\leq 32$ GB RAM, $\leq 30$ min walltime on partition `debug` or `batch`.
- Temporary staging files confined strictly to `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/tmp/`.
- No modifications to upstream databases, shared scripts, or other hypotheses.
