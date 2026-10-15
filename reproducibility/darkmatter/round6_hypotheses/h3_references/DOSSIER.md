# Evidence Dossier (v3): Hypothesis 3 — Bounded Reference-Release Comparison Under Exact Common-Support Length and Redundancy Matching

**Version:** 3.0 (Post-Review 1 Repair — Fully Executed & Verified on IBEX)  
**Date:** 2026-10-13  
**Assigned Local Directory:** `reproducibility/darkmatter/round6_hypotheses/h3_references`  
**Assigned Remote Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references`  
**Executed Slurm Jobs:** `53700569` (Top-1 vs. Reduced-DB Search Benchmark), `53704317` (Superseded v2 Pilot), `53707059` (v3 Repair Execution, 0 errors)  

---

## 1. Bounded Hypothesis, Estimand, and Scope Reconciliation

### 1.1 Background and Narrowed Scope from `round5/EXTENSIONS.md`
In §3.1 of the manuscript (`06_darkmatter/specificity_clean/results/`), Empty Quarter protein family representatives ($N = 106,042$ dark families in `FAMILIES_FINAL.txt`; $N = 106,042$ annotated control families in `CONTROL_FAMILIES.txt`) searched with `diamond blastp --sensitive --evalue 1e-3 -k 1` against full reference catalogues show:
- **SMAG (Ma et al. 2023, Zenodo 8223844; 40,039 soil MAGs, 129,550,340 proteins):** Dark hit rate **59.53%** (63,122 / 106,042; 40.47% absence); control hit rate **98.43%** (104,382 / 106,042; 1.57% absence).
- **GEM (Nayfach et al. 2021; 52,515 multi-biome MAGs, 122,978,191 proteins):** Dark hit rate **35.62%** (37,767 / 106,042; 64.38% absence); control hit rate **97.03%** (102,892 / 106,042; 2.97% absence).

Because `/ibex/user/velazqal/eqmags_31/smag/smag_db.dmnd` has restricted permissions (`chmod 700`) and nucleotide assemblies for the 52,515 GEM MAGs are not staged on IBEX (only published Prodigal `.faa.gz` files exist in `/ibex/scratch/projects/c2014/rob/metagenomics/databases/gem/faa/faa/`), `v3` tests an explicit **bounded hypothesis**:
> **Bounded Hypothesis (H3 v3):** Conditional on the MD5-verified public Zenodo Record 8223844 SMAG archive-order genome pool (212 complete MAGs from `mag.tar.gz.aa`–`ab`, recalled with Pyrodigal v3.7.1 `-p meta --no-stop-codon`) and the staged 52,515 published GEM Prodigal proteomes, after canonicalizing sequences (stripping trailing `*` and uppercasing), deduplicating exact canonical sequences within each reference pool, applying an explicit cross-release overlap rule, and enforcing **exact integer-length ($1\text{ aa}$) common-support quota matching** so that both arms have identical unique protein counts, identical total amino acid residues, and identical length distributions, does the SMAG reference subset yield a different paired hit rate than matched GEM reference subsets for Empty Quarter dark protein families under identical `diamond blastp` settings?

### 1.2 Explicit Scope Restrictions
1. **Archive-Order Convenience Conditioning:** The 212 complete SMAG MAGs extracted from `mag.tar.gz.aa` and `mag.tar.gz.ab` form a fixed archive-order convenience pool, not a random sample of all 40,039 SMAG genomes. Findings are conditioned on this verified pool and cannot be generalized to full catalogues without downloading all 373 Zenodo archive parts.
2. **Caller- and Composition-Conditional Estimand:** SMAG proteins are called with Pyrodigal v3.7.1 (`-p meta --no-stop-codon`), whereas GEM proteins come from the published GEM Prodigal v2.6.3 release. This design evaluates reference-release subset composition after exact sequence-count, length-distribution, and duplicate matching, conditional on published caller releases and subset taxonomic composition. It does **not** establish causal soil-biome adaptation, geographic endemism, or caller-independent biological effects.
3. **Strict Separation of Dark and Control Cohorts:** Annotated control families differ fundamentally in length distribution ($62.72\%$ $<100\text{ aa}$ vs $47.75\%$ in dark families) and evolutionary conservation. Controls are never used to calibrate or subtract technical bias from dark families, and `dark + ctrl` queries are never pooled.

---

## 2. Executed Repairs for Review 1 Critical Findings (Slurm Job `53707059`)

### 2.1 Safe Streaming Tar Extraction and Deterministic SMAG Provenance (`scripts/08_v3_extract_and_call_smag.py`)
- **Verified Primary Archive Parts (`results/v3_smag_zenodo_parts_manifest.tsv`):**
  - `mag.tar.gz.aa`: `104,857,600` bytes, MD5 `5fb8a52b97940d1d1c14aecda57ac5bc` (matches `zenodo_8223844.json`), SHA-256 `a627c9acda5bdac7d95d7f6f4e6501a5a795d912f256d0a3f98ffea468406f39`.
  - `mag.tar.gz.ab`: `104,857,600` bytes, MD5 `de5f1ecb60b8ceac42b616e2780b2d04` (matches `zenodo_8223844.json`), SHA-256 `0c395de38fd07ecbf25ce5503ba0d8090edc0c8b622a870c17eb93ed9c43c5e8`.
- **Streaming `tarfile` Byte-Size Verification (`results/v3_smag_tar_extraction_audit.tsv`):**
  - Streamed the ordered parts via Python's `tarfile.open(fileobj=..., mode="r|gz")` with path-traversal checks and exact header byte verification (`observed_bytes == member.size`).
  - **Members 1–212 (`MAG/SRR5579859rbin.24.fa` [member 1] through `MAG/SRR5437667rbin.15.fa` [member 212] in tar stream order in `results/v3_smag_tar_extraction_audit.tsv`; sorted alphabetically starting at `MAG/D3rbin.13.fa` [member 92] in `results/v3_smag_genome_manifest.tsv`):** All 212 members satisfied `observed_bytes == declared_bytes`, valid FASTA syntax, and valid IUPAC nucleotide alphabets (`status = COMPLETE_RETAINED`). All 212 `.fa` files are permanently retained in `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/genomes/smag_extracted/` (zero files deleted) with individual SHA-256 checksums recorded in `results/v3_smag_genome_manifest.tsv`.
  - **Member 213 (`MAG/SRR5855438rbin.1.fa`):** Explicitly caught and rejected (`declared_bytes = 12,783,848`, `observed_bytes = 6,946,816`, partial SHA-256 `ce0de439fa2a20809c53c895940a1a593145ca3dda4cd5dcc79c99b6fe289663`, `starts_with_gt = True`, `ends_with_newline = False`, `status = REJECTED_TRUNCATED_STREAM_END`, `detail = ReadError: unexpected end of data`).
- **Deterministic Pyrodigal v3.7.1 Gene Calling:**
  - Executed `pyrodigal v3.7.1 -p meta --no-stop-codon` across the 212 complete MAGs and sorted all output records deterministically by `(mag_id, gene_idx)` before writing:
    - `results/v3_smag_pool_canonical.faa`: **693,227** canonical proteins (`202,894,612` aa; `0` trailing `*`, `0` internal `*`, `0` lowercase/invalid chars; `231,281,034` bytes; SHA-256 `f3e6572942e83ab393167a0576ef6a252ccd23db7f9a8f2ba3b8335af75ac5ee`).
    - `results/v3_smag_pool_unique.faa`: **692,194** unique canonical proteins (`202,751,091` aa; `1,033` exact duplicates [`0.15%`] mapped in `results/v3_smag_duplicate_map.tsv`; `231,094,200` bytes; SHA-256 `4863eb745ab3065d0fa53077fe7cd735660ef80203088ef8ded5b813821e1023`).

### 2.2 Stop-Codon Canonicalization, Deduplication, Cross-Release Overlap, and Exact Length-Quota Matching (`scripts/09_v3_match_gem_and_search.py`)
- **GEM Candidate Pools Across 3 Deterministic Seeds (`results/v3_gem_genome_manifest_seed{20261013,20261014,20261015}.tsv`):**
  - For each seed $s \in \{20261013, 20261014, 20261015\}$, whole GEM `.faa.gz` files were sampled deterministically until the deduplicated canonical GEM pool reached $\ge 692,194$ unique proteins (matching the deduplicated SMAG pool size). Every sampled `.faa.gz` file is recorded with its byte size, SHA-256 digest, and stop-codon audit.
  - Across the three seeds, **94.70%–95.69% of published GEM Prodigal proteins ended with a literal `*` stop codon** (`678,136 / 709,463` in seed `20261013`; `685,307 / 716,158` in seed `20261014`; `673,565 / 711,230` in seed `20261015`), with `0` internal `*` characters.
  - After canonicalization (`seq.strip().rstrip("*").upper()`) and exact full-sequence SHA-256 deduplication, exact cross-release shared sequences (`results/v3_cross_release_shared_sequences.tsv`) are **250** (`0.036%`, seed `20261013`), **34** (`0.005%`, seed `20261014`), and **85** (`0.012%`, seed `20261015`)—correcting the false `0` (`0.00%`) overlap in v2 caused by unstripped `*` characters.
- **Exact Integer-Length ($1\text{ aa}$) Common-Support Quota Matching (`results/v3_matching_closure_audit.tsv`):**
  - For every integer amino acid length $L$, setting $q_s(L) = \min(N_{\text{SMAG, uniq}}(L), N_{\text{GEM, uniq}, s}(L))$ and sampling $q_s(L)$ unique canonical proteins per arm achieves **100% programmatic closure** (verified by fail-closed assertions in `09_v3_match_gem_and_search.py:329-338`):

| Run / Seed | Overlap Policy | GEM MAGs | GEM Pool Raw / Uniq (Redund. %) | Pre / Post Shared Seqs | Matched Unique Proteins (SMAG == GEM) | Matched Residues (SMAG == GEM) | Mean / SD / Median Length (aa) | Matched 5-Stratum Counts (`<100`, `100–149`, `150–249`, `250–399`, `>=400` aa) |
|---|---|---|---|---|---|---|---|---|
| `seed20261013` | Retained | 299 | 709,463 / 693,721 (2.22%) | 250 / 217 | **657,779 == 657,779** | **197,005,061 == 197,005,061** | 299.50 / 204.21 / 260.00 | `80,327`, `76,523`, `154,335`, `190,364`, `156,230` (both arms) |
| `seed20261013_no_overlap` | Excluded | 299 | 709,463 / 693,721 (2.22%) | 250 / **0** | **657,529 == 657,529** | **196,961,439 == 196,961,439** | 299.55 / 204.22 / 260.00 | `80,253`, `76,459`, `154,272`, `190,332`, `156,213` (both arms) |
| `seed20261014` | Retained | 308 | 716,158 / 696,211 (2.79%) | 34 / 33 | **660,049 == 660,049** | **197,201,755 == 197,201,755** | 298.77 / 203.92 / 259.00 | `81,304`, `77,118`, `155,348`, `190,076`, `156,203` (both arms) |
| `seed20261015` | Retained | 312 | 711,230 / 693,028 (2.56%) | 85 / 78 | **660,421 == 660,421** | **196,954,586 == 196,954,586** | 298.23 / 202.92 / 259.00 | `81,315`, `77,853`, `155,509`, `189,820`, `155,924` (both arms) |

---

## 3. Rough Heuristic Threshold Diagnostic for Post-Hoc Subsampling vs. Fresh Subset Search (`results/v3_top1_vs_dbsize_decomposition.tsv`)

In Job `53700569` (`scripts/03_pilot_subsample_search.py`), searching the 2,000-query panel against a 2,500-genome GEM subset (`5,815,177` proteins, `1,834,739,150` database letters) yielded 126 dark hits and 953 control hits, whereas post-hoc filtering of full-GEM (`122,978,191` proteins) top-1 rows (`gem_hits.m8`) recovered only 16 dark hits (110 missed; 87.30% discrepancy) and 58 control hits (895 missed; 93.91% discrepancy). Reducing the searched database alters both `-k 1` secondary-hit masking and E-value scaling. As an ancillary **rough heuristic threshold diagnostic** (not an exact attribution to `-k 1` masking vs. database size), `results/v3_top1_vs_dbsize_decomposition.tsv` partitions the empirical subset hits at the nominal protein-count-rescaled threshold $E_{\text{sub}} \le 10^{-3} / (122,978,191 / 5,815,177) = 4.73 \times 10^{-5}$:
- **Dark Discrepancy ($n = 110$ queries):** **79 / 110 (71.82%)** satisfy $E_{\text{sub}} \le 4.73 \times 10^{-5}$ (all 79 had a different top-1 subject in full GEM), while **31 / 110 (28.18%)** fall in $4.73 \times 10^{-5} < E_{\text{sub}} \le 10^{-3}$ (including 5 queries that had no hit in full GEM).
- **Control Discrepancy ($n = 895$ queries):** **885 / 895 (98.88%)** satisfy $E_{\text{sub}} \le 4.73 \times 10^{-5}$, while **10 / 895 (1.12%)** fall in $4.73 \times 10^{-5} < E_{\text{sub}} \le 10^{-3}$ (`0` queries lacked a full-GEM hit).
- **Methodological Caveat:** This protein-count E-value split cannot identify the causal mechanism for individual queries because Karlin-Altschul / DIAMOND E-values depend on effective searched amino-acid residue space ($n' = N_{\text{letters}} - N_{\text{seqs}} \cdot l_{\text{eff}}(m)$ with query-length edge correction), composition-based adjustments/masking, and block/seed search heuristics rather than raw sequence counts. The existing TSV column names in `results/v3_top1_vs_dbsize_decomposition.tsv` (`strict_top1_masked_e_le_4p73e5`, `borderline_dbsize_rescaled_e_gt_4p73e5`, `strict_top1_masked_share_pct`, `borderline_dbsize_rescaled_share_pct`) are retained historical diagnostic labels, not validated per-query mechanism classifications. This diagnostic is ancillary to the primary v3 experiment (§4), where matched SMAG and GEM databases have identical protein counts, identical total amino-acid residues, and identical exact integer-length distributions.

---

## 4. Paired Search Results Under Exact Length and Redundancy Matching

### 4.1 Separate Cohort Summaries: Stratified Population-Weighted vs. Finite-Panel (`results/v3_matched_cohort_summary.tsv`)

Using exact population stratum weights $W_h = N_h / N$ and inclusion probabilities $\pi_h = n_h / N_h$ ($n_h = 200$ per stratum; $N = 106,042$ per cohort), the primary Horvitz-Thompson / stratified population-weighted marginal estimates and unweighted finite-panel ($n = 1,000$) estimates are:

| Run Label | Hit Rule | Overlap Policy | Cohort | Pop-Weighted SMAG Hit % (SE) | Pop-Weighted GEM Hit % (SE) | Pop-Weighted $\hat{\Delta}_{\text{wt}}$ % (FPC SE) | Pop-Weighted 95% CI (%) | Pop-Weighted Wald $p$-val | Pop-Weighted Disc. Ratio | Finite-Panel ($n=1000$) Both / SMAG-only ($b$) / GEM-only ($c$) / Neither | Panel Exact McNemar $p$-val | Panel $b/c$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `seed20261013` | `E<=1e-3` | Retained | **Dark** | **15.63%** (1.46%) | **5.33%** (0.99%) | **+10.30%** (1.47%) | **[+7.42, +13.18]** | $2.38 \times 10^{-12}$ | **6.34×** | 25 / **127** / **13** / 835 (15.20% vs 3.80%) | $1.14 \times 10^{-24}$ | **9.77×** |
| `seed20261014` | `E<=1e-3` | Retained | **Dark** | **15.52%** (1.46%) | **4.23%** (0.88%) | **+11.29%** (1.46%) | **[+8.43, +14.14]** | $9.21 \times 10^{-15}$ | **8.68×** | 22 / **128** / **13** / 837 (15.00% vs 3.50%) | $6.30 \times 10^{-25}$ | **9.85×** |
| `seed20261015` | `E<=1e-3` | Retained | **Dark** | **15.75%** (1.47%) | **5.69%** (0.99%) | **+10.06%** (1.41%) | **[+7.30, +12.81]** | $8.15 \times 10^{-13}$ | **6.46×** | 26 / **124** / **19** / 831 (15.00% vs 4.50%) | $4.44 \times 10^{-20}$ | **6.53×** |
| `seed20261013_no_overlap` | `E<=1e-3` | Excluded | **Dark** | **15.44%** (1.45%) | **5.33%** (0.99%) | **+10.11%** (1.46%) | **[+7.24, +12.98]** | $5.15 \times 10^{-12}$ | **6.24×** | 25 / **125** / **13** / 837 (15.00% vs 3.80%) | $3.77 \times 10^{-24}$ | **9.62×** |
| `seed20261013_id30` | `E<=1e-3 & id>=30%` | Retained | **Dark** | **14.72%** (1.44%) | **4.99%** (0.98%) | **+9.74%** (1.46%) | **[+6.88, +12.59]** | $2.28 \times 10^{-11}$ | **6.32×** | 21 / **111** / **12** / 856 (13.20% vs 3.30%) | $3.03 \times 10^{-21}$ | **9.25×** |
| `seed20261013` | `E<=1e-3` | Retained | **Ctrl** | **85.41%** (1.81%) | **78.50%** (2.06%) | **+6.90%** (1.76%) | **[+3.46, +10.34]** | $8.44 \times 10^{-5}$ | **3.56×** | 876 / 55 / 14 / 55 (93.10% vs 89.00%) | $6.92 \times 10^{-7}$ | **3.93×** |
| `seed20261014` | `E<=1e-3` | Retained | **Ctrl** | **85.58%** (1.80%) | **79.08%** (2.06%) | **+6.50%** (1.95%) | **[+2.67, +10.33]** | $8.77 \times 10^{-4}$ | **2.63×** | 879 / 51 / 22 / 48 (93.00% vs 90.10%) | $9.14 \times 10^{-4}$ | **2.32×** |
| `seed20261015` | `E<=1e-3` | Retained | **Ctrl** | **85.89%** (1.78%) | **80.20%** (2.01%) | **+5.69%** (1.85%) | **[+2.07, +9.32]** | $2.08 \times 10^{-3}$ | **2.46×** | 880 / 51 / 24 / 45 (93.10% vs 90.40%) | $2.44 \times 10^{-3}$ | **2.12×** |
| `seed20261013_no_overlap` | `E<=1e-3` | Excluded | **Ctrl** | **85.41%** (1.81%) | **78.43%** (2.06%) | **+6.97%** (1.76%) | **[+3.53, +10.42]** | $7.20 \times 10^{-5}$ | **3.59×** | 875 / 56 / 14 / 55 (93.10% vs 88.90%) | $4.30 \times 10^{-7}$ | **4.00×** |
| `seed20261013_id30` | `E<=1e-3 & id>=30%` | Retained | **Ctrl** | **85.06%** (1.82%) | **78.32%** (2.06%) | **+6.74%** (1.76%) | **[+3.30, +10.19]** | $1.26 \times 10^{-4}$ | **3.35×** | 869 / 56 / 17 / 58 (92.50% vs 88.60%) | $5.27 \times 10^{-6}$ | **3.29×** |

- **Monte Carlo Reference-Draw Sensitivity (3 Primary Seeds `20261013`–`20261015`):**
  - **Dark Cohort:** Population-weighted $\hat{\Delta}_{\text{wt}}$ averages **+10.55 percentage points** (range **+10.06% to +11.29%**; population-weighted discordant ratio range **6.34× to 8.68×**; finite-panel discordant ratio range **6.53× to 9.85×**).
  - **Control Cohort:** Population-weighted $\hat{\Delta}_{\text{wt}}$ averages **+6.36 percentage points** (range **+5.69% to +6.90%**), driven almost entirely by the `<100` aa stratum ($W_1 = 0.6272$), whereas strata $\ge 150\text{ aa}$ show near-identical recovery ($94.5\%$ vs $93.0\text{--}94.5\%$ in $150\text{--}249\text{ aa}$; $100.0\%$ vs $98.5\text{--}100.0\%$ in $\ge 250\text{ aa}$).

### 4.2 Per-Length-Stratum Breakdown for Primary Seed `20261013` (`results/v3_matched_length_stratified_results.tsv`)

| Cohort | Length Stratum (aa) | $n_h$ | $N_h$ | Inclusion Prob $\pi_h$ | Pop Weight $W_h$ | Both | SMAG Only ($b_h$) | GEM Only ($c_h$) | Neither | SMAG Hit % | GEM Hit % | $\hat{\delta}_h$ % (FPC SE) | Exact McNemar $p$-val | Discordant Ratio ($b_h/c_h$) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Dark** | `<100` | 200 | 50,634 | 0.003950 | 0.477490 | 10 | **24** | 7 | 159 | 17.00% | 8.50% | **+8.50%** (2.72%) | $3.33 \times 10^{-3}$ | **3.43×** |
| **Dark** | `100–149` | 200 | 19,720 | 0.010142 | 0.185964 | 2 | **23** | 0 | 175 | 12.50% | 1.00% | **+11.50%** (2.25%) | $2.38 \times 10^{-7}$ | `undefined` |
| **Dark** | `150–249` | 200 | 21,301 | 0.009389 | 0.200873 | 6 | **25** | 1 | 168 | 15.50% | 3.50% | **+12.00%** (2.40%) | $8.05 \times 10^{-7}$ | **25.00×** |
| **Dark** | `250–399` | 200 | 9,782 | 0.020446 | 0.092246 | 3 | **27** | 2 | 168 | 15.00% | 2.50% | **+12.50%** (2.52%) | $1.62 \times 10^{-6}$ | **13.50×** |
| **Dark** | `>=400` | 200 | 4,605 | 0.043431 | 0.043426 | 4 | **28** | 3 | 165 | 16.00% | 3.50% | **+12.50%** (2.59%) | $4.65 \times 10^{-6}$ | **9.33×** |
| **Ctrl** | `<100` | 200 | 66,508 | 0.003007 | 0.627185 | 136 | 24 | 7 | 33 | 80.00% | 71.50% | **+8.50%** (2.72%) | $3.33 \times 10^{-3}$ | **3.43×** |
| **Ctrl** | `100–149` | 200 | 14,789 | 0.013524 | 0.139464 | 162 | 20 | 2 | 16 | 91.00% | 82.00% | **+9.00%** (2.25%) | $1.21 \times 10^{-4}$ | **10.00×** |
| **Ctrl** | `150–249` | 200 | 15,346 | 0.013033 | 0.144716 | 181 | 8 | 5 | 6 | 94.50% | 93.00% | **+1.50%** (1.79%) | $5.81 \times 10^{-1}$ | **1.60×** |
| **Ctrl** | `250–399` | 200 | 6,907 | 0.028956 | 0.065135 | 197 | 3 | 0 | 0 | 100.00% | 98.50% | **+1.50%** (0.85%) | $2.50 \times 10^{-1}$ | `undefined` |
| **Ctrl** | `>=400` | 200 | 2,492 | 0.080257 | 0.023500 | 200 | 0 | 0 | 0 | 100.00% | 100.00% | **0.00%** (0.00%) | $1.00 \times 10^{0}$ | `undefined` |

*(Full 50-row stratum table across all 5 runs is retained in `results/v3_matched_length_stratified_results.tsv`.)*

---

## 5. Scientific Interpretation and Optional Catalogue-Scale Extensions

1. **Exact Length and Redundancy Matching Closes the Technical-Size Confound Within the Bounded Pool:** Even when both reference arms are canonicalized, deduplicated to `0.00%` internal redundancy, and matched to the exact integer amino-acid length distribution (`657,779` to `660,421` unique proteins and `196.95M` to `197.20M` amino acid residues per arm), the 212-MAG SMAG subset recovers a population-weighted **15.52%–15.75%** of dark families compared with **4.23%–5.69%** in the matched GEM subsets ($\hat{\Delta}_{\text{wt}} = +10.06\%$ to $+11.29\%$, all $p \le 8.15 \times 10^{-13}$).
2. **Persistence Across Full-Length Protein Strata and Hit-Rule / Overlap Sensitivities:** Unlike annotated controls—where discordant recovery is confined to $<150\text{ aa}$ queries and vanishes for long proteins ($\ge 150\text{ aa}$)—dark-family excess recovery in the SMAG subset is +10.00% to +13.00% across all medium and long strata ($150\text{--}249$, $250\text{--}399$, and $\ge 400\text{ aa}$) across all three seeds, when shared exact sequences are excluded (`seed20261013_no_overlap`, $\hat{\Delta}_{\text{wt}} = +10.11\%$), and when $\ge 30\%$ identity is enforced during DIAMOND search (`seed20261013_id30`, $\hat{\Delta}_{\text{wt}} = +9.74\%$).
3. **Strict Inferential Boundary:** These results establish a bounded, caller- and pool-conditional reference-composition difference between the verified 212-MAG SMAG archive-order pool (called with Pyrodigal v3.7.1) and length/redundancy-matched GEM Prodigal subsets. They do **not** prove a causal soil-biome adaptation, geographic endemism, or caller-independent effect, nor do they generalize to unexamined parts of the 373-part SMAG archive.
4. **Optional Catalogue-Scale Additions (Outside Bounded Viability Requirements):**
   - Downloading all 373 Zenodo 8223844 parts (~39 GB) to sample randomly across all 40,039 SMAG MAGs.
   - Fetching external nucleotide assemblies for GEM MAGs to re-call GEM with Pyrodigal v3.7.1.
   - Stratifying GEM by `ecosystem == "Soil"` vs non-soil biomes with taxonomic and 95% ANI / MMseqs2 dereplication controls.
