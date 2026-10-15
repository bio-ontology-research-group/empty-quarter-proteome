# Experimental Design (v3): Hypothesis 3 — Bounded Reference-Release Comparison Under Exact Common-Support Length and Redundancy Matching

**Version:** 3.0 (Post-Review 1 Repair — Written Prior to v3 Execution)  
**Date:** 2026-10-13  
**Assigned Local Directory:** `reproducibility/darkmatter/round6_hypotheses/h3_references`  
**Assigned Remote Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references`  
**Superseded Versions Preserved:** `DESIGN_v2_superseded.md`, `DOSSIER_v2_superseded.md`, `INPUTS_v2_superseded.tsv`, `scripts/01`–`07`, and all v2 `results/` and `logs/`  

---

## 1. Scientific Scope and Bounded Estimand

### 1.1 Historical Context and Narrowed Scope from `round5/EXTENSIONS.md`
In §3.1 of the manuscript (`06_darkmatter/specificity_clean/results/`), Empty Quarter protein family representatives ($N = 106,042$ dark families in `FAMILIES_FINAL.txt`; $N = 106,042$ length- and prevalence-matched annotated control families in `CONTROL_FAMILIES.txt`) were searched with `diamond blastp --sensitive --evalue 1e-3 -k 1` against full reference catalogues:
- **SMAG (Ma et al. 2023, Zenodo 8223844):** 129,550,340 proteins across 40,039 soil MAGs; dark hit rate **59.53%** (40.47% absence), control hit rate **98.43%** (1.57% absence).
- **GEM (Nayfach et al. 2021):** 122,978,191 proteins across 52,515 multi-biome MAGs; dark hit rate **35.62%** (64.38% absence), control hit rate **97.03%** (2.97% absence).

Because local access to `/ibex/user/velazqal/eqmags_31/smag/smag_db.dmnd` is restricted (`chmod 700`) and nucleotide assemblies for the 52,515 GEM MAGs are not staged on IBEX (only published Prodigal `.faa.gz` files are staged in `/ibex/scratch/projects/c2014/rob/metagenomics/databases/gem/faa/faa/`), v3 defines an explicit, falsifiable **bounded estimand**:
> **Bounded Hypothesis (H3 v3):** Conditional on the MD5-verified public Zenodo Record 8223844 SMAG archive-order genome pool (212 complete MAGs from `mag.tar.gz.aa`–`ab`, recalled with Pyrodigal v3.7.1 `-p meta --no-stop-codon`) and the staged 52,515 published GEM Prodigal proteomes, after canonicalizing sequences (stripping trailing `*` and uppercasing), deduplicating exact canonical sequences within each reference pool, applying an explicit cross-release overlap rule, and enforcing **exact integer-length ($1\text{ aa}$) common-support quota matching** so that both arms have identical unique protein counts, identical total amino acid residues, and identical length distributions, does the SMAG reference subset yield a different paired hit rate than matched GEM reference subsets for Empty Quarter dark protein families under identical `diamond blastp` settings?

### 1.2 Explicit Inferential Boundaries
1. **Convenience Conditioning, Not Whole-Catalogue Random Sampling:** The 212 SMAG MAGs extracted from the first two ordered archive parts (`mag.tar.gz.aa`, `mag.tar.gz.ab`) constitute a fixed archive-order convenience pool, not a simple random sample of all 40,039 SMAG genomes.
2. **Caller- and Taxonomy-Conditional Comparison:** SMAG proteins are called with Pyrodigal v3.7.1 (`-p meta --no-stop-codon`), whereas GEM proteins come from the published GEM Prodigal v2.6.3 release. The contrast measures reference-release subset composition conditional on published/pipeline callers and subset taxonomic composition; it does **not** isolate a causal soil-biome adaptation, geographic endemism, or a caller-independent biological effect.
3. **No Control Calibration of Technical vs. Biological Shares:** Annotated control families differ fundamentally from dark families in conservation and length distribution ($62.72\%$ of controls vs $47.75\%$ of dark families are $<100\text{ aa}$). Control results are reported strictly as a separate positive-sensitivity reference cohort and are **never** used to calibrate or subtract technical bias from dark families, nor pooled with dark families.
4. **Rough Heuristic Diagnostic for Post-Hoc Subsampling vs. Fresh Subset Search Discrepancy:** Comparing historical full-database (`123.0M` proteins) top-1 rows filtered to a 2,500-genome GEM subset against a fresh DIAMOND search of that 2,500-genome (`5.82M` proteins, `1.83B` letters) subset (`scripts/03_pilot_subsample_search.py`, Job `53700569`) reflects **both** `-k 1` secondary-hit masking and database-size E-value rescaling. Partitioning the discrepancy at the nominal protein-count-rescaled threshold $E_{\text{sub}} \le 10^{-3} / 21.15 = 4.73 \times 10^{-5}$ (`results/v3_top1_vs_dbsize_decomposition.tsv`, where $71.8\%$ of dark and $98.9\%$ of control discrepancies have $E_{\text{sub}} \le 4.73 \times 10^{-5}$ and $28.2\%$ / $1.1\%$ fall in $4.73 \times 10^{-5} < E_{\text{sub}} \le 10^{-3}$) is a **rough heuristic threshold diagnostic**, not an exact per-query attribution to `-k 1` masking vs. database size, because effective residue search space, composition adjustments, and DIAMOND search heuristics differ (and the TSV column names `strict_top1_masked_*` / `borderline_dbsize_rescaled_*` are historical diagnostic labels rather than validated mechanism classifications). Physical subset searches are therefore mandatory, and in the primary matched experiment (§2.2) both arms are matched to identical protein counts, total residues, and integer-length distributions.

---

## 2. Pre-Declared Protocol and Repairs for Review 1 Critical Findings

### 2.1 Safe Streaming Tar Extraction and Deterministic SMAG Gene Calling (`scripts/08_v3_extract_and_call_smag.py`)
1. **Input Verification:** Verify byte sizes, MD5 checksums (against `reviews/darkmatter-2026-10-13/round6/public_source_checks/zenodo_8223844.json`), and SHA-256 digests of `mag.tar.gz.aa` (`104,857,600` bytes) and `mag.tar.gz.ab` (`104,857,600` bytes). Record in `results/v3_smag_zenodo_parts_manifest.tsv`.
2. **Streaming Tar Integrity Gate:** Read the ordered concatenated byte stream directly via Python's `tarfile.open(fileobj=..., mode="r|gz")` (no shell `tar -xz ... || true`).
   - Reject any non-regular file or unsafe path (`startswith("/")`, `".."` in path).
   - For every member $m$, record 1-based `member_order_idx`, `tar_member_name`, and `declared_bytes = m.size`.
   - Read member payload bytes and compare `observed_bytes = len(data)` against `declared_bytes`.
   - Any member with `observed_bytes != declared_bytes`, a `tarfile.ReadError`/`EOFError`, or invalid FASTA syntax is **rejected** and logged in `results/v3_smag_tar_extraction_audit.tsv` (retaining complete members 1–212 from `MAG/SRR5579859rbin.24.fa` through `MAG/SRR5437667rbin.15.fa` and specifically rejecting the truncated 213th member `MAG/SRR5855438rbin.1.fa`, `declared_bytes = 12,783,848`, `observed_bytes = 6,946,816`, `status = REJECTED_TRUNCATED_STREAM_END`).
   - Retain all 212 complete extracted genome `.fa` files permanently in `genomes/smag_extracted/` on IBEX (no deletion), recording their `nucleotide_sha256`, `contig_count`, and `total_bp`.
3. **Deterministic Pyrodigal Gene Calling and Canonicalization:**
   - Call genes on each retained `.fa` using `/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal` (v3.7.1) with `-p meta --no-stop-codon`.
   - Collect per-genome results and **sort deterministically by `(mag_id, gene_idx)`** before writing `results/v3_smag_pool_canonical.faa` and `results/v3_smag_pool_unique.faa`.
   - Audit every protein for trailing `*`, internal `*`, lowercase characters, and non-standard characters. Canonicalize via `seq.strip().rstrip("*").upper()`.
   - Record complete per-genome provenance (including `raw_faa_sha256` and `canonical_faa_sha256`) in `results/v3_smag_genome_manifest.tsv`.

### 2.2 Canonicalization, Deduplication, Cross-Release Overlap, and Exact Length-Quota Matching Across 3 Seeds (`scripts/09_v3_match_gem_and_search.py`)
1. **Sequence Canonicalization & Stop-Codon Audit:**
   - All GEM `.faa.gz` sequences read from `/ibex/scratch/projects/c2014/rob/metagenomics/databases/gem/faa/faa/` are audited for trailing `*` (present on ~95.6% of published GEM Prodigal proteins), internal `*`, and case before canonicalizing via `seq.strip().rstrip("*").upper()`. Any sequence containing an internal `*` is logged and excluded.
2. **Repeated Deterministic Candidate-Pool Draws (3 Seeds: `20261013`, `20261014`, `20261015`):**
   - For each seed $s \in \{20261013, 20261014, 20261015\}$, shuffle the sorted list of 52,515 GEM `.faa.gz` filenames using `random.Random(s)` and draw whole GEM genomes sequentially until the **deduplicated canonical GEM sequence count** reaches at least the deduplicated canonical SMAG sequence count ($\ge 692,194$ unique canonical proteins).
   - Emit a full file-level manifest `results/v3_gem_genome_manifest_seed{s}.tsv` enumerating every sampled GEM `.faa.gz` file, its byte size, `.faa.gz` SHA-256 digest, raw protein count, trailing-`*` count, internal-`*` count, canonical residue count, and new unique canonical sequences contributed.
   - *Interpretation of Repeated Draws:* Each seed $s$ draws a fresh deterministic sample of GEM genomes against the fixed 212-MAG SMAG archive-order pool. These 3 draws quantify reference-subset selection and matching sensitivity (Monte Carlo variation across GEM genome draws), **not** independent biological replication of whole catalogues.
3. **Within-Pool Deduplication and Cross-Release Overlap Policy:**
   - Deduplicate both the SMAG pool and each seed's GEM candidate pool by exact full-sequence SHA-256 / MD5 digest of the canonical amino acid string (`seq_canon`), retaining the first deterministic sequence ID (`<mag_id>::<mag_id>_<gene_idx>` for SMAG; `<gem_id>` in sorted filename + record order for GEM).
   - Audit exact canonical sequence overlap $O_s = S_{\text{uniq}} \cap G_{\text{uniq}, s}$ and write all shared sequences and their SMAG/GEM IDs to `results/v3_cross_release_shared_sequences.tsv`.
   - **Primary Overlap Policy:** Retain each release's deduplicated pool prior to common-support length matching (shared sequences reflect genuine cross-release sequence sharing and can only contribute concordant `Both Hit` outcomes).
   - **Sensitivity Overlap Policy:** Also execute an **overlap-excluded** construction (`S_uniq \ O_s` vs `G_uniq,s \ O_s`) for seed `20261013` to verify programmatically whether removing exact cross-release shared sequences alters discordant counts $(b, c)$.
4. **Exact Integer-Length ($1\text{ aa}$) Common-Support Quota Matching:**
   - For every exact integer amino acid length $L \in \{1, \dots, L_{\max}\}$, define the common-support quota for seed $s$:
     $$q_s(L) = \min\big(N_{\text{SMAG, uniq}}(L),\; N_{\text{GEM, uniq}, s}(L)\big)$$
   - Using `random.Random(s)`, deterministically sample exactly $q_s(L)$ unique canonical proteins of length $L$ from $S_{\text{uniq}}$ and $q_s(L)$ unique canonical proteins of length $L$ from $G_{\text{uniq}, s}$, written in deterministic sorted ID order to `v3_smag_matched_seed{s}.faa` and `v3_gem_matched_seed{s}.faa`.
   - **Programmatic Closure Guarantee:** By construction, for every seed $s$, `v3_smag_matched_seed{s}.faa` and `v3_gem_matched_seed{s}.faa` have:
     - **Identical unique protein count:** $\sum_L q_s(L)$ (`0` protein difference).
     - **Identical total amino acid residues:** $\sum_L L \cdot q_s(L)$ (`0` residue difference).
     - **Identical internal redundancy:** `0.00%` duplicates in both arms.
     - **Identical length distribution:** identical mean, SD, median, quantiles, and counts in every length stratum (`<100`, `100–149`, `150–249`, `250–399`, `>=400` aa).

### 2.3 Pre-Declared Hit Rules and Physical DIAMOND Searches
1. **Primary Hit Rule (`E <= 1e-3`):**
   - `diamond makedb` and `diamond blastp -q pilot_query_2000.faa -d {db} --sensitive --evalue 1e-3 -k 1 -p {threads}` (DIAMOND v2.0.15.153).
   - Chosen for exact parameter compatibility with upstream `specificity_clean/results/{smag,gem}_hits.m8`.
2. **Pre-Declared Secondary Sensitivity Hit Rule (`E <= 1e-3` AND `pident >= 30%` during search):**
   - Executed as a separate physical DIAMOND search passing `--id 30` (`diamond blastp -q pilot_query_2000.faa -d {db} --sensitive --evalue 1e-3 --id 30 -k 1 -p {threads}`) on the seed `20261013` matched databases so DIAMOND retrieves the best alignment meeting $\ge 30\%$ sequence identity directly during search rather than post-hoc filtering `-k 1` rows.

### 2.4 Cohort Separation, Stratum Analysis, and Design-Based Population Weighting
1. **Strict Cohort Separation:** Dark queries ($n = 1,000$ representing $N = 106,042$ dark families) and annotated control queries ($n = 1,000$ representing $N = 106,042$ control families) are analyzed strictly separately. **No pooled `dark + ctrl` row is computed.**
2. **Per-Stratum Estimands ($h \in \{1,\dots,5\}$; $n_h = 200$ per stratum):**
   - Report $a_h$ (`both`), $b_h$ (`smag_only`), $c_h$ (`gem_only`), $d_h$ (`neither`), stratum hit rates $\hat{p}_{A,h} = (a_h+b_h)/n_h$ and $\hat{p}_{B,h} = (a_h+c_h)/n_h$, paired difference $\hat{\delta}_h = (b_h - c_h)/n_h$ with finite-population-corrected standard error $\text{SE}(\hat{\delta}_h) = \sqrt{(1 - n_h/N_h) s^2_{d,h} / n_h}$, exact two-sided binomial McNemar $p$-value (`binomtest(b_h, b_h + c_h, 0.5)`), and discordant ratio $b_h / c_h$ (recorded as `undefined` whenever $c_h = 0$, never `inf`).
3. **Stratified Population-Weighted Marginal Estimands (Primary Cohort Summary):**
   - Using exact population stratum counts $N_h$ from `FAMILIES_FINAL.txt` and `CONTROL_FAMILIES.txt` ($\pi_h = n_h / N_h$, $W_h = N_h / N$):
     - Dark $N_h$: `<100`: $50,634$ ($W_1 = 0.477490$), `100–149`: $19,720$ ($W_2 = 0.185964$), `150–249`: $21,301$ ($W_3 = 0.200873$), `250–399`: $9,782$ ($W_4 = 0.092246$), `>=400`: $4,605$ ($W_5 = 0.043426$).
     - Control $N_h$: `<100`: $66,508$ ($W_1 = 0.627185$), `100–149`: $14,789$ ($W_2 = 0.139464$), `150–249`: $15,346$ ($W_3 = 0.144716$), `250–399`: $6,907$ ($W_4 = 0.065135$), `>=400`: $2,492$ ($W_5 = 0.023500$).
   - Compute Horvitz-Thompson / stratified population-weighted hit rates $\hat{P}_A^{\text{wt}} = \sum_h W_h \hat{p}_{A,h}$, $\hat{P}_B^{\text{wt}} = \sum_h W_h \hat{p}_{B,h}$, paired difference $\hat{\Delta}_{\text{wt}} = \sum_h W_h \hat{\delta}_h$, design-based standard error $\text{SE}(\hat{\Delta}_{\text{wt}}) = \sqrt{\sum_h W_h^2 (1 - n_h/N_h) s^2_{d,h}/n_h}$, 95% CI, and population-weighted discordant ratio $\hat{\psi}_{\text{wt}} = [\sum_h W_h (b_h/n_h)] / [\sum_h W_h (c_h/n_h)]$.
   - Unweighted finite-panel totals ($n=1,000$ equal-quota panel) are reported separately with explicit labeling.
4. **Informative Outcomes in Both Directions:**
   - **Positive ($\hat{\Delta}_{\text{wt}} > 0$):** Demonstrates that within this bounded, caller-conditional comparison of the 212-MAG SMAG archive-order pool vs matched GEM proteomes, excess dark-family recovery persists after exact common-support matching of unique protein count, residue volume, and integer-length distribution.
   - **Null or Reversed ($\hat{\Delta}_{\text{wt}} \le 0$):** Demonstrates that exact deduplication and length-stratum matching eliminate or reverse the apparent recovery excess in this subset. Both outcomes are scientifically informative and viable without requiring statistical significance.
