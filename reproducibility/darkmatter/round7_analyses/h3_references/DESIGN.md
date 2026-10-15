# Experimental Design: Hypothesis 3 (Round 7 Full Analysis)

## Bounded Reference-Release Comparison Under Broad Random Sampling, True Caller Harmonization, and Exact Common-Support Length Matching

**Date:** 2026-10-14  
**Assigned Local Directory:** `reproducibility/darkmatter/round7_analyses/h3_references/`  
**Assigned Remote Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/`  
**Authors / Operators:** Antigravity (Gemini Flash 3.8 High) for Robert Hoehndorf  
**Target Manuscript:** Nature Communications (Dark-matter manuscript, Section 3.1)  

---

## 1. Scientific Context and Purpose

In §3.1 of the manuscript (`06_darkmatter/specificity_clean/`), Empty Quarter protein family representatives were searched against two major public reference catalogues:
1. **SMAG (Soil Metagenome-Assembled Genomes; Ma et al. 2023, Nature Communications 14:7318; Zenodo Record 8223844):** 40,039 MAGs, 129,550,340 proteins; dark hit rate **59.53%** (40.47% absence), control hit rate **98.43%** (1.57% absence).
2. **GEM (Genomes from Earth's Microbiomes; Nayfach et al. 2021, Nature Biotechnology 39:499–509; NERSC portal):** 52,515 multi-biome MAGs, 122,978,191 proteins; dark hit rate **35.62%** (64.38% absence), control hit rate **97.03%** (2.97% absence).

### Key Weaknesses in Round 6 Pilot:
1. **Convenience Prefix Conditioning:** Round 6 evaluated only the first 212 complete SMAG MAGs extracted from the first two Zenodo archive parts (`mag.tar.gz.aa` and `mag.tar.gz.ab`), which represented an archive-order convenience pool rather than a random sample of the full 40,039-genome catalogue.
2. **Un-Harmonized Gene Callers:** SMAG was called with Pyrodigal v3.7.1 (`-p meta --no-stop-codon`), whereas GEM was conditioned on published Prodigal v2.6.3 proteomes because GEM nucleotide assemblies had not yet been located.
3. **Finite Query Sample:** Round 6 searched a stratified panel of 2,000 queries (1,000 dark + 1,000 control) rather than the complete 212,084-query catalogue.

### Scope of Round 7 Full Analysis:
1. **Broad Random Sampling Across Full Releases:**
   - SMAG: Draw a broad, reproducibly random sample of $M = 1,500$ genomes uniformly at random without replacement from the complete 40,039 catalogue MAGs (sampling frame established by `Supplementary Data 2.xlsx`, seed `20261014`).
   - GEM: Draw matching candidate pools from the complete 52,515 GEM catalogue across multiple random seeds (`20261014`, `20261015`, `20261016`).
2. **True Caller Harmonization:**
   - Locate and download the nucleotide contigs for the sampled GEM cohort directly from the public NERSC repository (`https://portal.nersc.gov/GEM/genomes/fna/`).
   - Call genes on BOTH sides using the EXACT same Pyrodigal v3.7.1 binary and parameters (`-p meta --no-stop-codon`).
   - Also evaluate the published GEM Prodigal v2.6.3 proteomes as an end-to-end gene-calling and matching sensitivity. Changing callers also changes the eligible protein pool and length quotas, so this contrast does not isolate a pure caller effect.
3. **Exact Common-Support 1-aa Length Quota Matching:**
   - Canonicalize all sequences (`strip().rstrip('*').upper()`), audit trailing/internal `*`, uppercase, deduplicate within-pool (`0.00%` internal redundancy).
   - Enforce exact integer-length ($1\text{ aa}$) common-support quotas:
     $$q_s(L) = \min\big(N_{\text{SMAG, uniq}}(L),\; N_{\text{GEM, uniq}, s}(L)\big)$$
   - Both databases guaranteed identical protein counts, identical total residues, identical length distributions.
4. **Full Catalogue Query Census:**
   - Search ALL 212,084 queries (106,042 dark + 106,042 annotated controls) in `query_all.faa`.
   - Eliminates finite query sampling variance ($\text{SE}_{\text{query}} = 0$).
   - Report separate dark and control arms.

---

## 2. Sampling Frame, Pre-Declared Seeds, and Estimands

### 2.1 Pre-Declared Sampling Specifications
- **SMAG Sampling Frame:** Exactly 40,039 unique MAG IDs from `Supplementary Data 2.xlsx` (saved in `results/smag_catalogue_40039_ids.txt`).
- **SMAG Cohort Size:** $M = 1,500$ MAGs sampled uniformly at random with pre-declared seed `20261014`. Exact replay: `sorted(random.Random(20261014).sample(sorted(frame), 1500))`. Independent replay and direct comparison with every `user_genome` in the original strict-OOXML Supplementary Data 2 pass in `round7/root_checks/audit_h3_sampling.py` and `h3_sampling_audit.json`.
- **GEM Sampling Frame:** Exactly 52,515 unique MAG IDs from `databases/gem/faa/faa/` (saved in `results/gem_catalogue_52515_ids.txt`).
- **GEM Cohort Size:** Exactly 2,055 genomes per seed, selected before matching using the fixed-count rule `ceil(4,500,000 * 1.05 / 2,300)`. This is not an adaptive stopping rule at 4.5 million observed proteins.
- **GEM Cohort Seeds:**
  - Primary Seed: `20261014` (both Pyrodigal v3.7.1 harmonized calling and published Prodigal v2.6.3 proteomes).
  - Sensitivity Seeds: `20261015` and `20261016`, each with both harmonized Pyrodigal and published Prodigal proteomes. Harmonized sensitivity draws were added before any full SMAG comparison result was available; exact sampled genome IDs are shared across callers. The fixed SMAG draw is reused, so these runs assess GEM-draw and matching sensitivity, not SMAG-draw uncertainty.

### 2.2 Comparison Estimands
For each reference pair (SMAG vs GEM), searched against all $N = 106,042$ dark queries and $N = 106,042$ control queries:
1. **Marginal Hit Rates:**
   - $\hat{P}_{\text{SMAG}} = \frac{a + b}{N}$
   - $\hat{P}_{\text{GEM}} = \frac{a + c}{N}$
2. **Paired Difference:**
   - $\hat{\Delta} = \hat{P}_{\text{SMAG}} - \hat{P}_{\text{GEM}} = \frac{b - c}{N}$
   - Where $a$ = Both hit, $b$ = SMAG only, $c$ = GEM only, $d$ = Neither hit.
   - Optional exact two-sided binomial McNemar test on discordant pairs $(b, c)$ describes a model-based paired-query null assuming independent discordances; related families can violate this assumption. Exact finite-cohort counts and differences are primary, not nominal significance.
   - Discordant ratio: $\psi = b / c$ (formatted as `undefined` if $c = 0$).
3. **Per-Length-Stratum Breakdown:**
   - Strata: `<100`, `100–149`, `150–249`, `250–399`, `\ge 400` aa.
   - Denominator in each stratum is the exact catalogue stratum count $N_h$.
4. **Sensitivity Estimands:**
   - **Caller Harmonization:** SMAG (Pyrodigal v3.7.1) vs GEM (Pyrodigal v3.7.1) vs GEM (Prodigal v2.6.3).
   - **Cross-Release Overlap:** Overlap-retained vs overlap-excluded matched databases.
   - **Identity Threshold:** DIAMOND search with `--id 30` vs default E-value threshold ($E \le 10^{-3}$).

---

## 3. Inferential Boundaries and Non-Causal Framing

1. **Reference-Draw Sensitivity, Not Biological Replication:**
   - Random reference draws from GEM (seeds `20261014`, `20261015`, `20261016`) evaluate sensitivity to reference-subset selection and matching noise; they are NOT independent biological replicates of the Earth's biomes.
2. **Release-Specific Observational Differences:**
   - Any difference in hit rate between SMAG and GEM reflects the specific genomic contents, biome representations, and assembly characteristics of the two published catalogues.
   - We NEVER claim a causal soil-biome adaptation, geographic endemism, or catalogue-wide selection regime from reference hit-rate differences.
3. **Query Census (SE = 0):**
   - Searching all 212,084 prespecified recurrent-family and annotated-control representatives eliminates finite query subsampling uncertainty for that enumerated cohort. It is not a census of all 43 million protein representatives. Reference sampling, assembly/gene-calling errors, heuristic search and dependence/generalization uncertainties remain.

---

## 4. Slurm Execution and Resource Budget

- **Environment & Binaries:**
  - Python: `/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python`
  - Pyrodigal v3.7.1: `/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal`
  - DIAMOND v2.0.15: `/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond`
- **Slurm Configuration:**
  - Partition: `batch`, Account: `pi-hohndor`
  - Download jobs: recorded Slurm workers use urllib HTTPS streaming, strict public checksums and automatic retry/backoff. Later auxiliary workers resume validated HTTP ranges; exact jobs and concurrency changes are recorded in round7 STATE.json and ROOT_AUDIT_NOTES.md.
  - Gene calling: 32 CPU cores, time limit 12h.
  - Matched search: 32 CPU cores, time limit 12h.
