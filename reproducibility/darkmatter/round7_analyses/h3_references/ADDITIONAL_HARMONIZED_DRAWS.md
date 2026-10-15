# H3 Reference Catalogue Specificity: Additional Harmonized Comparator Draws

**Date:** 2026-10-14  
**Author:** Gemini Flash 3.8  
**Hypothesis:** H3 (Reference Catalogue Specificity & Comparator Harmonization)  
**Assigned Directory:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/`  
**Active Job ID:** `53742584` (`h3_gem_draws`)

---

## 1. Executive Summary & Scientific Design Frame

### 1.1 Fixed-Count Sampling Rule (2,055 Genomes per Draw)
In accordance with the pre-declared protocol in `DESIGN.md` and clarified in `round7/ROOT_AUDIT_NOTES.md`:
- The candidate GEM comparator draws are selected from the 52,515-genome GEM catalogue (`gem_catalogue_52515_ids.txt`) using pseudo-random shuffles seeded by integer dates (`20261014`, `20261015`, `20261016`).
- The sampling rule selects a **fixed sample size of 2,055 genomes** per draw, derived as:
  $$\lceil \frac{4,500,000 \times 1.05}{2,300} \rceil = 2,055 \text{ genomes}$$
- It is **not** an adaptive stopping rule evaluated dynamically during streaming; the sample frame is fixed, deterministic, and identical between the published Prodigal evaluation and the harmonized Pyrodigal evaluation.

### 1.2 Nature of the Estimands: Reference-Draw Sensitivity vs. Biological Replicates
- Seeds `20261014`, `20261015`, and `20261016` evaluate **Monte Carlo reference-draw sensitivity** across different random draws of 2,055 genomes from the global GEM catalogue.
- They are **not biological replicates** and cannot be treated as independent experimental units for inferential pooling; they assess the stability of the SMAG vs. GEM contrast under stochastic variation in reference representation.

### 1.3 Methodological Scope: End-to-End Caller-Plus-Matching Sensitivity
- Comparisons between the harmonized Pyrodigal cohort (SMAG Pyrodigal vs. GEM Pyrodigal) and the published Prodigal cohort (SMAG Pyrodigal vs. GEM Prodigal) assess **end-to-end caller-plus-matching sensitivity**.
- Because exact 1-aa integer length quota matching (`04_exact_length_matching.py`) conditions on the length distribution of predicted proteins—which differs subtly between Pyrodigal v3.7.1 and Prodigal v2.6.3—the matched protein sequence content varies between caller specifications.
- This contrast reflects the realistic, end-to-end sensitivity of the pipeline to reference caller choice and length quota matching, rather than an isolated, causal gene-caller effect holding reference sequence content strictly constant.

---

## 2. Implementation Architecture (`scripts/06_harmonize_additional_draws.py`)

The additional harmonized draws workflow operates under strict non-destructive constraints:

### 2.1 Immutability of Existing Primary Artifacts
- **Primary Seed 20261014**: Primary retained genecalls in `genecalls/gem/` (`.raw.faa`, `.gff`, `.id_map.tsv`) and the extended manifest `gem_harmonized_manifest_seed20261014.with_coordinates.tsv` are treated as immutable baseline artifacts.
- **Shared Genomes**:
  - Seed 20261015 shares 65 genomes with seed 20261014.
  - Seed 20261016 shares 76 genomes with seed 20261014.
  - Seed 20261015 and seed 20261016 share 75 genomes with each other.
- **Reuse without Overwrite**: For any genome already called and verified in seed 20261014, the existing raw FAA, GFF, and ID-to-coordinate mapping files are reused directly. Pyrodigal is **never** re-run on these genomes, avoiding race conditions and preventing GFF header discrepancies (where temporary input path strings embedded by Pyrodigal can alter GFF SHA256 hashes even if predicted sequences are identical).

### 2.2 Assembly Download from NERSC Portal
- For the 3,894 genomes not in seed 20261014, nucleotide contigs (`.fna.gz`) are downloaded concurrently from the public NERSC repository (`https://portal.nersc.gov/GEM/genomes/fna/{gid}.fna.gz`).
- Downloads use atomic temporary files (`.tmp_<pid>_<ts>`) with automatic exponential-backoff retries (up to 5 attempts) and gzip readability validation.

### 2.3 Harmonized Pyrodigal Calling
- **Caller Binary**: `/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal` (v3.7.1).
- **Biologically Identical Flags**: `-p meta --no-stop-codon -f gff`.
- **Outputs Retained in `genecalls/gem/`**:
  - `{gid}.raw.faa`: raw protein translations with original Prodigal headers.
  - `{gid}.gff`: gene coordinates with exact start, end, strand, and contig attributes.
  - `{gid}.id_map.tsv`: mapping table linking canonical protein IDs (`{gid}::{gid}_{gene_idx}`) to original contigs, start/end coordinates, and strand.

### 2.4 Artifact Generation per Seed
For both `seed20261015` and `seed20261016`:
1. `gem_harmonized_pool_canonical_seed{seed}.faa`: full canonical protein pool (sorted by genome ID).
2. `gem_harmonized_pool_unique_seed{seed}.faa`: deduplicated pool of unique canonical proteins.
3. `gem_harmonized_duplicate_map_seed{seed}.tsv`: mapping of duplicate sequences to canonical representatives.
4. `gem_harmonized_manifest_seed{seed}.tsv`: legacy 11-column manifest matching seed 20261014 format.
5. `gem_harmonized_manifest_seed{seed}.with_coordinates.tsv`: extended manifest recording `raw_faa_path`, `raw_faa_sha256`, `gff_path`, `gff_sha256`, `id_map_path`, and `id_map_sha256`.
6. `gem_harmonized_protein_coordinate_mapping_seed{seed}.tsv`: aggregated coordinate mapping table across all 2,055 genomes in the draw.
7. `gem_harmonized_draws_summary.json`: JSON summary recording exact counts, residue totals, redundancy rates, and SHA256 hashes of all generated artifacts.

---

## 3. Downstream Candidate Updates

Candidate scripts have been created under `reviews/darkmatter-2026-10-13/round7/h3_additional_draws/` for root review and promotion:

1. **`04_exact_length_matching.py`**:
   - Incorporates `harmonized_seed20261015` and `harmonized_seed20261016` pairs alongside `harmonized_seed20261014`, `harmonized_seed20261014_no_overlap`, `published_seed20261014`, `published_seed20261015`, and `published_seed20261016`.
   - Produces matched FASTA pairs (`matched_smag_*.faa` and `matched_gem_*.faa`).
2. **`05_diamond_search_full_census.py`**:
   - Incorporates paired DIAMOND searches for `harmonized_seed20261015` and `harmonized_seed20261016`.
   - Strictly preserves all signature-pinned caching (`.dmnd.sig.json`, `.m8.sig.json`), fail-closed semantics on corrupt or unverified files, valid zero-hit result support, and the exact 212,084 full query census validation.

---

## 4. Execution Tracking

- **Slurm Job ID**: `53742584`
- **Script**: `scripts/06_harmonize_additional_draws.py`
- **Batch Wrapper**: `scripts/06_harmonize_additional_draws.sbatch`
- **Output Log**: `logs/gem_draws_53742584.out`
- **Error Log**: `logs/gem_draws_53742584.err`
- **Status JSON**: `STATUS_ADDITIONAL_DRAWS.json`
- **Exit Status**: Completed successfully (`0:0`), total runtime: 1,592.0s (~26.5 minutes).
- **Search Dependency Recommendation**: When promoting candidate matching/search scripts, update job `53738997` dependency to include `53742584` (`afterok:53738996:53738968:53742584`).

---

## 5. Completed Empirical Results & Provenance Checksums

### 5.1 Harmonized Comparator Summary Table

| Metric | Primary Seed 20261014 | Harmonized Seed 20261015 | Harmonized Seed 20261016 |
| :--- | :--- | :--- | :--- |
| **Genomes (fixed count)** | 2,055 | 2,055 | 2,055 |
| **Selection Rule** | `ceil(4.5M*1.05/2300)` | `ceil(4.5M*1.05/2300)` | `ceil(4.5M*1.05/2300)` |
| **Overlap with Seed 20261014** | — (baseline) | 65 genomes (3.16%) | 76 genomes (3.70%) |
| **Canonical Proteins** | 4,824,709 | 4,693,590 | 4,845,330 |
| **Canonical Residues (aa)** | 1,523,280,088 | 1,482,250,212 | 1,529,947,999 |
| **Unique Canonical Proteins** | 4,306,926 | 4,220,573 | 4,337,716 |
| **Unique Residues (aa)** | 1,376,624,176 | 1,347,953,920 | 1,386,748,080 |
| **Internal Duplicates** | 517,783 | 473,017 | 507,614 |
| **Sequence Redundancy Rate** | 10.73% | 10.08% | 10.48% |
| **Gene Caller** | Pyrodigal v3.7.1 | Pyrodigal v3.7.1 | Pyrodigal v3.7.1 |
| **Caller Flags** | `-p meta --no-stop-codon` | `-p meta --no-stop-codon` | `-p meta --no-stop-codon` |

### 5.2 Artifact Provenance & SHA256 Checksums

#### Seed 20261015 Artifacts
- **Canonical Pool FASTA**: `results/gem_harmonized_pool_canonical_seed20261015.faa`  
  `SHA256: 52d1d55594ee179fd58faeeb3ef6ef96fa7ead1d85bb458139478ea5bbac0215`
- **Unique Pool FASTA**: `results/gem_harmonized_pool_unique_seed20261015.faa`  
  `SHA256: f40af369dd4c6327ed030cac26bd06751f88e7d3db51618ecdc27397d1b75db3`
- **Duplicate Map**: `results/gem_harmonized_duplicate_map_seed20261015.tsv`  
  `SHA256: 716c4d14a2257dd2cedbf890b3adcd77ef1a3e522c7f27c9e6e7b37d25d338c2`
- **Legacy Manifest**: `results/gem_harmonized_manifest_seed20261015.tsv`  
  `SHA256: a92ae52252dd536f08dc10f5b0a808db847e326b3b9c31e5695bf3e7fffeaf0a`
- **Extended Manifest (with coordinates)**: `results/gem_harmonized_manifest_seed20261015.with_coordinates.tsv`  
  `SHA256: 925f0ece06e34563dd0013f9635d145691d9fa97cd0abfe0cb23df6654510fc8`
- **Coordinate Mapping**: `results/gem_harmonized_protein_coordinate_mapping_seed20261015.tsv`  
  `SHA256: fa2a546cc47d1373aec2fcf1b03d01d07a0b7eb766036906a470c8003e3fc2b0`

#### Seed 20261016 Artifacts
- **Canonical Pool FASTA**: `results/gem_harmonized_pool_canonical_seed20261016.faa`  
  `SHA256: 18dc6026c28f1b5a8ac284fa1477a8eafa4cd5834fe433c416553aadd7437a9c`
- **Unique Pool FASTA**: `results/gem_harmonized_pool_unique_seed20261016.faa`  
  `SHA256: d227531a46b1e6ca47d20a3b436afed7e9e06ca0e9eff007215cb19f687b473f`
- **Duplicate Map**: `results/gem_harmonized_duplicate_map_seed20261016.tsv`  
  `SHA256: 8116f78d18005e9d5d5e26ca28c6dc647559d5c06c61f76c47612249fcd2a5a4`
- **Legacy Manifest**: `results/gem_harmonized_manifest_seed20261016.tsv`  
  `SHA256: a44a5c8e8967cb35ea38679a90998b127c6f051aa801e04d08d5e05760020bd6`
- **Extended Manifest (with coordinates)**: `results/gem_harmonized_manifest_seed20261016.with_coordinates.tsv`  
  `SHA256: 5e46c10186ecb913fb0c7b39c478e0aa10ecbae26522c24e5fa6bc335a8d50c4`
- **Coordinate Mapping**: `results/gem_harmonized_protein_coordinate_mapping_seed20261016.tsv`  
  `SHA256: 0cee87bf5cb94dabe00e75491f746e5580dd96d34fac6480beccd0d57195bb21`

