# RQ-F2 — Historical exploratory MIBiG-relative BGC clustering

**Status.** This documents an archived workflow and its operational no-MIBiG-member statistic. Because 90.1% of
the EQ antiSMASH regions are contig-edge fragments, the result is a fragmentation-sensitive upper bound and is
not promoted to the active manuscript. Here “no-reference GCF” means only that a cluster contains no member of
the named MIBiG release under the stated feature and threshold rule; it is not a biological novelty call.

## Input data

- **EQ BGCs:** `$M/06_metabolism/antismash/<MAG>/*.region*.gbk` — **119,090** antiSMASH-8.0.4
  region GBKs across the **15,469** single-sample HQ+MQ soil MAGs (one antiSMASH dir per MAG).
  antiSMASH was run with `--minimal --genefinding-tool prodigal` (job `86_antismash_array.sh`);
  the region GBKs carry prodigal CDS + `/translation` + `gene_kind`, which is all BiG-SLICE
  needs (it re-scans biosynthetic Pfam/sub-Pfam HMMs on the CDS itself).
- **Reference:** **MIBiG 3.1** experimentally-characterized BGCs (`mibig_gbk_3.1`, **2,502**
  GBKs; `dl.secondarymetabolites.org/mibig/mibig_gbk_3.1.tar.gz`). MIBiG is the reference set
  that contributes characterized references to the wider biosynthetic landscape. Co-clustering against MIBiG is
  not equivalent to querying the much larger BiG-FAM predicted-BGC collection.

## Tool choice — BiG-SLICE (not BiG-SCAPE)

119,090 BGCs is far beyond a tractable BiG-SCAPE all-vs-all (≫20k). **BiG-SLICE v2.0.2**
(kautsar2020bigslice) was designed for >1.2 M BGCs: it converts each BGC to a biosynthetic
sub-Pfam feature vector and clusters with BIRCH (super-linear), then assigns membership to GCF
centroids. We **co-cluster EQ + MIBiG in one Mode-1 run** so MIBiG BGCs land in the same GCF
space; a GCF with ≥1 EQ BGC and **0** MIBiG members is recorded operationally as a no-MIBiG-member GCF.

- Tool: `bigslice` 2.0.2 (pip), conda env `bigslice` under
  `/ibex/scratch/projects/c2014/rob/miniconda3/envs/bigslice` (python 3.9, hmmer 3.3.2,
  numpy 1.23, scikit-learn 1.1, **pyhmmer 0.10.15**).
- HMM models: `download_bigslice_hmmdb` → `bigslice-models.2022-11-30` (biosynthetic_pfams +
  sub_pfams), staged at `…/envs/bigslice/bin/bigslice-models` and passed via
  `--program_db_folder`.
- Clustering threshold **T = 0.4** (BiG-SLICE v2 default for normalized features; the v2 GCF
  cutoff). BIRCH, deterministic.
- `-t 32`, `--mem=300G`, `--partition=batch --account=pi-hohndor`, walltime 2 days.

### Two compatibility fixes (documented for reproducibility)

1. **pyhmmer pin.** bigslice 2.0.2 leaves `pyhmmer` unpinned; pip pulled 0.12.1 whose
   `hit.name` is `str`, crashing bigslice's `….name.decode()`. Pinned **pyhmmer<0.11**
   (→ 0.10.15), which restores the bytes API. Validated: full run finishes.
2. **MIBiG version-string patch.** MIBiG 3.1 GBKs carry `##antiSMASH-Data## Version :: False`.
   bigslice 2.0.2 only parses GBKs whose antiSMASH major version ∈ {5,6,7,8}, so it silently
   skipped *all* MIBiG files. Patched the version line to `6.1.1`
   (`sed -i 's/Version      :: False/Version      :: 6.1.1/' mibig_gbk_3.1/*.gbk`); the MIBiG
   `subregion aStool="mibig"` branch then parses correctly (2,502 BGCs, CDS + `/label` class).
   This only changes a metadata string, not sequences or features.
- A scipy/HiGHS `GLIBCXX` issue on the login node is avoided in the job by
  `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`.

## Pipeline

1. `bigslice_input/datasets.tsv` (tab-sep, empty taxonomy column):
   `eq → $M/06_metabolism/antismash`, `mibig → …/mibig_gbk_3.1`. BiG-SLICE globs
   `<folder>/**/*.gbk` and keeps only names matching `^BGC[0-9]{7}$`, `*.cluster[0-9]+`, or
   `*.region[0-9]+` — so the per-MAG full-genome `<MAG>.gbk` files are auto-ignored; genome
   identity = parent dir (= MAG id).
2. `bigslice -i bigslice_input bigslice_out --threshold 0.4 -t 32 --program_db_folder <models>`
   → biosyn_pfam scan → sub_pfam scan → feature extraction → BIRCH GCF models → membership
   assignment → `bigslice_out/result/data.db`.
3. `87_bgc_novelty_analyze.py data.db results/` → parses the SQLite result DB. Reference =
   any dataset whose name contains `mibig`; everything else = EQ. Assigns each BGC to its
   best-rank GCF, counts EQ/MIBiG members per GCF, computes novel-GCF fractions overall and
   per biosynthetic class (`bgc_class → chem_subclass → chem_class`).
4. `make_fig_bgc_novelty.py` (local) plots % novel GCFs overall + by class from the pulled
   tables.

## Validation

End-to-end smoke test on 22 EQ region GBKs (1 MAG) + 12 MIBiG GBKs via `srun` (batch
partition): full pipeline completed (RUN_FINISHED) and `87_bgc_novelty_analyze.py` produced
correct per-class tables (Terpene/RiPP/NRP/Polyketide/Saccharide/Other), confirming both the
clustering and the novelty parsing before the full submission.

## Scripts (copies on IBEX at `$M/scripts/ibex/`)

- `87_bgc_gcf_novelty.sbatch` — submission + clustering + analysis driver.
- `87_bgc_novelty_analyze.py` — result-DB → novelty tables/summary.
- `make_fig_bgc_novelty.py` — figure (run locally after pulling `results/`).

## Caveats

- antiSMASH `--minimal` detects BGC regions but skips KnownClusterBlast; novelty here is from
  GCF membership vs MIBiG, not from antiSMASH's own known-cluster hits.
- **107,316/119,090 (90.1%)** of EQ regions sit on contig edges. Fragmentation can split regions and distort
  feature vectors, systematically inflating the no-MIBiG-member fraction. A complete-region sensitivity run is
  required before any conservative novelty estimate.
- MIBiG (~2.5k characterized BGCs) is a *characterized*-biosynthesis reference; a full BiG-FAM
  query (≈1.2 M predicted BGCs) would estimate no-hit rates relative to a much broader named database
  and is the planned extension (BiG-SLICE Mode-2 `--query` against the BiG-FAM v1 model with
  `--no_normalize_feature`).

**References:** kautsar2020bigslice (BiG-SLICE), kautsar2021bigfam (BiG-FAM),
blin2023antismash (antiSMASH), terlouw2023mibig (MIBiG 3.0). Keys in `manuscript/references.bib`
(kautsar2020bigslice + terlouw2023mibig added by this RQ).
