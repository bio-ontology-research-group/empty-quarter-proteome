# Pipeline guide — Empty Quarter MAGs

**What this is.** The forensic record of how the paper's numbers were produced: tool, version, exact
parameters, input, output, job ID. It is the transcription source for Methods M.1–M.6.

> 🛑 **DO NOT REWRITE THIS GUIDE.** It is finished and pushed. On 2026-08-17 three sessions ran the
> rewrite concurrently; a parallel rewrite duplicates finished work and a force-push would destroy it.
> **Add** to it, or tell Alejandra.

**Rewritten 2026-08-17** (session 7a) from the thirteen `_clean` workbench folders, verified on disk.
The June 2026 version is superseded in full — it predates every 2026-08 recompute and documented
`M/scripts/`, which is dead PI code.

> ⚠⚠ **METHODS NUMBERING CHANGED 2026-08-17 (session 7b).** This guide's Methods pointers ("M.7 is
> already written", "M.8 and M.10 are the better source", etc.) use the **OLD** numbering. The
> section-specific block was renumbered into section order: old M.7→**M.9**, M.8→**M.10**, M.9→**M.7**,
> M.10→**M.11**, M.11→**M.8**; M.12–M.14 unchanged. Mapping table and rationale in `CORRECTIONS.md`
> 2026-08-17. Read every "M.x" below through that table.

> ⚠ **`notes/PIPELINE_GUIDE.pdf` is the June 2026 build and now contradicts this file.** It was not
> regenerated with this rewrite. **Regenerate or delete it — Alejandra's call.** Until then, the `.md`
> is the only current version; do not transcribe Methods from the PDF.

---

## 0. How to read this

**Trees.**

| tree | path | holds |
|---|---|---|
| workbench `M` | `/ibex/scratch/projects/c2014/rob/metagenomics` | all compute: scripts, intermediates, results |
| paper repo `R` | `/ibex/project/c2014/alelopezv/repos/empty-quarter-mags` | manuscript, cited tables/figures, notes |
| raw data `P` | `/ibex/scratch/projects/c2014/EmptyQuarter_Data` | reads, assemblies |

**The `_clean` convention.** Every current analysis lives in one self-contained folder
`M/<area>/<theme>_clean/` with a `README.md` (STATUS block: verified numbers, job IDs, what is
settled), `scripts/` and `results/`. **A `_clean` folder is authoritative; anything outside one is
not.** Enumerate them with `find $M -maxdepth 2 -type d -name '*_clean'` — there are **thirteen**. A
`-maxdepth 1` sweep over hand-picked top-level dirs misses `04_functional/repko_clean` and
`04_protein_catalog/darkmatter_clean`, both load-bearing.

**Verification marks.**

- **V** — read off the script, log, or conda metadata on disk during this rewrite.
- **R** — README/prose claim only, not re-derived from the executed command.
- **⚠ inherited** — a consequential parameter that was **never passed**, so a tool default applied.
  These are the ones that bite in Methods.
- **UNVERIFIED** — not yet confirmed on disk. Named honestly, with the next step. Do **not**
  transcribe an UNVERIFIED row into Methods without checking it first.

**Precedence for any number:** `R/notes/CORRECTIONS.md` > `R/notes/CATALOGUE_NUMBERS.md` >
everything else. `reproducibility/claims.tsv` is the PI's machine mirror and is allowed to lag.

---

## 1. Pipeline map

```
reads ──fastp──> QC reads ──MEGAHIT──> 549 assemblies
                                          │
        ┌─────────────────────────────────┴──────────────────────┐
        │ SemiBin2 single-sample (pretrained)  SemiBin2 self-supervised │
        └─────────────────────────────────┬──────────────────────┘
                          DAS_Tool 1.1.7 (merges the two)  ──> 14,426 refined bins
                                          │  CheckM2 1.1.0 → 4,225 HQ / 8,209 MQ / 1,992 LQ
                                          │
                          dRep 3.6.2 ANI95 (12,433 in) ──> 5,438 species representatives
                                          │
        ┌─────────────────────────────────┼─────────────────────────────┐
   GTDB-Tk 2.7.2/r232              Bakta 1.12.0 + eggNOG 2.1.13    CoverM 0.7.0
   (51 phyla)                       → rep_ko_v2.tsv               (rep-level, 548)
                                          │
   assemblies ──pyrodigal 3.7.1──> 753M proteins ──MMseqs2 linclust──> 449M reps
                                          │  mmseqs vs UniRef50 (dark = no hit)
                                          ├──> 69.9M dark ──MAG-gate──> 4,292,446 genome-resolved
                                          │                                   │
                        ┌─────────────────┴──────────────┐          ┌─────────┴────────┐
              ProstT5 + FoldSeek (§1.3)         mmseqs 50/80 (§1.4) │  106,042 families │
              → 17,221 fold families            + CoverM member-level (551)
```

---

## 2. Cohorts — the denominators, disambiguated once

These near-identical numbers are the likeliest transcription error. They are **not**
interchangeable.

| cohort | N | what it is | used by |
|---|---|---|---|
| canonical soil samples | **551** | headline N, 5 expeditions | everywhere |
| assemblies | **549** | 551 minus 1 failed, 1 empty | M.2 |
| DAS_Tool input samples | **542** | **V** `05_refine/sample_list.txt` — binned by *both* runs | M.2 |
| samples producing DAS_Tool bins | **539** | **V** counted on disk (see §4) | M.2 |
| samples contributing HQ+MQ MAGs | **543** | union: single 540 ⊆ self 543 | M.2 |
| gene-prediction basis | **552** | assembly files pyrodigal ran on; Methods-level basis only | M.4, §1.2 |
| **rep-level** CoverM cohort | **548** | canonical ∩ profiled; `51Sr1`, `F5Dr1`, `F13PRr1` failed | §2.1–2.4, **M.5** |
| **member-aware** CoverM cohort | **551** | complete; separate run vs family member genes | §1.4, §2.5, M.9 |
| CoverM profiles produced | **569** | **V** dirs with a `run_provenance.json`; all fresh recomputes | §2.1–2.4 provenance |
| `profiles_full_recompute/` dirs | **573** | **V** 569 + 4 abandoned `.tmp` staging dirs — a directory count, **not** a profile count and **not** the study N | nothing; do not quote |
| amplicon-paired | **448** | canonical ∩ CoverM-profiled ∩ 16S-profiled | §2.1, M.11 |
| XRF chemistry / climate | **248** / **453** | profiled samples with soil chemistry / climate | §2.3, §2.4 |
| env matrix (§2.5) | **547** (XRF 251, climate 456) | the shipped §2.3 env table | §2.5 |
| site-paired compartment set | **44** | of 45 Surface+Deep sites, both profiled | §2.3 |
| site **labels** | **61** labels = **60** sites | **V** — 59 numbered transect sites + 2 trip1-only aliases (62, 64) of *named* locations. **The site count is 60** — see the retraction below | §3.1 |

⚠ **Two CoverM runs exist; both are current; they are different runs.** Rep-level maps reads to the
5,438 representatives (→548 usable). Member-aware maps reads to the 536,437 dark-family member genes
(→551). Never quote one cohort for the other's result.

⚠⚠ **RETRACTED 2026-08-17 (session 7b). This guide previously said “fix the manuscript to 61”. That
was wrong — **“60 sites” is correct and is the project’s settled number.** The metadata does carry 61
distinct numeric labels (IDs 1–64, gaps at 56, 61, 63), but the companion 16S ecology paper and data
descriptor both state **60** primary sampling sites on a 1,043 km transect plus **10 named**
special-interest locations, and the descriptor establishes by coordinate matching that identifiers
**61–64 appear only on Trip 1 and denote named locations, not additional physical sites**. In our 551
samples those extra labels are **62 and 64, trip1 only, 4 samples**. So: **60 sites; 61 field-notebook
labels.** M.1 states this. Do not “correct” 60 to 61 again.

---

## 3. M.1 — Sampling and sequencing  ⚠ LARGEST GAP

**V** from `M/metadata/sample_metadata.tsv` (551 rows + header):

| axis | value |
|---|---|
| samples | **551** |
| sites | **60** (61 numeric labels — see §2) |
| compartments | Surface **117** · Deep **166** · PlantRoot **268** |
| expeditions (`trip`) | trip1 **210** · trip2 **37** · trip3 **128** · trip4 **85** · trip5 **91** |
| columns available | `sample, site, compartment, replicate, trip, trip_num, assembly_dir, contigs, R1, R2, SE, has_contigs, has_reads` |

**UNVERIFIED — not present in the workbench and required by M.1:** DNA extraction kit, library
preparation, sequencing platform and read layout, sequencing depth, sequencing centre, sampling
depths in cm per compartment, expedition dates, ENA/SRA accessions, XRF instrument and method.

**Next step:** these belong to the companion RAK-PR data paper, not this workbench. Get them from
Alejandra or the RAK-PR repo — do **not** invent them, and do not infer a platform from read length.

---

## 4. M.2 — Assembly, binning and refinement

**Upstream source of truth = the Snakemake workflow at `/ibex/project/c2014/alelopezv/EQ_WGS`**
(`assembly.smk`, `binning.smk`, `rules/{assembly,binning}/*.rule`, `config/config.yml`, per-rule
conda pins in `envs/*.yml`). It supersedes the standalone `*.sh` templates in `wgs_assemblies/`.
The **pretrained single-sample** binning is a separate PI re-bin in the workbench, **not** part of
this Snakemake.

| stage | tool | version | parameters as run | in → out | mark |
|---|---|---|---|---|---|
| read QC | fastp | **1.0.1** (`envs/fastp.yml`) | `--correction --qualified_quality_phred 20 --unqualified_percent_limit 40 --length_required 50 --complexity_threshold 20`, paired + unpaired outputs, HTML/JSON reports | reads → QC reads | **V** rule `qc_reads` |
| assembly | MEGAHIT | **1.2.9** (`envs/megahit.yml`; seqkit 2.8.2) | `--presets meta-large --min-contig-len 500 --min-count 2 --prune-level 1 --low-local-ratio 0.3`, `-1/-2/-r` (paired + orphan SE), 16 threads, 256 GB | QC reads → **549** assemblies | **V** rule `assembly` + `config.yml` |
| contig rename | rule `rename_contigs` | — | `sed 's/^>[^_]*_\([^ ]*\).*/>{sample}_\1/'` over MEGAHIT `final.contigs.fa` → contigs become **`<sample>_<N>`**; MEGAHIT's `k141_` prefix and all trailing `flag=/multi=/len=` metadata are discarded. The rule then **deletes** `final.contigs.fa` and `intermediate_contigs/` | → `<sample>_contigs.fa` | **V** `rules/assembly/rename_contigs.rule` |
| assembly QC | QUAST (rule `quast`) | **5.2.0** (`envs/quast.yml`) | `quast.py <assembly> -o <dir> --threads N` — **no other flags; all QUAST defaults** | → `assemblies/qc/<sample>/report.txt` | **V** `rules/assembly/qc_assembly.rule` |
| binning B | SemiBin2 | **2.2.1** (conda `semibin2`) | `SemiBin2 single_easy_bin --input-fasta <contigs> --input-bam <bam> --self-supervised -m 1000 --threads 24` (BAM from `run_bwa`) | assemblies → **21,616** HQ+MQ | **V** rule `semibin` |
| binning A | SemiBin2 | **2.2.1** | `minimap2 -a -x sr -t 16 <contigs> <R1> <R2> \| samtools sort` then `SemiBin2 single_easy_bin --environment global --sequencing-type short_read -t 16` — **`-m` not set** | assemblies → **16,333** HQ+MQ | **V** `M/scripts/01_semibin2_single.sh` (PI re-bin, outside the Snakemake) |
| MAG QC input | CheckM2 DB | **`uniref100.KO.1`** — `M/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd`, 3,082,500,605 bytes, dated 2025-02-08 | `--database_path $CHECKM2DB` | — | **V** path resolves and file present |
| refinement | **DAS_Tool** | **1.1.7** | `--bins single.contigs2bin.tsv,self.contigs2bin.tsv --contigs <sample>_contigs.fa --labels single,self --search_engine diamond --proteins NULL --write_bin_evals FALSE --write_bins TRUE` | **542** samples in (`05_refine/sample_list.txt`) → **539** produced bins → **14,426** refined bins | **V** from `05_refine/<sample>/<sample>_DASTool.log`, 2026-06-04 |
| MAG QC | CheckM2 | **1.1.0** (conda `checkm2`) | DB version UNVERIFIED | 14,426 → 4,225 HQ / 8,209 MQ / 1,992 LQ | **V** version; `_checkm2_refined/quality_report.tsv` |

⚠ The Snakemake also contains `metabat2`, `maxbin2`, `concoct`, `vamb`, `magscot` and a `drep` block
with `completeness: 75`. **None of these is the paper's path.** The paper uses the two SemiBin2
binnings, DAS_Tool refinement, and the separate `05_refine/_drep` run at completeness 50 (§5). Do not
transcribe the Snakemake's `drep` config into M.3.

⚠⚠ **The two binnings are NOT symmetric, and M.2 must say so.** They differ in three ways beyond the
model:

| | binning A (pretrained) | binning B (self-supervised) |
|---|---|---|
| model | `--environment global`, no training, no `model.pt` | `--self-supervised`, trains a per-sample `model.pt` |
| coverage mapper | **minimap2** `-x sr`, **paired only, no SE** | **BWA-MEM**, paired **and** SE reads |
| min contig | **`-m` not set** → SemiBin2's adaptive rule resolves to **~2500 bp** on these fragmented soil assemblies | **`-m 1000`** |
| source | `M/scripts/01_semibin2_single.sh` | `EQ_WGS/rules/binning/semibin.rule` |

Consequence: the pretrained track binned only contigs ≥~2500 bp, dropping roughly a quarter of each
assembly that the self-supervised track kept. That is why the two tracks recover different bins and
why DAS_Tool merging them is not redundant.

⚠ `M/scripts/` is superseded as a **live pipeline**, but `01_semibin2_single.sh` is the genuine
historical **provenance** of binning A. Cite it as provenance; do not treat the tree as runnable.

⚠⚠ **TERMINOLOGY — locked decision, do not drift (2026-08-02).** **BOTH** binning runs are
**single-sample**: every sample is binned independently, neither is a co-assembly. The pipeline's
`single_sample` / `self_supervised` labels (file prefixes `ss__` / `selfsup__`, MAG-name tokens
`__single_` / `__self_`) are **misleading**. The real axis is SemiBin2 **pretrained (global) model vs
self-supervised training**. **Always write "pretrained-model vs self-supervised", never "single-sample
vs self-supervised."** M.2 must phrase it that way.

**DAS_Tool funnel — V, both counted on disk:** **542** samples in (`05_refine/sample_list.txt` = the
samples binned by *both* runs) → **539** produced DAS_Tool bins
(`find $M/05_refine -maxdepth 2 -name '*_DASTool_contig2bin.tsv' | wc -l` = 539, and 539
`*_DASTool_bins` dirs). The 3-sample gap is a funnel step, not a discrepancy: 3 DAS_Tool inputs
yielded no refined bins. Related sample-level counts: single **540** ⊆ self **543** (single-only 0,
self-only 3), union **543** = samples contributing HQ+MQ MAGs.

⚠ **Two corrections to the M.2 `\todo` stub in `main.tex`:**
1. **DAS_Tool DID run** — confirmed by per-sample outputs and logs. It merges the two binnings; that
   is what turns the 37,949-bin HQ+MQ union into the 14,426 refined bins. A project memory note
   claiming "both binnings are plain SemiBin2, no DAS_Tool" is **wrong**.
2. **GUNC did not run on the catalogue** — not on the 14,426 refined bins nor the 5,438
   representatives, so M.2 must not describe a catalogue-wide chimerism screen. It **did** run on
   the ten order-gap candidates (GUNC 1.1.1 vs ProGenomes 3,
   `taxonomic_novelty_evidence/legacy_candidates_gunc/`, script `rebuild_20260729/scripts/118_*`):
   all ten pass, CSS 0.01–0.32. The earlier flat "GUNC never ran" came from a `-maxdepth 2` search;
   the artefacts sit deeper. Usable as a rebuttal on the novelty candidates.

**MIMAG quality filter (from `CATALOGUE_NUMBERS.md`, R):** completeness ≥50, contamination ≤10
**inclusive**. Two bins sit at exactly 10.0 — they pass `≤10` but fail strict MIMAG `<10`. **State
the inclusive rule explicitly in M.2/M.3.**

---

## 5. M.3 — Dereplication and taxonomy

| stage | tool | version | parameters as run | in → out | mark |
|---|---|---|---|---|---|
| dereplication | dRep `dereplicate` | **3.6.2** (conda `drep`) | see full table below | **12,433** clustered (of 12,434 passing) → **5,438** reps | **V** all params from `_drep/out/log/cluster_arguments.json`, run 2026-06-07 |

**dRep parameters — V, complete** (`M/05_refine/_drep/out/log/cluster_arguments.json`):

| parameter | value | parameter | value |
|---|---|---|---|
| `S_ani` (secondary, species) | **0.95** | `completeness` | **50.0** |
| `S_algorithm` | **fastANI** | `contamination` | **10.0** |
| `P_ani` (primary, Mash) | 0.9 | `length` (min genome) | 50,000 |
| `cov_thresh` | 0.1 | `MASH_sketch` | 1000 |
| `coverage_method` | larger | `clusterAlg` | average |
| `SkipMash` / `SkipSecondary` | False / False | `processors` | 32 |
| `genomeInfo` | `_drep/genomeInfo.csv` (CheckM2 values supplied; dRep did **not** re-run CheckM) | `run_tertiary_clustering` | False |

Representative-scoring weights: completeness 1 · contamination **5** · N50 0.5 · centrality 1 ·
strain-heterogeneity 1 · size 0. ⚠ `checkM_method = lineage_wf` appears in the argument dump but is
inert — quality came from the supplied `genomeInfo.csv`. Do not report CheckM1/lineage_wf in M.3.
| taxonomy (**current**) | **GTDB-Tk 2.7.2** | reference **r232** | `classify_wf --genome_dir .../legacy_taxonomic_screen_r232/inputs/genomes --out_dir ... --extension fa --cpus 64 --pplacer_cpus 32 --keep_intermediates --write_single_copy_genes` | 5,438 → 51 named phyla | **V** `rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/gtdbtk_2.7.2_r232/run/gtdbtk.log`, 2026-07-31, job **49670976** |
| taxonomy (**superseded**) | GTDB-Tk 2.4.1 | reference **r226** | `classify_wf --genome_dir .../dereplicated_genomes --out_dir .../_gtdbtk_refined --extension fa --cpus 32 --pplacer_cpus 16 --skip_ani_screen` | 5,438 | **V** `05_refine/_gtdbtk_refined/gtdbtk.log`, 2026-06-04 |

**V — 5,438 representatives** counted directly:
`find $M/05_refine/_drep/out/dereplicated_genomes -maxdepth 1 -name '*.fa*' | wc -l` = 5438.

⚠⚠ **The registry cites the wrong path for the R232 numbers.** `CATALOGUE_NUMBERS.md` §1.1 gives
`_gtdbtk_refined/` as the source for "GTDB-Tk 2.7.2 / R232", but that directory is the **r226 / v2.4.1**
run. The manuscript's R232 census actually comes from the PI's rebuild, derived by
`rebuild_20260729/scripts/derive_r232_manuscript_census.py` (driven by
`129_derive_r232_manuscript_census.sbatch`, `--expected-legacy 5438 --expected-lq 304`) from
`coassembly_novelty/legacy_taxonomic_screen_r232/summary`, checksum-gated on
`LEGACY_TAXONOMIC_SCREEN_SUMMARY_COMPLETE` + `taxonomic_screen_summary.sha256`.
**M.3 must cite the R232 path, and the registry source column should be corrected.**

⚠ Note the flag difference between the two runs: the r226 run passed `--skip_ani_screen`; the R232
run did not. They are not the same procedure.

### ⚠⚠ The paper mixes GTDB releases — UNADJUDICATED, and M.10 already asserts otherwise

**§1.1's census is R232. Every downstream taxonomic label is r226.** Eight `_clean` scripts read
`05_refine/_gtdbtk_refined/gtdbtk.{bac120,ar53}.summary.tsv`, which this rewrite established is the
**GTDB-Tk 2.4.1 / r226** run:

| script | section | what it labels |
|---|---|---|
| `darkfam_basis_clean/scripts/fam_9_phylum.sh` | **§1.4** | phylum attribution — 97.9% single-phylum, Actinomycetota 32% / Thermoplasmatota 19% |
| `osmolyte_clean/scripts/osmo_1_census.py`, `osmo_6_turnover.py`, `osmo_7_robustness.py` | **§2.3** | phylum labels incl. the **Halobacteriota** salt-turnover result |
| `uv_survival_clean/scripts/uv_1_census.py` | **§2.3** | phylum labels |
| `resistome_clean/scripts/amr_1_census.py` | **§2.4** | phylum labels |
| `resistome_clean/scripts/amr_3_pathogens.py` | **§2.4** | **genus and species** matching — the pathogen census |
| `bgc_clean/scripts/bgc_1_census.py` | **§2.4** | phylum labels |

⚠ **This is probably harmless at phylum rank but is NOT obviously harmless at genus/species rank.**
**1,898 representatives changed classification between r226 and R232.** The pathogen census
(`amr_3_pathogens.py`) matches **genus and species** names against fixed pathogen panels — exactly the
ranks where reclassification bites — and it is the basis of the "2 named pathogen species of 5,438"
claim.

⚠ **Methods M.10 (now **M.11**) stated the pathogen census "matches GTDB R232 genus and species
assignments of the representatives".** On this evidence that sentence was **wrong as written**: the script
reads the r226 table. **RESOLVED 2026-08-17 (7b): the sentence was corrected to r226** and now points at
M.3's release disclosure. The remaining question is a decision for §2.4's owner, not a gap — re-run the
census against R232 and move §2.4's "2 named pathogen species" number, or retain r226 with the
disclosure.

⚠ Also mislabelled: `darkfam_basis_clean/README.md` describes its phylum attribution as
"R232 at 05_refine/_gtdbtk_refined" — that path is r226. The **numbers** are fine; the **release
label** is not.

**Not adjudicated here.** This is a finding, not a decision — it belongs to Alejandra and the section
owners. **M.3 must state which release each label came from**, and 7b must not paper over it by
writing "R232" throughout.

**First-empty-rank census — how it is derived** (the M.3 stub names it, so M.3 must describe it).
Script `rebuild_20260729/scripts/derive_r232_manuscript_census.py`, driven by
`129_derive_r232_manuscript_census.sbatch` with `--expected-legacy 5438 --expected-lq 304`, reading
`legacy_taxonomic_screen_r232/summary` → `manuscript_census/`. It walks
`RANKS = (domain, phylum, class, order, family, genus, species)` per representative and assigns each
to the **first rank with no named assignment**, reporting depths
`(order, family, genus, species, fully_assigned)`. The classes are therefore **mutually exclusive** —
a genome counted at "genus" is not also counted at "species". Current values: **0/0/0** at
domain/phylum/class · **11** order (6 HQ) · **142** family (64 HQ) · **1,343** genus (530 HQ) ·
**3,755** species (1,555 HQ) · **187** assigned through species. The step is checksum-gated
(`LEGACY_TAXONOMIC_SCREEN_SUMMARY_COMPLETE` + `sha256sum -c taxonomic_screen_summary.sha256`) and
fails closed. **1,898** representatives changed classification between r226 and R232.

⚠ M.3 must state the census is mutually exclusive, or "11 novel orders" reads as cumulative.

⚠ `rebuild_20260729` is **not** blanket-ignorable. The "ignore the rebuild folder" rule applies to
its *coassembly numbers* only. Its `legacy_taxonomic_screen_r232/`, `picrust2_current_20260803/coverm_current/`
and `protein_novelty/references/` are live, paper-facing infrastructure consumed by path.

---

## 6. M.4 — Gene calling, annotation and the catalogues

### 6.1 Assembly protein catalogue (§1.2 basis)

| stage | tool | version | parameters | in → out | mark |
|---|---|---|---|---|---|
| gene calling | pyrodigal | **3.7.1** (conda `pyrodigal`; also stamped `version=pyrodigal.v3.7.1` in the run's own GFF) | `pyrodigal -p meta -i <asm> -a <faa> -o <gff> -f gff -j 4 --no-stop-codon` — **meta mode**, translation **table 11** (per the emitted GFF header), **no min-gene-length**, **closed-ends off** so contig-end genes are kept as `partial=10/01/11` | **552** assemblies → **753,305,533** proteins | **V** verbatim from `M/scripts/ibex/22_pyrodigal_assembly_array.sh`; count **V** in `04_protein_catalog/_pooled_record_count.txt` |
| clustering | MMseqs2 linclust | **18.8cc5c** (conda `mmseqs2`) | **`--min-seq-id 0.9 -c 0.8 --cov-mode 0`** | → **449,265,306** non-redundant reps | **V** (`M/scripts/ibex/26b_linclust_all.sh`, `26_linclust_delta.sh`) |

### 6.2 Sequence darkness — the paper's central definition

Folder `M/04_protein_catalog/darkmatter_clean/`, scripts `darkclean_1` … `darkclean_7`.

| stage | script | tool | parameters as run | in → out | mark |
|---|---|---|---|---|---|
| chunk + search | `darkclean_3a_chunk.sh`, `darkclean_3b_uniref.sh` | `mmseqs search` vs **UniRef50**, MMseqs2 **18.8cc5c** (conda `mmseqs2`) | **`-s 7.5 -e 1e-3 --max-seqs 10 --threads 16 --split-memory-limit 200G`** — pure E-value, **no coverage and no identity filter**; 120-task array | uncovered-rep chunks → per-chunk hits | **V** verbatim from the script |
| combine | `darkclean_4_combine.sh` | — | union of **two** search legs: `delta_vs_uniref50` + `03_uncovered_vs_uniref50` | → `04_catalog_dark/dark_rep_ids.txt.gz`, **69,938,628 dark (15.6%)** | **V** logic |
| length gate | `darkclean_6_dark_lengths.sh` | — | ≥100 aa | → **16,971,287 / 315,606,983 = 5.4%** confident dark | R counts |
| MAG gate | `darkclean_5_mag.sh` | — | `contig_to_bin.tsv.gz` (**both** binnings, unrefined) × `results/catalogue_provenance.tsv` (**37,949** HQ+MQ bins) | → `05_mag_dark/mag_dark_summary.tsv`, **4,292,446 genome-resolved dark reps (9.89%)** | **V** logic |

Query DB: `M/04_structure/uniref50/uniref50`. Search env: conda `mmseqs2`. The script notes it
mirrors the earlier `27_uniref50_search_delta.sh` — same parameters, which is why the two legs are
poolable.

⚠ **The dark call is MEMBER-level, not representative-level.** Verbatim from `darkclean_4`: *"A
cluster is DARK iff it has a searched member and NO hit member."* A cluster with any hit member is
not dark. M.4 must not describe this as "the representative had no hit".

⚠ **The MAG gate is PER-PROTEIN.** Verbatim from `darkclean_5`: *"A cluster is 'MAG-mapped' if ANY
member is on an HQ+MQ-binned contig"*, checked on **every** member — because a cluster's members sit
on different contigs, some binned and some not, so the representative alone is not representative.

⚠ **"MAG-mapped" means protein-to-contig MEMBERSHIP, not read mapping.** This is the one place in
the paper where "mapped" does not mean reads. M.4 must not phrase it ambiguously. Rationale recorded
in the script: the 753M assembly proteins include unbinned/fragment/noise contigs, so only proteins
on a binned HQ+MQ contig are genome-resolved.

⚠ The 4,292,446 genome-resolved dark reps are the **single shared input object** for §1.3
(structural) and §1.4 (families). The gate exists deliberately, to exclude unbinned noise.

### 6.3 Per-representative KO catalogue — `rep_ko_v2.tsv`

Folder `M/04_functional/repko_clean/`. **STATUS: DONE 2026-08-12.**

| step | script | tool | version + parameters | in → out | job |
|---|---|---|---|---|---|
| 1a annotate | `repko_1_annotate.sh` | Bakta | **1.12.0** (conda `bakta`), DB **`/ibex/reference/KSL/bakta/6.0/db`** — `bakta --db $BAKTA_DB --output … --prefix … --threads 16 --skip-plot --force` | the **124** never-annotated raw bins | **50376640** (124/124) ✅ |
| 1b orthology | `repko_1_annotate.sh` | eggNOG-mapper | **2.1.13** (conda `eggnog`), DIAMOND **2.0.15**, data dir **`/ibex/reference/KSL/eggnog-mapper/2.1.12`** — `emapper.py -i <faa> --output … --output_dir … --data_dir $EGGNOG_DB --cpu 16 --itype proteins --override` | → `KEGG_ko` assignments | **50378049** (123/124) + resubmit **50379295** ✅ |
| 2 submap | `repko_2_submap.py` | — | contig-subset construction for `_sub` reps | **1,010/1,010** corrected | ✅ |
| 3 build | `repko_3_build.py` | — | assemble final table | → **`results/rep_ko_v2.tsv`** | ✅ |

**Result: 5,438 rows · 0 unannotated (was 168) · 4,580 identical to the original · 690 reps shed at
least one spurious KO (8,134 total).** — **V** from the README STATUS block.

⚠ **emapper inherited defaults.** No `--sensmode`, `--tax_scope`, `-m` or e-value was passed, so
eggNOG-mapper 2.1.13's defaults applied (DIAMOND search, `sensitive` mode, automatic taxonomic
scope). M.4 should say the defaults were used rather than imply tuned settings. Note the tool is
**2.1.13** while the reference data directory is named **2.1.12** — that is the on-disk state; quote
both, do not silently harmonise them.

⚠ **Use `rep_ko_v2.tsv`. Never the PI's `rep_ko.tsv`** (`M/04_functional/refs/rep_ko.tsv`, produced by
`scripts/eco_probes/gather_rep_ko.py`, hohndor 2026-06-11). The original is **never modified** by
policy; both files exist side by side. §2.1–§2.5 all use v2. Consequence of the fix: §2.2 moved
1,822→1,860 guild reps and 18.0%→18.4% abundance.

**Three defects fixed — V, from the full 5,438 reconciliation:**

| reps | state | cause |
|---|---|---|
| 4,304 | annotated, contig-identical to raw bin | ✅ correct, untouched |
| **966** | annotated `_sub` reps that inherited the **PARENT's** annotation | ⚠ superset — 915 retain ≥95% of parent, **51 retain <95%** (worst 25%). Measured on one rep: 422 of 1,521 KOs (28%) spurious |
| **44** | unannotated — `48PRr1__SemiBin_N`, no `single_`/`self_` token | ❌ the resolver splits on that token and fails; **the annotation existed all along** |
| **80** | unannotated — raw bin was LQ pre-refinement, never on the HQ+MQ list | ❌ refinement promoted them to HQ/MQ afterwards; nothing re-annotated them |
| **44** | unannotated `_sub` whose parent was also never annotated | ❌ parent-fallback had nothing to fall back to |

**Root cause:** annotation ran on the **raw, pre-refinement** HQ+MQ bins
(`03_annotation/mag_list.txt`, 15,469 rows) while the catalogue is built from **post-refinement**
bins — refinement both created new `_sub` bins and promoted LQ bins. The two 44-rep sets are
disjoint (0 overlap); 168 unannotated total.

⚠ Note the 44 `48PRr1` reps are the *same* naming exception documented in §11 — the plain-token MAG
name is not cosmetic, it silently broke a join.

⚠ **KEY LESSON for any future annotation job:** eggNOG-mapper needs **`--dbmem`** on WekaFS. Without
it, emapper random-reads the 41 GB `eggnog.db` and concurrent tasks deadlock in D-state on
`wekafs_lock_helper`.

**Gene-count denominators — not interchangeable:** **19,186,853** protein-coding genes (Bakta calls,
the proteome-consistent denominator, → 55.6% no-KO ceiling) vs **20,445,320** (CheckM2's own calls,
~6% higher, → 58.4%). Quote 55.6% and say which denominator.

---

## 7. M.5 — Abundance estimation  ⚠ THE INHERITED-DEFAULT SECTION

M.5 is currently an unfilled `\todo` stub yet it governs every community number in §2.1–2.4. It must
define "detected" once, and that definition is **not** a chosen threshold.

### 7.1 Rep-level run (PI-owned; §2.1–2.4)

Path: `M/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/`
— **PI-owned, no `_clean` equivalent, consumed by path. Do not rebuild.**

**The command, literally — V**, from `rebuild_20260729/scripts/run_current_coverm_profile.py:193-210`
and confirmed byte-for-byte against a recorded `run_provenance.json` `/command` entry:

```
coverm genome \
  --reference               <prepared>/reps.mmi \
  --minimap2-reference-is-index \
  --genome-definition       <prepared>/genome_def.tsv \
  --coupled                 <sample>_R1.fastq.gz <sample>_R2.fastq.gz \
  -m relative_abundance -p minimap2-sr -t 16
```

| field | value |
|---|---|
| tool | **CoverM 0.7.0** (`envs/eco_probes/bin/coverm`) — **V**, conda metadata |
| reference | `prepared/reps.mmi` (61.9 GB index) from `04_functional/refs/reps.fa` (20.94 GB) = the **5,438** representatives |
| genome definition | `prepared/genome_def.tsv` (contig → genome map) |
| cohort | **569 profiles, 548 usable** — see the count reconciliation below |

⚠ **That is the whole command. Everything else is a CoverM default.** Four consequential parameters
were **never passed**:

| parameter not passed | inherited default | consequence |
|---|---|---|
| **`--min-covered-fraction`** | **10%** (genome mode) | **defines "detected"** — see below |
| `--contig-end-exclusion` | 75 bp | ignores contig termini when computing coverage |
| `--min-read-percent-identity` | 0 (off) | no identity filter on alignments |
| `--min-read-aligned-percent` | 0 (off) | no alignment-length filter |

⚠⚠ **Write this into M.5 explicitly: "detected" = relative abundance > 0, which under the unset
`--min-covered-fraction` means ≥10% of the genome covered.** It is an **inherited default, not a
chosen threshold**, and must be stated as such. Point M.8 and M.11 at that one definition so
§2.1–2.4's absence language becomes a defined term. The same applies to the 75 bp end-exclusion —
report it as a default, not a choice.

**Unmapped-fraction handling — M.5 must state this** (the stub asks for it explicitly). Median
unmapped read fraction is **~51%**, and it is **balanced between Surface and Deep** (54% in both), so
paired compartment contrasts are not artefacts of catalogue coverage. Two denominators are used
throughout Part 2 and are **not interchangeable**:

| denominator | meaning | used for |
|---|---|---|
| **whole-community** | share of *all* reads; unmapped organisms treated as non-carriers | the conservative bound |
| **mapped-community** | share of the mapped fraction only | comparisons within what the catalogue represents |

M.8 and M.10 already use both (e.g. AMR carriage median 15% whole vs 33.5% mapped; BGC 48% vs 99%).
Root-adjacent comparisons use the mapped denominator. M.5 should define the pair once and let M.7–M.11
refer to it.

**Profile count reconciliation — V, all three counted on disk.** Say **"548 usable of 569 profiled"**;
573 is a directory count, not a profile count.

| count | what it is |
|---|---|
| **573** | directories under `profiles_full_recompute/` — includes **4 abandoned `.tmp` staging dirs** |
| **569** | directories carrying a `run_provenance.json` = **the real profile count** |
| **548** | usable canonical soil cohort quoted by §2.1–2.4 |

The four staging dirs are exactly the four failed CoverM array tasks:
`.51Sr1.tmp.50030409.0`, `.F5Dr1.tmp.50063405.0`, `.F13PRr1.tmp.50051922.0`,
`.SPositive.tmp.50089048.0` (the last is the positive control, not a soil sample). 573 − 4 = 569. ⚠ The
paper-facing three failures are `51Sr1`, `F5Dr1`, `F13PRr1`; `SPositive` is a control and was never in
the canonical cohort.

✅ **All 569 profiles are FRESH recomputes.** Verified:
`grep -ho '"profile_source"' … | sort | uniq -c` → **569 ×
`new_prebuilt_index_mapping_full_cohort_after_equivalence_rejection`**, zero reused. The source string
records *why*: the legacy-equivalence check was **rejected**, which is precisely why the whole cohort
was recomputed.

⚠ The script does contain a `validated_legacy_exact_reference_reuse` route (lines 181–191) that copies
a legacy `abund.tsv`/`coverm.log` on a reference match and fails closed on sha256 mismatch — but it is
an **unexercised code path for this cohort**, not a caveat on the numbers. Do not describe these
profiles as a mix of fresh and reused runs; that would imply a provenance weakness the data does not
have.

Downstream detection constants: `DETECT = 0.0` in `{osmo,uv,amr,bgc}_2_community.py`, with a
`DETECT2 = 0.01` sensitivity variant alongside; `MIN_SPP = 15` in `osmo_4_phylum.py:65`.
⚠ Corrected 2026-09-17: this previously cited `refbias_2_magko.py:82` for a `covf >= 0.5` presence
rule. No coverage-fraction test exists anywhere in `refbias_clean/`; line 82 is `elif v > 0:` and
presence there is CoverM's never-changed 10 % breadth default.

### 7.2 Member-aware run (§1.4, §2.5)

| field | value |
|---|---|
| script | `darkfam_basis_clean/scripts/fam_4_coverm_array.sh` |
| command | `coverm ... -m tpm covered_fraction`, 551-task array |
| reference | the **536,437** dark-family member CDS |
| presence rule | **any member gene covered fraction ≥ 0.5** (`fam_8_aggregate.sh`) |
| family abundance | **sum of member TPM** |
| cohort | **551** (complete) |
| job | **50440132**; aggregation **50480972**, refreshed **50531118** |

⚠ **The `covf ≥ 0.5` presence rule is length-biased — demonstrated, not inferred.** Median
representative length falls from 145 aa at 1–3 sites to **37 aa** at 56–61 sites; swapping the rule
for `TPM>0` reverses the relation (76 → 175 aa). **Any distribution claim from `fam_presence.tsv`
must be length-stratified.**

**§1.4 is the only affected claim and it survives.** Stratified <60 aa → ≥250 aa: median **201 → 55**
of 551 samples; all-five-expedition persistence **89.7 → 67.1%**; ≥4 expeditions ≥90.8% in every
stratum; all four single-sample families are <60 aa, none ≥60 aa. So `98/551`, `4` and
`79.4%/94.8%` **stand as computed** — 98 is composition-weighted and spans 201–55, so quote it with
that floor. The anti-relic-DNA inference is *strengthened*.

**Unaffected:** §2.5 gates on `TPM>0`; §1.4's other legs are assembly/annotation-based; §2.1–2.4 use
genome-level mapping — tested, median samples detected by genome size is 7/35/43/44/42 from <1.5 Mb
to ≥5 Mb, flat above 2.5 Mb and *lower* for small genomes, so there is no size analogue.
Audit scripts: `specificity_clean/scripts/{presence_rule_length_test,s14_persistence_by_length,mag_detection_size_test}.py`.

---

## 8. M.6 — Structural analysis (§1.3)

Folder `M/04_structure/00_structure_clean/`. **Complete** (2026-08-18, job 50685374): K = 1,031 of
2,838 confident models, control floor 101 of 3,013. The multiplier is now reported standardised on
model confidence (`05_validate/stratified_floor.tsv`, 2026-09-17), not crude. The earlier
"stage 05c still running (50603810)" header was stale and is removed.

| step | script | tool + version | parameters as run | in → out | job |
|---|---|---|---|---|---|
| 01 prep | `01_prep.sh` | — | **no length gate here** — splits the 4,292,446 MAG-dark reps into **2,898,693** with a reusable tier1/tier2 3Di and **1,393,753** needing prediction | → `reuse_ids.txt` / `missing_ids.txt` | **V** counted on disk |
| 02 3Di | `02_prostt5.sh` | ProstT5 via **FoldSeek GPU build `8dc75c74`** | `foldseek createdb <chunk.fa> <queryDB> --prostt5-model $STRUCT/prostt5_weights/prostt5-f16.gguf --gpu 1 --threads 8` — **no batch or length flags**; this `createdb` IS the 3Di generation, there is no separate step | → `queryDB` + `queryDB_ss` | **V** verbatim |
| 03a/b search | `03a_prep_query.sh` (**this is where the ≥100 aa gate lives**: MAG-dark ∩ len≥100 → **1,466,943** `target_ids.txt`, **V**), `03b_prep_chunks.sh`, `03b_search.sh` | FoldSeek **GPU build `8dc75c74ad0eddab73cfd905963d13bf74dc012b`** (`M/tools/foldseek/bin/foldseek`) — **NOT** 10-941cd33; an earlier version of this guide had that wrong | `--alignment-type 2 -s 9.5 -e 1e-3 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G` | 215 chunks of ~7k vs AFDB v6 + PDB | **50342556** (207/215) + **50542223** (8 requeued) — **CLOSED 215/215** |
| 03c aggregate | `03c_candidates.sh` | — | retrieval decomposition | → hit/no-hit tables | **50559990** ✅ |
| 04 cluster | `04_cluster.sh` | **FoldSeek `cluster`** (GPU build `8dc75c74`) — **NOT MMseqs2**; an earlier version of this guide had that wrong | `foldseek cluster <merged> <clu> <tmp> --single-step-clustering --cluster-mode 2 -s 4` (**no `--min-seq-id`, no `-c`, no `--cov-mode`**), then admission **≥3 members** applied in shell on `cluster_sizes.txt` (`awk '$1>=3'`), not by a tool flag; each member is already a 90%-id cluster rep | 997,646 no-hit → **901,638 clusters** → **17,221 families** | **50560204** ✅ |
| 05a rep prep | `05a_prep.sh` | — | re-shard for backfill | → 28 shards | 50560375, then **50573919** ✅ |
| 05b fold | `05b_fold.sh` (A100), `05b_fold_v100.sh` (V100) | ESMFold `facebook/esmfold_v1` via **HuggingFace `transformers` 4.40.2 / torch 2.2.2** — **not `fair-esm`** | half-precision ESM trunk, `set_chunk_size(64)` A100 / `(32)` V100, sequences **truncated to 1,022 aa**, **no `num_recycles`**; pLDDT = mean of `plddt[0,:,1]` written to `plddt_<shard>.tsv`; 28 shards | **13,770 folded** + **3,451 reused** = **17,221** models | **50573938** (1–14) + **50589898** (15–28), 28/28 ✅ |
| 05c screen | `05c_screen.sh` | FoldSeek **10-941cd33** (`module load`) — the **only** step on the packaged build, chosen to reproduce the PI's `167_postfold` invocation exactly | `--alignment-type 1 -e 10 -s 9.5 --max-seqs 50` | 20,905 queries = 13,770 folds + 3,451 reuse + **3,684 controls** | **50603810 RUNNING** |

⚠⚠ **TWO FoldSeek BUILDS, deliberately — M.6 must say so.** Steps 02, 03a, 03b and 04 all run
`$M/tools/foldseek/bin/foldseek`, version **`8dc75c74ad0eddab73cfd905963d13bf74dc012b`** (a git-commit
build id, not a release tag), because `02_prostt5.sh` records that the conda/module foldseek **crashes on
GPU**. Only **05c** uses `module load foldseek/10-941cd33`, to reproduce the PI's calibrated `167_postfold`
invocation byte-for-byte. This is defensible and is what M.6 states: the two stages are separate
measurements (predicted-3Di retrieval vs full-coordinate validation), and inside 05c the targets and the
3,684 controls pass through the **identical** build and command line, so the control floor is build-matched
to the thing it calibrates. Verified 2026-08-17 by reading `FS=`/`module load` in every script.

**Screen result (banked at 03c):** 1,466,943 searched → AFDB-any 463,316 · PDB-any 232,793 · dual
226,812 · AFDB-only 236,504 · PDB-only 5,981 · **no structural hit 997,646 = 68.0%**. The 851,393
singleton and 33,024 pair clusters carry **no claim**.

**Reference databases — resolved; the June guide was wrong.**

| DB | release | entries | path |
|---|---|---|---|
| AFDB | **v6 / UniProt 2025_03** | 241,070,489 | `rebuild_20260729/structural_calibration/references/afdb_uniprot_v6_uniprot_2025_03/database/afdb` |
| PDB (03b screen) | **250101** | 324,204 | `M/04_structure/fsdb_pdb` |
| PDB (05c validation) | **250101** — same snapshot | — | `references/pdb_20260730/foldseek_pdb` |

⚠ `afdb_pad` (AFDB **v4**) is used **nowhere** in this chain; the June guide pointed at it and was
wrong. ⚠ The two PDB paths are provably the **same data** — identical `.version` md5
(`aefd75e0a5d6acbe8a2b7791b53eb479`, `pdb100.tar.gz`) and `PDB_DATE 250101`. `pdb_20260730` is the
reference **freeze date, not a PDB release**. **M.6 states one PDB release for the whole chain.**

**Gates and call rule — M.6 must state all four:** length **≥100 aa** AND mean pLDDT **≥0.70** (of
the 13,770 folded, mean pLDDT 0.535; **2,452 = 17.8%** clear the gate); primary call = dual
**query-normalised** TM **< 0.5** against both references; alignment-normalised TM reported as a
**sensitivity companion** — never one normalisation alone.

**Known-fold controls.** **3,684** blinded controls (`KF000001…`) built by the PI: real PDB chain
sequences or deterministic crops simulating Prodigal partial genes (partiality models 00/01/10/11),
length- and partiality-matched to dual-novel targets, ESMFold-folded into `control_folds/shard_001..032`;
manifest `current_reference_control_manifest_20260802.tsv` with per-model sha256. All are
`fold_status=ok` and **co-queries in 05c** under the byte-identical invocation.

⚠ **The control floor M.6 states must be OURS, pending 50603810** — not the PI's 11.9% overall /
3.4% confidence-qualified (101/3,013, Wilson 95% 2.8–4.1%), which are on his query set.

⚠ **Never quote as ours:** the PI's **608** (old May tier-1 basis, 0.87% query overlap); and
**12,751**, which is **not a model count** but the ID-level overlap between his 68,527 candidates and
our 1,466,943 queries.

**Do not cite:** 50253325 (cancelled) · 50560000/50560011/50560036 (failed, script bugs) ·
50560586/50574010 (cancelled before starting, never ran) · 50603505 (ENOSPC) · 50603687 (duplicate).

⚠ **Job-ID drift between sources.** The workbench README STATUS block is dated 2026-08-16 and lists
the second 05b array as `50574010` while 05b was still running. That array was cancelled before
starting and re-issued as **50589898**; 05b then completed 28/28. **Use 50573938 + 50589898.** If the
README still shows 50574010 when 7b reads it, the README is the older snapshot — session 06 owns it
and is updating it as 05c lands.

⚠ **The 03c decomposition is a predicted-3Di screen and overstates novelty by design** — that is what
the ESMFold validation bounds. The 68.0% is not the paper's novelty figure; the calibrated,
control-floored number from 05c is. Do not quote 68.0% as the result.

---

## 9. Downstream analyses

### §1.4 Dark families — `M/06_darkmatter/darkfam_basis_clean/`

| step | script | tool | parameters | job |
|---|---|---|---|---|
| cluster | `fam_1_cluster.sh` | MMseqs2 easy-cluster | **`--min-seq-id 0.5 -c 0.8 --cov-mode 0`** | **50421060** |
| admission | `fam_2_recurrence.sh` | — | **≥3 diverged members AND ≥3 canonical samples AND ≥2 sites** | **50439896** |
| member CDS | `fam_3_cds.sh` | — | → `darkfam_members.fna` **536,437** | **50435064** |
| abundance | `fam_4_coverm_array.sh` | CoverM 0.7.0 | `-m tpm covered_fraction`, 551 tasks | **50440132** |
| control | `fam_5_control.sh` | MMseqs2 | HQ+MQ gate, full-dark exclusion, **length-stratified** sampling | **50439953** |
| dN/dS | `fam_6c_dnds.sh` | mafft + **`Bio.codonalign.cal_dn_ds`** (Biopython NG86, **not** codeml) | families with ≥3 member CDS; drops CDS <90 nt, proteins <25 aa or with internal stop, families left with <3 usable seqs, pairs with <30 usable codon sites | **50480238** (dark), **50480461** (control) |
| robustness | `fam_7_robustness.sh` | — | 30% recluster | **50480205** |
| aggregate | `fam_8_aggregate.sh` | — | presence `covf≥0.5`, family TPM = Σ members | **50480972** → **50531118** |
| phylum | `fam_9_phylum.sh` | — | contig→bin→species join | — |

**Universe: 106,042 families / 536,437 member CDS.** Excluded: 2,784 paralog-only, 782,616 young
(1–2 variants, no claims). **SETTLED, do not re-litigate:** the universe definition, the matched
control as the only defensible dN/dS comparison, the presence rule. B4/FUBAR retired.

### §2.1 Reference bias — `M/06_reference_bias/refbias_clean/`

PICRUSt2 **v2.4.1** predictions, checksum-pinned in `picrust2_current_20260803/inputs/` (10,543 KOs
× 1,270 profiles + per-profile NSTI, `INPUTS_COMPLETE`; `ko_metagenome_unstrat.tsv` sha256
`c71c2256…`, byte-identical to the companion amplicon repo's pinned input). Jobs **50586647** +
**50586699**, all 5,438 reps resolved, 0 failures. Ceiling from `rep_ko_v2`. Methods **M.11 is
already written.**

### §2.2 Trace gases — `M/06_atmospheric/lifeonair_clean/`

On 5,438 reps / `rep_ko_v2`; abundance over the 548 cohort. Methods **M.7 is already written.**

| step | script | tool | parameters as run | mark |
|---|---|---|---|---|
| 1 crude census | `lifeonair_1_census.sh` | awk over `repko_clean/results/rep_ko_v2.tsv` | KO presence: **hydB = K06281**, **coxL = K03520** | **V** |
| 2 typing | `lifeonair_2_type_guild.sh` | **DIAMOND** 2.0.15 (`$CONDA/eggnog/bin/diamond`) | `blastp -d $REFS/hyddb --ultra-sensitive --max-target-seqs 1 --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore` | **V** |
| 2 trees | same | MAFFT + FastTree (`$CONDA/busco/bin/`) | → `classify/hydB.nwk`, `coxL.nwk` | **V** |
| 3–4 | `lifeonair_3_abundance.py`, `lifeonair_4_energetics.py` | — | abundance over 548; energetic balance | — |

**V — the tool is DIAMOND, not HMMER**, confirming the June guide on this point. Step 1 is explicitly
a *"CRUDE KO CENSUS … FIRST PASS only: KO presence over-counts, because it does not distinguish"* the
groups — the 1h/1l and form-I calls come from the step-2 DIAMOND typing against HydDB, not from the
KO census. **M.7 must not describe the guild counts as a KO census.** Marker proteins are extracted
from the per-sample annotation tree at
`P/processed/sequencing-results/sites/<trip>/site<N>/<sample>/wgs/annotation/`. Nothing is
re-annotated; only the typing is new compute.

### §2.3 / §2.4 Stress, BGC, resistome — `M/06_metabolism/*_clean/`

| folder | tool | version | parameters as run | mark |
|---|---|---|---|---|
| `bgc_clean` | antiSMASH | **8.0.4** (`module load antismash/8.0.4`) | `antismash --minimal --genefinding-tool prodigal --cpus 8 --output-dir <mag> --output-basename <mag> <fna>`; array `1-272%60`, BATCH=20, 8 CPU / 24 G / 10 h. **`--minimal` = core cluster detection only** (no downstream comparative modules) | **V** `bgc_antismash_array.sh` |
| `resistome_clean` | AMRFinderPlus | **4.2.7**, DB **2026-05-15.1** (conda `amrfinder`) | `amrfinder -p <faa> --plus --name <mag> --threads 4`; array `1-37`, BATCH=150, 4 CPU / 8 G / 6 h | **V** `amr_array.sh` + `results/VERSIONS.txt` |
| `resistome_clean` | VFDB setB via MMseqs2 | — | id ≥40%, cov ≥70%, e ≤1e-10 | R — PI compute, reused |
| `uv_survival_clean` | KO panel on `rep_ko_v2` | — | **28 KOs across 8 functional groups** = 29 rows of `06_metabolism/uv_survival/uv_markers.tsv` minus `K01859`, excluded in code | **V** `uv_1_census.py` |
| `osmolyte_clean` | KO panel on `rep_ko_v2` | — | **13 KOs** (inline `KO_GENE` dict) | **V** `osmo_1_census.py` |
| `pathway_check_clean` | triangulation | — | antiSMASH × KO Jaccard **0.921**; trehalose KO × gapsmith **0.795** | R |

⚠ **`--plus` is passed; `--organism` is NOT.** `--plus` is what yields the STRESS/metal layer that
§2.4 reports (arsenic, mercury) alongside the AMR classes — the type column separates them. Omitting
`--organism` means **no organism-specific point-mutation detection and no organism blacklisting**;
resistance is called from acquired-gene hits only. State that in M.10 — it is an inherited default,
not a choice.

⚠ **`K01859` is deliberately excluded from the UV panel.** It is labelled "phrA
photolyase/cryptochrome" in `uv_markers.tsv` but is chalcone isomerase (EC 5.5.1.6) per
`kofam_ko_list`. Verified 2026-08-13: 1 carrier in 5,438, which also carries `K01669` → zero effect
on any count. The 28-KO panel size **confirms** the "28-KO" figure already written in M.8.

⚠ Both census scripts read taxonomy from `05_refine/_gtdbtk_refined/gtdbtk.bac120.summary.tsv`
(`osmo_1_census.py:27`, `uv_1_census.py:28`). Per §5 that path is **r226**, not R232 — so §2.3's
phylum labels rest on r226. Harmless at phylum rank if the labels agree, but **M.3 must say which
release each downstream label came from.** Not yet adjudicated.

Statistics: partial Spearman controlling ranked lat/lon, Benjamini–Hochberg **within** the salt and
climate blocks; n=248 chemistry, n=453 climate. Methods **M.8 and M.10 are already written** — they
are the better source for these two sections than this guide.

### §2.5 Dark-matter ecology — `M/06_darkmatter/darkfam_ecology_clean/`

`eco_{1,2,3}*.py` read **`fam_abundance.tsv`** (not the presence matrix) and gate on **non-zero TPM in
≥25 subset samples** (`MIN_PRESENT = 25`). PERMANOVA **999** permutations; Borcard variation
partitioning **499**. NumPy **1.26** / SciPy **1.11**. Job **50575481**. Shipped 32e24ae / 4cdfcab.
⚠ The gate was described three ways; the correct wording is **"non-zero abundance"**.

### §3.1 Specificity — `M/06_darkmatter/specificity_clean/`  ✅ compute CLOSED

All three legs complete. **Numbers are registered by session 05 in `CATALOGUE_NUMBERS.md` §3.1 —
that registry is authoritative; the values below are here for parameter context only.**

| leg | tool | parameters as run | scope | job |
|---|---|---|---|---|
| global/deserts | **DIAMOND v2.0.15** | `blastp --sensitive --evalue 1e-3 -k 1 -f 6 qseqid sseqid pident evalue` | mgnify90 717,738,164 · gem 122,978,191 · spire 269,496,580 · deserts_combined 147,531,178 | 50594992/50594993 |
| rarefaction | DIAMOND + `scripts/raref_curve.py` | one pass, `blastp --sensitive --evalue 1 -k 500 --max-hsps 1` vs mgnify90; curve computed **analytically** from it | size-matched expectation per catalogue | **50613400** ✅ 3h52m |
| depth | MMseqs2 | `search -s 7.5 -e 1e-3 --max-seqs 10` vs full `mgy90` **2026_07** | 1,661,142,385 seqs | **50618838** ✅ 1h37m |

Query = **212,084** seqs in one pass: 106,042 dark family reps **plus** 106,042 length/size-matched
annotated controls. deserts_combined = atacama 6,193,359 + desert_autotroph 7,200,701 + arid
30,569,562 + antarctic_colddesert 49,206,550 + mojave 54,361,006.

⚠ **The comparators are three different KINDS and M.6/M.12 must say so.** The five desert DBs are
SPIRE per-study **gene calls** (`SAMN*.psa_megahit.prodigal.faa.gz` — assembly-level, binned or not).
The SPIRE DB is `spire_representative_genomes/specI_v4_*.fa.gz` — dereplicated **representative
genomes**, gene-called locally with prodigal `-p meta` into 107,008 per-genome parts, then
`diamond makedb` streamed through a FIFO (index 121.8 GB at `/ibex/user/velazqal/eqmags_31/spire/`).

⚠ An e-value/DB-size confound was **measured**: 2,544 families (2.40%) hit only because the DBs were
split. Rarefaction is therefore **e-value-corrected** and **no subsampled databases were built** —
deliberately, because building them would reintroduce the very size-dependent threshold being
corrected for. The curve is analytic, derived from a single DIAMOND pass.

**Rarefaction — what it licenses, for M.12.** The curve gives, for any database size N, the absence a
*random N-sequence slice of MGnify90* would produce, so each real catalogue is scored against its
size-matched expectation. This turns the cross-catalogue ordering from a descriptive remark into a
claim, so **word it as a claim, do not merely swap numbers in**.

⚠⚠ **Get the sign right — this is the sentence M.12 will transcribe, and it is easy to invert.**
**Gap = expected − observed.** A **positive** gap means observed absence is *below* the size-matched
expectation: the catalogue leaves **less** of the dark proteome unexplained than its size alone
predicts, i.e. it contributes relevant content beyond raw depth. Worked example:
`deserts_combined` observed 42.60 against an expected 66.93 → **+24.33**, so the desert catalogues
explain 24 points *more* of our dark proteome than a size-matched random MGnify90 slice would.
A **negative** gap is the reverse — `atacama` observed 95.61 against expected 90.94 → **−4.67**,
leaving *more* unexplained than size predicts.
**Do not write "a positive gap explains more absence"** — that states the opposite of the result.

⚠ **The `mgnify90` self-validation is what licenses the rest and should lead the description:** at
p = 1 the model reproduces the observation to **0.15 points** (52.68 observed vs 52.53 expected).

⚠ **`-k 500` cap:** its bias direction **flipped** under the gap framing — truncation now inflates the
expectation, hence the gap. Magnitude is negligible ((1−p)^500 ≈ 1.3% at the smallest anchor). Any
older docstring saying otherwise is stale.

⚠ **Grid-anchor footnote:** the `spire` expectation is read at the curve's **271,000,000** grid anchor
while the real index holds **269,496,580** sequences. The difference is ~0.1 points and changes no
conclusion, but if the gap table is reproduced in the paper, either footnote the anchor or relabel it
so the arithmetic is internally exact — a referee checking the numbers will notice.

⚠⚠ **`atacama` must be flagged wherever the gap is tabulated.** Its *control* gap is **−29.97**: a
single 35-sample study is a narrow, deeply redundant slice of protein space and is worse than a random
draw even for conserved proteins. The gap metric assumes a random draw and one narrow study violates
that badly. Read against its own −30 control penalty, its −4.67 dark gap is actually favourable.
**The gap is interpretable for the large, broad catalogues only.**

**Depth leg — a robustness check, not a table row.** Of the 30,479 dark families with no DIAMOND hit
anywhere, only **1,760 (5.8%)** gain a hit against the complete 1.66 G-sequence MGnify90 2026_07;
**28,719 remain unmatched = 27.08%** of the 106,042 universe (28.74% before). Control: 435 → 149 hit →
286 unmatched = 0.27%. The asymmetry is the point: the deeper reference resolves **34.3%** of the
control's residual absences but only **5.8%** of the dark set's.
⚠ Keep the **MMseqs2-vs-DIAMOND** caveat — different tool from the primary legs. This is a robustness
leg, **never a `tab:specificity` row, and never unioned with the DIAMOND absences.**

---

## 10. Dead code and superseded objects

| object | status |
|---|---|
| `M/scripts/` | **superseded PI code** — this tree is what made the June guide stale. Not live pipeline. |
| `darkfam_basis_clean/scripts/basis_{1,2,3}*.sh` | **never run**; the strict re-search was replaced by the rebuild. Kept for the record. |
| June **23,382**-family dark object | superseded by the 106,042 universe; banner-marked, never deleted |
| PI's `rep_ko.tsv` | superseded by **`rep_ko_v2.tsv`**; never overwritten, both exist |
| `afdb_pad` (AFDB v4) | not used anywhere in the current structural chain |
| `05_refine/_gtdbtk_refined` | **r226 / v2.4.1** — superseded for taxonomy by the R232 screen |
| `rebuild_20260729` coassembly **numbers** | ignored — but the folder is **not** blanket-dead (see §5) |

**Superseded → current** (full table in `CATALOGUE_NUMBERS.md`):

| superseded | current | domain |
|---|---|---|
| union-dRep 9,253 species · R226 · 56 phyla | **5,438** · R232 · **51 phyla** | taxonomy |
| MAG-Bakta proteome 51.9M → 2.84M dark | 449M reps → 69.9M dark; **4.29M** genome-resolved | proteome |
| single_sample 15,469 / 15,477 (strict `<10`) | **16,333** (MIMAG `≤10`) | binning tiers |
| structural "89% novel" / "10,007 folds" | calibrated search + controls (§1.3) | structural |
| June 23,382 · dN/dS 0.391 · 52/551 | **106,042** · 0.425 vs 0.388 · 98/551 | dark families |
| C1 on 15,469 / 436: 56.7% · recall 93% | **55.6%** · **86.0%** on 5,438 / 448 | reference bias |
| B1a structural calibration 608/1,875 | superseded — **never quote 608** | structural |

---

## 11. Naming conventions

| pattern | binning | catalogue |
|---|---|---|
| `<sample>__single_SemiBin_N` | single-sample (pretrained) | **current** 5,438 line |
| `<sample>__self_SemiBin_N` | self-supervised | **current** 5,438 line |
| `<sample>__SemiBin_N` (plain) | older single-sample run | **old** 15,469 line |

⚠ **Verified exception:** the 44 representatives of sample `48PRr1` are plain `48PRr1__SemiBin_N`
yet **are** current catalogue — that sample has only self-supervised binning, so no token was added.
Rule: plain-named **and** `48PRr1` = current; plain-named and any other sample = old. Ground truth:
`M/05_refine/_drep/out/dereplicated_genomes/`.

`_sub` (e.g. `25PRr2__self_SemiBin_2_sub`) = contig-subset representative; its KO calls come from
`repko_2_submap.py`'s contig-subset construction, not the parent superset. One filter-passing bin was
dropped inside dRep (`25PRr2__self_SemiBin_2_sub`), which is why 12,433 is the count `Cdb.csv`
membership analyses see.

---

## 12. Environment and operational caveats

- ⚠ **The project filesystem is full and it has already killed a job.**
  `/ibex/scratch/projects/c2014` (= the `/ibex/project/c2014` symlink) sits at **144 TB of 145 TB,
  ~1.6 TB free, shared**. Job **50603505** (structural 05c) died on **ENOSPC** — FoldSeek could not
  write AFDB temp files. **The working pattern, already in `05c_screen.sh`:** point job tmp at
  `/ibex/user/velazqal` (1.5 TB, 5% used) and move only compact final outputs back. Any new heavy
  stage must do the same. Project-level capacity problem, on the PI list — the whole structural
  pipeline is only ~6 GB.
- `ls` times out on directories with >100k files. Use `find <dir> -maxdepth 1 | wc -l`. A broad
  `find` across the whole workbench also times out — go per-directory.
- Conda envs live at `/ibex/scratch/projects/c2014/rob/miniconda3/envs/`. Versions in this guide were
  read from `conda-meta/`, which is authoritative; README prose is not.

**Tool version quick reference (all V from conda-meta or `module load`):**

| tool | version | source |
|---|---|---|
| fastp | 1.0.1 | `EQ_WGS/envs/fastp.yml` |
| MEGAHIT | 1.2.9 | `EQ_WGS/envs/megahit.yml` |
| seqkit | 2.8.2 | `EQ_WGS/envs/megahit.yml` |
| CoverM | 0.7.0 | conda `eco_probes` |
| MMseqs2 | 18.8cc5c (conda `mmseqs2`); 16.747c6 in `eggnog` env | conda |
| DIAMOND | 2.0.15 | conda `eggnog` |
| eggNOG-mapper | 2.1.13 | conda `eggnog` |
| Bakta | 1.12.0 | conda `bakta` |
| pyrodigal | 3.7.1 | conda `pyrodigal` |
| CheckM2 | 1.1.0 | conda `checkm2` |
| dRep | 3.6.2 | conda `drep` |
| DAS_Tool | 1.1.7 | conda `das_tool` + run log |
| SemiBin | 2.2.1 | conda `semibin2` |
| GTDB-Tk | 2.7.2 (current) / 2.4.1 (superseded) | run logs |
| FoldSeek (05c validation only) | 10-941cd33 | `module load` |
| FoldSeek (02/03a/03b/04) | **8dc75c74ad0eddab73cfd905963d13bf74dc012b** | `M/tools/foldseek/bin/foldseek version` |
| transformers / torch (ESMFold) | 4.40.2 / 2.2.2 | pip dist-info in conda env `esmfold` |
| antiSMASH | 8.0.4 | `module load` |
| AMRFinderPlus | 4.2.7 | conda `amrfinder` |
| BUSCO | 6.0.0 | conda `busco` |
| inStrain | 1.10.0 | conda `instrain` |

---

## 13. What is still UNVERIFIED

Honest gap list. **Do not transcribe these into Methods without checking first.**

| gap | needed by | next step |
|---|---|---|
| All M.1 sampling/sequencing facts (kit, platform, read layout, depth, accessions, XRF instrument, sampling depths, expedition dates) | **M.1** | **STILL OPEN.** Ask Alejandra / the RAK-PR data-paper repo. Not in this workbench. Written into M.1 as one explicit `\todo` listing exactly these. |
| **The paper mixes GTDB releases** — §1.1 is R232, but 8 scripts across §1.4/§2.3/§2.4 read r226, incl. the **genus/species** pathogen census | **§2.4 owner** | **HALF CLOSED (7b).** M.3 now discloses which release each label came from, and M.11 (ex-M.10) was corrected R232 → **r226**. What remains is a *decision, not a gap*: re-run the pathogen census against R232 and move §2.4's number, or retain r226 with the disclosure. §2.4 is shipped, so 7b did not decide it. |
| §1.3 control floor + final calibrated fraction | §1.3 prose | pending job **50603810** (session 06 will report). **Not an M.6 gap** — M.6 is written; the method is fully specified on disk independent of the result. |
| §3.1 rarefaction + depth results | §3.1 | **CLOSED** — jobs 50613400 and 50618838 landed; numbers registered in `CATALOGUE_NUMBERS.md` §3.1. |

**Closed 2026-08-17, after the first gap list was written** — contig renaming scheme + QUAST (§4),
antiSMASH `--minimal` mode, AMRFinderPlus DB `2026-05-15.1` with the `--plus`/no-`--organism`
finding, and the KO panel sizes (28 UV / 13 osmolyte), all now **V** in their stage tables.

**Closed 2026-08-17 by session 7b, verified on disk before transcription** (details in
`CORRECTIONS.md` 2026-08-17, second entry):
- **pyrodigal parameters** — `-p meta -f gff -j 4 --no-stop-codon`, table 11, no min length, closed-ends off.
- **ProstT5 invocation** — `foldseek createdb --prostt5-model prostt5-f16.gguf --gpu 1`, no batch/length flags.
- **ESMFold** — HF `transformers` 4.40.2 port of `esmfold_v1`, chunk 64/32, 1,022-aa truncation, no recycles.
- **Structural clustering** — `foldseek cluster --single-step-clustering --cluster-mode 2 -s 4`, **not MMseqs2**.
- **The two FoldSeek builds** — GPU `8dc75c74` for 02/03/04, module `10-941cd33` for 05c only.
- **539 vs 542** — reconciled as a funnel, not a discrepancy: **542** samples submitted to DAS_Tool
  (`05_refine/sample_list.txt`) → **539** produced refined bins, 3 produced none. The registry row was
  relabelled accordingly. The HQ+MQ-filtered per-run counts (540 ⊆ 543) count a different thing and are
  not quoted in the paper.
- **≥100 aa query gate** — lives in `03a_prep_query.sh`, not `01_prep.sh`; substance unchanged.

The remaining gaps are not artefacts of missing information: the eight extraction agents assigned to
this rewrite were killed mid-run by an account-level API limit, so everything here was verified
directly in the main session instead. What is left needs either a pending job, a decision, or
Alejandra.

## 2026-10-13 review additions (append-only)

The three-round dark-matter review work uses IBEX root `06_darkmatter/p4_review_20261013/` with `scripts/`, `results/`, and `logs/`; every compute job runs in Slurm account pi-hohndor/partition batch. Local portable compact driver: `bash reproducibility/darkmatter/run_current.sh`. Read its README for the distinction between compact reproduction, original remote analyses and incomplete raw-data provenance. The historical whole-paper verifier is preserved; current review assertions are separate and do not certify the old registry.

Completed jobs:53678755CDScovariates;53678768supersededapproximatePAF;53679350explicitbasealignment;53678944AntiFam;53680096coordinateverifiedoppositeoverlap;53680202AntiFamsummary/cultureloci;53679366selectionbootstrap/rawFoldseek;53679245siteaggregate;53679246environmentcapture. Failed attempts and their diagnoses are retained in review records. `upstream_code_manifest.tsv` identifies162preserved original clean-pipeline source files. `review_inputs.sha256.tsv` pins compact evidence; rawPDBs/FASTAs/PAFsremainatdocumentedIBEXpaths.

Final follow-up53680886completed the all-DNA-supported culture locus audit and cluster90-to-admitted-family mapping. Large-file hash job53680819completed1,494unique paths(292.8GB),zeroerrors; `ibex_files_complete.tsv.gz`fillsall1,497large-file references in the clean-analysis inventory. Original inventories and failed-attempt logs are retained. No raw-read file was changed.

## 2026-10-13 round-4 audit additions

Gemini4ArgonHigh via agy reviewed the frozen manuscript and code. Root adjudicated
findings and ran IBEX job53695688 (40seconds). Corrected soil-position and partial-
correlation outputs, the rRNA parser audit, inputs and source hashes are under
`reproducibility/darkmatter/round4_corrections/`; the job script is `scripts/rerun.py`.
No historical input or source was overwritten. `round4_checks.py` adds40read-only
checks to the49prior revision checks. The full compact driver passes and compiles
28pages. See the round4 response for reviewer inaccuracies and unresolved provenance.
