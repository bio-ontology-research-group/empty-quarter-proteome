# IBEX working-directory inventory — Empty Quarter MAGs

**Path:** `M=/ibex/scratch/projects/c2014/rob/metagenomics`
**Project:** Empty Quarter MAGs paper (genome/protein/metabolism companion to the amplicon + phages papers)
**Date:** 2026-06-30
**Author of this doc:** handover inventory for Alejandra (velazqal). Read-only survey (`ls`/`du`/`head`).
**Total size:** ~14.6 TB across 38 top-level subdirectories.

> **Historical inventory only.** Sizes, job states, database releases and labels below are a point-in-time scratch
> survey from 2026-06-30, not current manuscript evidence. Current catalogue definitions and promoted claims live
> in `binning-and-annotation-scope.md`, `reproducibility/claims.tsv` and `manuscript/main.tex`.

This documents what each IBEX scratch subdirectory contained at the survey date, its approximate size and its
then-observed pipeline state. A copy was placed at `$M/INVENTORY.md`. The active paper deliberately labels more
than one data universe; this inventory does not decide which results are current or promoted.

> **Caveat on "scratch":** this is IBEX `/scratch`, subject to purge policy. Canonical, kept data live under
> `/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/`. Heavy GPU/structural work (ProstT5/FoldSeek for the
> 753M assembly catalog) lives on ORIX `/mnt/data/u/hohndor/eq-prostt5/`. Treat anything here as reproducible
> intermediate unless promoted to `EmptyQuarter_Data` or the paper repo `results/`.

## Three historical workbench pillars
1. **Functional & structural dark matter** — 51.9M MAG proteins → 23.4M reps → 2.84M UniRef50 no-hit →
   reference-unassigned candidate models and recurrent candidate families (historical `novelfold` naming).
   (`04_structure`, `06_darkmatter`, `06_esmfold_novelfold`)
2. **Novel/extreme metabolism ("life on air")** — atmospheric H₂/CO oxidation, V-nitrogenase, desiccation/UV
   survival, rhodopsins, BGCs. (`05_metabolism`, `06_atmospheric`, `06_metabolism`, `06_aquatic_scavenging`, `06_smetana`)
3. **Reference bias / PICRUSt2 failure** — MAG-encoded vs predicted gene content. (`05_picrust_validation`)

Plus reviewer-proofing: contamination/pathogen controls (`06_resistome_safety`), global comparison (`databases`),
abundance (CoverM). **Phages + mobilome (`06_phage`, `06_mobilome`) belong to the SISTER phages repo** but the
compute physically lives here — they are the two biggest space consumers.

---

## Summary table

| Dir | Size | Purpose | Pillar / RQ | Status |
|-----|-----:|---------|-------------|--------|
| `01_qc` | 0 | empty placeholder (QC done in canonical tree) | infra | empty/stale |
| `02_taxonomy` | 14G | GTDB-Tk r226 — single_sample MAGs | classification | done |
| `02_taxonomy_self_supervised` | 5.3G | GTDB-Tk r226 — self_supervised MAGs | classification | done |
| `03_annotation` | 28M | driver MAG lists + array cursors (Bakta/EggNOG) | annotation driver | done |
| `04_drep_self_vs_single` | 407G | dRep 95% ANI on raw union (37,085 MAGs) → 9,211 reps; 3,386 novel self-sup | dedup/scope | done |
| `04_drep_merged` | — | fresh merged dRep of 37,949 raw bins (velazqal, job 47940799) | dedup | running |
| `04_functional` | 20G | per-sample functional/refine intermediates (Dr/PR/Sr) | annotation | done |
| `04_protein_catalog` | 2.5T | **753M-protein assembly-level catalog** (all contigs, both binners) | dark-matter (everything arm) | done (not in paper numbers) |
| `04_structure` | 472G | **51.9M MAG proteins → reps → FoldSeek/ProstT5 structural novelty** | Pillar 1 (core) | done |
| `05_metabolism` | 89G | gapseq GEMs + community design + reference soils | Pillar 2 | done/partial |
| `05_picrust_validation` | 42M | MAG-encoded KO vs PICRUSt2 prediction | **RQ-C1** / Pillar 3 | done |
| `05_refine` | 392G | DAS_Tool per-sample refinement (14,426 refined bins) + `_drep/` → **5,438** dereplicated cross-binner MAGs | binning/dedup | done |
| `06_aquatic_scavenging` | 3.6G | gapsmith/FBA nutrient-scavenging models | Pillar 2 | done |
| `06_atmospheric` | 81M | atmospheric guild (hhyL/coxL/RuBisCO) census + integrate | **RQ-D1/D2/D3** | done |
| `06_darkmatter` | 182G | gapsmith dark-protein function assignment (C·V·T, LRO-cal) | Pillar 1 (core) | done |
| `06_esmfold_novelfold` | 11G | historical ESMFold + FoldSeek re-screen of reference-unassigned candidate models vs AFDB/PDB | **RQ-A3/B1** | done |
| `06_metabolism` | 259G | marker census + antiSMASH BGC + **bgc_novelty (BiG-SLiCE)** + rhodopsin + UV | Pillar 2 / **RQ-F2** | done |
| `06_mobilome` | 4.9T | AMRFinder/defense/plasmid/CRISPR/inStrain — **sister phages repo** | (phages paper) | done (not this paper) |
| `06_novelty` | 0 | empty placeholder | — | empty/stale |
| `06_pangenome` | 13G | species pangenome openness / gene-family enrichment | **RQ-C2** (job 47843640) | done |
| `06_phage` | 3.5T | geNomad/VIBRANT/PhaBOX/vConTACT2/inStrain — **sister phages repo** | (phages paper) | done (not this paper) |
| `06_resistome_safety` | 347M | intrinsic resistome, metal STRESS, pathogen + VFDB census | reviewer-proofing (AMR/safety) | done |
| `06_smetana` | 11G | community metabolic modelling / cross-feeding (16-sample pilot) | **RQ-G1** | done (pilot) |
| `07_aquatic_paper` | 89G | analyses for a separate aquatic/enrichment paper | (other paper) | active |
| `07_nitrogen` | 0 | empty placeholder | Pillar 2 (planned) | empty/stale |
| `08_dormancy` | 0 | empty placeholder | Pillar 2 (planned) | empty/stale |
| `09_kraken2` | 0 | empty placeholder | — | empty/stale |
| `databases` | 1.8T | **reference DBs** (iphop, spire, mgnify90, dram, deserts, gem, imgvr…) | infra / global comparison | reference (not output) |
| `dcgm` | 2.7M | DCGM GPU-monitor logs from ProstT5/FoldSeek/ESMFold jobs | infra/logs | stale |
| `envs` | 4K | one conda env yml (`metabolism-env.yml`) | infra | minimal |
| `family6_validation` | 3.0M | deep validation of one reference-unassigned candidate family (historical “novel fold” label; ESMFold/hhsearch/signalp…) | Pillar 1 (case study) | done |
| `handover` | 424M | handover material — `literature/` (cited PDFs) | handover | partial |
| `logs` | 8.6G | Slurm logs across all jobs | infra | active |
| `metadata` | 6.1M | sample/MAG lists, CheckM2 tables, dRep reps, crosswalks | infra (core) | done |
| `results` | 25M | CheckM2 full tables (single + self_supervised) | infra | done |
| `results_agg` | 29M | **curated cross-MAG result tables** (mags/ + phages/) | promotion staging | done |
| `scripts` | 93M | **188 pipeline scripts** (the reproducibility record) + `ibex/`, `analysis/` | infra (core) | active |
| `scripts_tmp` | 76K | ~16 scratch scripts (gtdbtk self-sup, foldseek microchunk, DAS_Tool) | infra | superseded/scratch |
| `tools` | 3.8G | installed binaries: foldseek, gapsmith, p2rank, cobra-venv | infra | active |

---

## Per-directory detail

### Pipeline core (binning → annotation → taxonomy)

**`02_taxonomy` (14G)** — historical GTDB-Tk r226 classification of the **single_sample** MAGs.
`gtdbtk_input/`, `gtdbtk_input_chunks/`, `gtdbtk_out/`, `gtdbtk_out_chunks/`. It remains provenance
for the historical functional snapshot but is superseded for manuscript taxonomy by the canonical
GTDB-Tk 2.7.2/R232 ANI95 census: 5,438 HQ/MQ representatives, 51 named phyla, no domain/phylum/class
gap, 11 first order gaps (10 adjudicated with zero proposals; one low-MSA candidate in recovery). **Done; historical.**

**`02_taxonomy_self_supervised` (5.3G)** — Same, for the **self_supervised** MAGs (`gtdbtk_input_chunks/`,
`gtdbtk_out_chunks/`). Confirms self_supervised is fully classified. **Done.**

**`03_annotation` (28M)** — Not the annotations themselves (those are distributed in the canonical tree per
`binning-and-annotation-scope.md`); this holds the **driver lists + array cursors**:
`mag_list.txt` (2.6M, single_sample), `mag_list_self_supervised.txt` (3.8M), `…_all.txt` (15M),
`…_new.txt` (598K = 3,386 novel self-sup), `mag_list_single_sample_lq.txt` (5.4M LQ), `state/`. **Done.**

**`04_drep_self_vs_single` (407G)** — dRep at 95% ANI on the **union** of both binning strategies (37,085 HQ+MQ),
raw bins (no DAS_Tool step) → **9,211** reps (`out/dereplicated_genomes/`). `genomeInfo.csv`, `new_mags.tsv`,
`redundant_mags.tsv` (the 3,386 non-redundant self_supervised flag), `inputs/`, `out/`. Drives the scope statement. **Done.**

**`04_drep_merged`** — a **fresh merged dRep** of the 37,949 raw single+self bins, created **2026-07-02 by velazqal
(Alejandra)** — job **47940799** running as of writing (`out/dereplicated_genomes/` not yet populated). Distinct
from `04_drep_self_vs_single` and from the refined dRep below.

**`04_functional` (20G)** — Per-sample sub-dirs named by sample (`10Dr2`, `11PRr3`, … D=Dune, PR=interdune,
S=other compartments). Functional-annotation intermediates (Jun 2026). **Done.**

**`05_refine` (392G)** — Per-sample **DAS_Tool** bin refinement (`<sample>/<sample>_DASTool_bins/`, aggregating
single+self SemiBin2 → best non-redundant set per sample) → **14,426 refined bins** total (also in
`04_protein_catalog/contig_to_bin_refined.tsv.gz`). **`05_refine/_drep/`** dereplicates those refined bins (95% ANI,
12,433 HQ+MQ in the frozen historical membership table) → **5,438 representative MAGs** in
`_drep/out/dereplicated_genomes/` (named `<sample>__{single,self}_SemiBin_N.fa`). The current refined HQ+MQ total
is 12,434, so the one-bin membership non-closure must remain explicit; see `binning-and-annotation-scope.md`.
Scripts `50–54_*` in `scripts_tmp/`. **Done at survey time.**

### Pillar 1 — functional & structural dark matter

**`04_structure` (472G) — CORE.** The paper's structural dark-matter ladder.
- `pooled_proteins.faa` (18G) = the **51.9M single_sample MAG proteins** (the headline catalog).
- `cluster90_rep_seq.fasta` (7.7G), `cluster90_cluster.tsv` (2.8G), `cluster90_all_seqs.fasta` (19G) = mmseqs
  linclust90 reps.
- `foldseek_results/` (607M): `struct_nohit.fasta` (451M, reference-relative AFDB no-hit queries),
  `struct_nohit_ids.txt`, `afdb_hits.tsv` (→ symlink to `prostt5_foldseek_orix/merged_afdb_hits.tsv`),
  `c4_characterised_summary.tsv` (80.2M hit rows; only 9.25% of rescued queries hit a SwissProt target).
- `tier1_foldseek/`, `tier2_foldseek/`, `prostt5_foldseek_orix/`, `afdb_pad/`, `novelty_analysis/`,
  `fsdb_pdb*` (FoldSeek PDB ref DB). Many `tmp_foldseek_*` scratch dirs (purgeable). **Done.**

**`06_darkmatter` (182G) — CORE.** gapsmith dark-protein **function assignment** engine output (C·V·T scoring,
LRO-calibrated). Key: `scores_v3/candidate_scores_v3.tsv` (11G) + `discovery_top.tsv`, `pathway_assoc_*`,
`reliable_byOG.tsv`, `gap_equations.tsv`; plus `calibration*/`, `dossier/` (cand73 + json case files),
`esmfold/`, `esmfold_v3/`, `pfba/`, `reallydark/`, `orthogroups/`, `localization/`, `antismash/`,
`darkfam_ecology/`, `scripts/`. Versioned `scores`/`scores_v2`/`scores_v3` (v3 current). **Done.**

**`06_esmfold_novelfold` (11G)** — Historical ESMFold predictions + FoldSeek re-screen of recurrent,
reference-unassigned candidate representatives against AFDB and PDB. `best_afdb_tm.tsv`, `best_pdb_tm.tsv` (best TM-score per query),
`128_novelfold_foldseek_reconfirm.sh`, sharded `hits_afdb.m8.*`, `fs_tmp/`. Supports **RQ-A3/B1**. **Done.**

**`family6_validation` (3.0M)** — Deep single-family validation of one reference-unassigned candidate
("family 6"; historical “novel fold” label):
`esmfold/`, `foldseek/`, `hhsearch/`, `hmm/`, `hmm_screen/`, `msa/`, `multimer/`, `signalp/`, `synteny/`,
`seqs/`, `scripts/`. A worked case study for testing recurrence and structural evidence, not proof of a new fold.
**Done.**

**`04_protein_catalog` (2.5T)** — The **assembly-level "annotate everything" catalog** (753,305,533 proteins
across 552 samples, binned + unbinned, both strategies). NOT in the paper's headline numbers.
`assembly_proteins_pooled.faa.gz` (94G), `assembly_proteins.mmseqsdb` (130G) + `_h`/`.index`/`.lookup` (26G),
`existing_reps.mmseqsdb` (6.5G), `contig_to_bin_refined.tsv.gz` (34M, links both binners),
`delta_linclust90/`, `delta_vs_uniref50/`, `assembly_vs_existing_reps/`, `orix_*` (ProstT5 hand-off to ORIX).
`_pooled_record_count.txt` is the (empty/in-progress) record-count verification. **Done** (structural arm on ORIX).

### Pillar 2 — novel / extreme metabolism

**`05_metabolism` (89G)** — gapseq genome-scale models + community design. `gapseq_out/`, `communities/`
(`communities_index.tsv` + per-community `T1__siteNNN__{Dr,PR,Sr}/`), `reference_soils/` (comparator genomes for
global benchmarking), `hq_mag_list.tsv`, `coverage_samples.tsv`, `state/`, `tmp/`. **Done/partial** (superseded in
part by `06_smetana` cross-feeding).

**`06_atmospheric` (81M)** — "Life on air" guild census. Sub-dirs: `classify/`, `queries/`, `integrate/`
(`atmo_guild_per_mag.tsv`, `atmo_guild_by_phylum.tsv`, `atmo_operon_context.tsv`, `INTEGRATE_DONE`),
`carbon_fixation/`, `autotrophy_fba/`, `rubisco/` (`rbcL_form_per_protein.tsv`, alignment, clustering — the
CBB/RuBisCO form validation: known Calvin cycle, ~16% form-tightened), `operon_dark/`, `dark_recover_gas/`,
`cultures/`, `refs/`. Supports **RQ-D1/D2/D3**. **Done.**

**`06_metabolism` (259G)** — Broad metabolism layer.
- `marker_summary.tsv` / `marker_hits.tsv` — KO marker census (recA, otsA/B trehalose, coxSML CODH, hydB NiFe,
  rbcL RuBisCO, ectABCD ectoine, phrB photolyase, nifHDK — note nif very rare: 16–20 MAGs).
- `antismash/` — BGC predictions.
- `bgc_novelty/` (73M) — **BiG-SLiCE** clustering of BGCs vs MIBiG 3.1 (`bigslice_input/`, `bigslice_out/`,
  `mibig_gbk_3.1/`, `results/`). Supports **RQ-F2** (job 47874840).
- `rhodopsin/` — microbial rhodopsin census.
- `uv_survival/` (15M) — desiccation/UV markers: `uv_marker_summary.tsv`, `uv_per_mag.tsv`, `UV.md`,
  `phylum_uv_markers.tsv`. **Done.**

**`06_aquatic_scavenging` (3.6G)** — gapsmith/FBA nutrient-scavenging models. `results/` has finished reports
(`fba_scavenging_report.txt`, `scavenging_matrix.tsv`, `dem123.tsv`, `.analysis_done` flags), `gapsmith_models/`,
`scripts/`, `metadata/`, `_smoketest/`. **Done.**

**`06_smetana` (11G)** — Community metabolic modelling / cross-feeding pilot (16 communities). Has a good
`README.md`. `models_sbml/` (316 gapseq GEMs), `smetana_out/global/` (MIP/MRO), `fba_crossfeed/`
(`mag_exchange.tsv`, `edges.tsv`, `community_summary.tsv`), pilot design TSVs, cplex bench scripts (no MILP solver
on cluster → honest negative: guild-as-feeder REFUTED). Supports **RQ-G1**. **Done (pilot only).**

### Pillar 3 — reference bias

**`05_picrust_validation` (42M)** — MAG-encoded gene content vs PICRUSt2-from-16S prediction.
`mag_annotation_stats.tsv` (per-MAG protein/KO counts + CoverM relab), `mag_ko_count_by_sample.tsv`,
`mag_ko_presence_by_sample.tsv`, `mag_ko_relab_by_sample.tsv` (14M each). Supports **RQ-C1**. **Done.**

**`06_pangenome` (13G)** — Species pangenome openness + gene-family enrichment. `out/family_table.tsv` (212M),
`species_openness.tsv`, `species_members.tsv`, `enrichment_summary.json`; `work/`. Supports **RQ-C2** (job
47843640). **Done.**

### Reviewer-proofing — safety / contamination

**`06_resistome_safety` (347M)** — Intrinsic resistome + biosafety. Excellent `README.md`. Sub-dirs:
`amr_metal/` (STRESS metal resistance + XRF link, `metal_census.tsv`, `metal_xrf_corr.tsv`),
`pathogens/` (contamination control: only 7/15,469 MAGs are named pathogens — `pathogen_genus_census.tsv`,
`PATHOGENS.md`), `vfdb/` (virulence census, signal collapses with stringency; `vfdb_hits.tsv` ~346M),
`code/` (`phase_amr_metal.py`, `pathogen_census.py`, `vfdb_analyze.py`). Curated to paper repo
`results/{amr_metal,pathogens,vfdb}/`. **Done.**

### Sister phages repo (physically here, not this paper)

**`06_mobilome` (4.9T) — BIGGEST consumer.** AMRFinder, DefenseFinder, PADLOC, CCTyper, plasmidome,
inStrain, iPHoP, AMG networks, dN/dS, defense islands. `mobilome_refs.fna` (2.8G), `mag_index.tsv`,
`abundance*/`, `defense*/`, `plasmidome/`, `pangenome_flux/`, `flagship_struct/`. **Belongs to the phages paper.**

**`06_phage` (3.5T)** — geNomad, VIBRANT, VirSorter2, PhaBOX, CheckV, vConTACT2, inStrain, iPHoP, vOTU
catalog, BORG screen, global IMG/VR compare, DRAMv, terL tree. **Belongs to the phages paper.**

> AMR (intrinsic, host-genome) lives in `06_resistome_safety` (this paper). The mobile/defense/phage census is
> the sister repo's; aggregated tables for both are mirrored in `results_agg/` (below).

### Other / separate work

**`07_aquatic_paper` (89G)** — Analyses for a **separate aquatic/enrichment-culture paper** (numbered stages
`00_crosswalk` … `13_mobilome_contrast`: salinity, plasmid HGT, strain pangenome, endemicity, metabolic
complementarity, dark-matter context, phylo-trait). Active, **not the MAGs paper.**

### Infrastructure / reference

**`scripts` (93M) — the reproducibility record.** 188 numbered pipeline scripts (`00_setup_drep_env.sh` …
`130_mag_ko*`), grouped roughly: 00–09 binning/CheckM2/Bakta/EggNOG/GTDB-Tk; 05/06 pooling/clustering/UniRef50/
ProstT5/FoldSeek; 10–16 gapseq/geNomad/VirSorter2/VIBRANT/PhaBOX/CheckV; 20/30 CoverM/comparator; 40–62
tier1/tier2 fold clustering; 70–80 gapsmith/orthogroups/scoring/LRO-calibration/localization; 81–98
antiSMASH/marker-census/InterProScan/mobilome/defense/AMR/ESMFold; 99–131 size-diagnostics/rescore/iphop/
plasmid/global-specificity. Sub-dirs: `ibex/` (newer numbered jobs 120–131: foldseek, antidefense, UV survival,
defense-on-phage/plasmid, desert/global specificity), `analysis/` (`aggregate_mobilome.py` etc.), `gspa/`,
`eco_probes/`, `dcgm/`, `pilot_af2/`. Plus loose helpers (`build_metadata.py`, `build_communities.py`,
`lro_calibration.py`, `score_candidates.py`, `community_design.md`). **The canonical home is `workbench/scripts/`;
commit there.**

**`scripts_tmp` (76K)** — ~16 scratch scripts (self-sup GTDB-Tk prep, foldseek microchunking, tier1 aggregation,
DAS_Tool refine 50–54). Superseded staging; safe to treat as historical.

**`metadata` (6.1M)** — Core sample/MAG bookkeeping: `sample_metadata.tsv` (287K), `canonical_samples.txt`,
`checkm2_all_samples.tsv` (5.5M), `drep_reps.txt`/`drep_reps_sub.txt`, `self_supervised_samples*`,
`assemblies_*` audit set, `alejandra_sample_counts.tsv`. **Done.**

**`results` (25M)** — `checkm2_single_sample_full.tsv` (8.7M), `checkm2_self_supervised_full.tsv` (16M). **Done.**

**`results_agg` (29M)** — **Curated cross-MAG tables**, the promotion staging set.
`mags/`: `amr_class_summary.tsv`, `amr_long.tsv`, `bgc_per_mag.tsv`, `bgc_product_summary.tsv`,
`metabolism_marker_long.tsv` + per-phylum `coxL/hydB/nifH/rbcL` summaries.
`phages/`: `defense_per_mag.tsv`, `defense_systems_long.tsv`, `retrons.tsv`, `votu_catalog.tsv`,
`borg_*`, `votu_*_summary.tsv`. **Done.**

**`databases` (1.8T) — reference DBs, NOT project output.** `iphop_db` (572G), `spire` (273G), `mgnify90`
(221G), `dram` (188G), `deserts` (178G, other-desert genomes for EQ-specificity), `gem` (157G), `imgvr`
(154G), `gov2` (16K), `checkv` (6.4G), `checkm2` (2.9G), `mob_suite` (2.9G), `phabox_db_v2_2` (1.6G),
`busco_downloads` (49M), `cluster_one-1.0.jar`. The `deserts`/`gem`/`spire`/`imgvr`/`mgnify90` sets are the
**global-comparison** reviewer-proofing DBs. **Reference (re-downloadable).**

**`tools` (3.8G)** — Installed binaries: `foldseek/` (+ avx2/gpu tarballs 460M), `gapsmith/` & `gapsmith_src/`
(the dark-matter function engine), `p2rank/`, `cobra-venv/` (FBA python env). **Active.**

**`envs` (4K)** — Single conda spec `metabolism-env.yml`. (Most envs are managed elsewhere / in workbench.)

**`logs` (8.6G)** — Slurm stdout/stderr across all jobs. **Active.**

**`dcgm` (2.7M)** — DCGM GPU-monitor logs from ProstT5/FoldSeek/ESMFold GPU jobs (`hohndor_esmfold_v3`,
`hohndor_foldseek`, `hohndor_fs_*`, `hohndor_prostt5`). Diagnostic; **stale.**

**`handover` (424M)** — Handover staging. Currently only `literature/` (cited PDFs + `fetch_papers.py` +
`desert_arid_cold/` subset — genomad, atacama-virome, defense-islands, permafrost, etc.). **Partial** — this
`INVENTORY.md` extends it.

**Empty placeholders / stale (0 bytes):** `01_qc`, `06_novelty`, `07_nitrogen`, `08_dormancy`, `09_kraken2`.
Scaffolding for planned/relocated analyses (QC done in canonical tree; nitrogen/dormancy folded into
`06_metabolism`/`06_atmospheric`; kraken2 not pursued). Safe to ignore.

---

## Historical handover notes
- **Scope at survey time:** the 51.9 M-protein ladder came from the 15,469 single-sample-MAG snapshot, whereas
  current taxonomy and recurrent-family analyses use the frozen 5,438 ANI95 representatives. The 753 M
  assembly-level catalogue is a separate discovery universe. See `binning-and-annotation-scope.md` for the active
  separation rules.
- **Space:** 4 dirs hold ~11.7 TB of the ~14.6 TB: `06_mobilome` (4.9T) + `06_phage` (3.5T) are the **sister
  phages repo** (not this paper); `04_protein_catalog` (2.5T) is the everything-arm; `databases` (1.8T) is
  re-downloadable reference DBs. The MAGs paper's own core (`04_structure` 472G, `06_darkmatter` 182G,
  `06_metabolism` 259G) is comparatively small.
- **Versioning:** prefer `…_v3` / latest (`06_darkmatter/scores_v3`, `esmfold_v3`); earlier `scores`/`scores_v2`
  are superseded. `scripts_tmp`, `dcgm`, many `04_structure/tmp_foldseek_*` are scratch/purgeable.
- **Then-curated outputs:** `results_agg/mags/`, `06_resistome_safety/`, `06_smetana/README.md`,
  `06_metabolism/{marker_summary,uv_survival}`, `06_atmospheric/integrate/`, `06_pangenome/out/`. Promotion now
  requires a checksum-closed source/claim entry; scratch location alone is not sufficient.
