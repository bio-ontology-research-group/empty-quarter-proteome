# Dossier 01 — Catalogue and taxonomy (§1.1 · M.2 · M.3 · tab:scale)

`M = /ibex/scratch/projects/c2014/rob/metagenomics` · `R = notes/CATALOGUE_NUMBERS.md` · **V** = I opened the file this session; "guide-**V**" = verified in `notes/PIPELINE_GUIDE.md`, not re-checked here. Versions **V** in `/ibex/scratch/projects/c2014/rob/miniconda3/envs/*/conda-meta/` (`checkm2-1.1.0`, `semibin-2.2.1`, `das_tool-1.1.7`, `drep-3.6.2`) and `/ibex/project/c2014/alelopezv/EQ_WGS/envs/{fastp,megahit,quast}.yml` (`fastp=1.0.1`, `megahit=1.2.9` + `seqkit=2.8.2`, `quast=5.2.0`); dRep, DAS_Tool and both GTDB-Tk runs also stamp their own logs.

## Question
Does an extreme, hyperarid biome never before sampled at genome resolution at this scale harbour deep-branching lineages absent from reference taxonomy? Answer: no — every species representative is placed at domain, phylum and class; novelty starts at order and is concentrated at genus/species. That relocates the paper's novelty claim from taxonomy to gene content, structure and metabolism.

## Basis — the denominator chain, with the exact rule at each step
| N | what it is | exact rule | file | V |
|---|---|---|---|---|
| **551** | canonical soil samples (headline) | rows of the metadata table = soil samples with complete metadata, 5 expeditions | `M/metadata/sample_metadata.tsv` (551 rows + header) | R:25 |
| **549** | assemblies | 551 − 1 failed − 1 empty. **The two are now named (V, this session):** `51Sr1` = `has_contigs=False`, no contigs file at all (failed); `F5Dr1` = contigs file present but **0 bytes** (empty). All 549 others are non-empty | `M/metadata/sample_metadata.tsv` cols `contigs`/`has_contigs` (550 True / 1 False) + `stat` on each path | R:26 |
| **542** | DAS_Tool input samples | samples binned by **both** SemiBin2 runs (intersection) | `M/05_refine/sample_list.txt`, `wc -l`=542 | **V** |
| **539** | samples yielding refined bins | ≥1 DAS_Tool bin produced; 3 inputs yielded none (a funnel, not a discrepancy) | `find M/05_refine -maxdepth 2 -name '*_DASTool_contig2bin.tsv'`=539 | **V** |
| **14,426** | refined bins | all DAS_Tool bins, one row each | `M/05_refine/_checkm2_refined/quality_report.tsv`, 14,426 data rows | **V** |
| 4,225 / 8,209 / 1,992 | HQ / MQ / LQ | HQ = comp ≥90 **and** contam ≤5; MQ = comp ≥50 and contam ≤10; else LQ | same file; my awk on cols 2,3 reproduces all three exactly | **V** |
| **12,434** | HQ+MQ | HQ + MQ, contamination ceiling **inclusive** (≤10) | same file | R:60 |
| **12,433** | dRep clustered | 12,434 minus one bin dRep dropped internally: `25PRr2__self_SemiBin_2_sub` | `_drep/out/data_tables/Cdb.csv`; `reproducibility/drep_membership_reconciliation.tsv` | **V** |
| **5,438** | species representatives | one per ANI95 secondary cluster (2,259 HQ + 3,179 MQ) | `_drep/out/data_tables/Wdb.csv`; `dereplicated_genomes/*.fa` | R:61 |
| **37,949** | HQ+MQ union of the two *unrefined* binnings | 16,333 pretrained (4,030 HQ + 12,303 MQ, of 48,424 bins) + 21,616 self-supervised | `M/results/checkm2_single_sample_full.tsv` (awk **V**) · `M/03_annotation/mag_list_self_supervised.txt` (`wc -l`=21,616 **V**) | R:272 |

**Why 37,949 exists separately from 5,438.** It is not a superseded catalogue. It is the *pre-refinement* union, kept for one live purpose only: §1.2's genome-resolved protein census maps assembly proteins onto contigs in those 37,949 bins (`contig_to_bin.tsv.gz` covers 100 % of them). R:275 — *"This is the only live use of the union; taxonomy/species counts come from the refined 5,438 catalogue (§1.1), not from here."* Never quote it as a MAG count.

**Sites: 60, not 61.** Metadata carries **61 distinct numeric labels** = 59 numbered transect sites + 2 trip-1-only aliases of *named* locations (62 and 64, 4 samples). Both companion 16S papers say 60 primary sites + 10 named locations. A change to 61 was made and reverted (CORRECTIONS 2026-08-17). §3.1/M.12 legitimately count **61 sampled locations** in occupancy statistics — different quantity, different word.

## How it ran
| step | tool · version | parameters as run | script | logs confirm? |
|---|---|---|---|---|
| read QC | fastp **1.0.1** | `--correction --qualified_quality_phred 20 --unqualified_percent_limit 40 --length_required 50 --complexity_threshold 20`; paired **and** orphan outputs kept | `EQ_WGS/rules/assembly/…` rule `qc_reads` | guide-**V** |
| assembly | MEGAHIT **1.2.9** | `--presets meta-large --min-contig-len 500 --min-count 2 --prune-level 1 --low-local-ratio 0.3`, `-1/-2/-r` | `EQ_WGS` rule `assembly` + `config/config.yml` | guide-**V** |
| contig rename | — | `sed` rewrites headers to `<sample>_<N>`; MEGAHIT `k141_` prefix and `flag=/multi=/len=` discarded; `final.contigs.fa` deleted | `rules/assembly/rename_contigs.rule` | guide-**V** |
| assembly QC | QUAST **5.2.0** | all defaults | `rules/assembly/qc_assembly.rule` | guide-**V** |
| binning A (**pretrained**) | SemiBin2 **2.2.1** | `minimap2 -a -x sr -t 16 $CONTIGS $R1 $R2 \| samtools sort` → `SemiBin2 single_easy_bin -i -b -o --environment global --sequencing-type short_read -t 16` — **`-m` never passed** | `M/scripts/01_semibin2_single.sh` (PI re-bin, outside the Snakemake) | **V** verbatim |
| binning B (**self-supervised**) | SemiBin2 **2.2.1** | BWA-MEM (paired **+** orphan) → `SemiBin2 single_easy_bin --tmpdir … -m 1000 --input-fasta --input-bam --self-supervised --threads 24` | `EQ_WGS/rules/binning/semibin.rule`; `min_len: 1000` at `config.yml:128` | **V** verbatim |
| refinement | DAS_Tool **1.1.7** | `--bins single.contigs2bin.tsv,self.contigs2bin.tsv --labels single,self --contigs <sample>_contigs.fa --search_engine diamond --proteins NULL --write_bins TRUE --threads 16 --score_threshold 0.5 --duplicate_penalty 0.6 --megabin_penalty 0.5 --max_iter_post_threshold 10` | per-sample, one run per sample | **V** `M/05_refine/22PRr1/22PRr1_DASTool.log`, 2026-06-04; params identical in 10Dr2, F41PRr1, V51PRr2. `score_threshold`, `duplicate_penalty`, `megabin_penalty` are all DAS_Tool **defaults** — nothing was tuned |
| MAG QC | CheckM2 **1.1.0** | DB `uniref100.KO.1.dmnd` (3,082,500,605 B, dated 2025-02-08) | — | guide-**V**; version **V** in conda-meta |
| dereplication | dRep **3.6.2** `dereplicate` | `S_ani 0.95` · `P_ani 0.9` · `S_algorithm fastANI` · `completeness 50.0` · `contamination 10.0` · `length 50000` · `cov_thresh 0.1` · `coverage_method larger` · `clusterAlg average` · `MASH_sketch 1000` · weights comp 1 / contam **5** / N50 0.5 / centrality 1 / strain-het 1 / size 0 · `genomeInfo` supplied from CheckM2 | — | **V** `_drep/out/log/cluster_arguments.json`, run 2026-06-07 |
| taxonomy (**current, R232**) | GTDB-Tk **2.7.2**, ref **r232** | `classify_wf --genome_dir … --extension fa --cpus 64 --pplacer_cpus 32 --keep_intermediates --write_single_copy_genes`; ANI pre-screen **on** (skani v0.3.2); internal gene calls Prodigal V2.6.3 | — | **V** `M/rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/gtdbtk_2.7.2_r232/run/gtdbtk.log`, 2026-07-31, job 49670976. Run covers **5,742** genomes (201 ANI-screened + 5,541 marker-identified) = our 5,438 + 304 LQ salvage |
| taxonomy (**superseded, r226**) | GTDB-Tk **2.4.1**, ref **r226** | `classify_wf … --cpus 32 --pplacer_cpus 16 --skip_ani_screen` | — | **V** `M/05_refine/_gtdbtk_refined/gtdbtk.log`, 2026-06-04 |

**Where the paper's taxonomy table comes from.** Not `_gtdbtk_refined/` (that is r226). The census is derived by `M/rebuild_20260729/scripts/derive_r232_manuscript_census.py`, driven by `129_derive_r232_manuscript_census.sbatch --expected-legacy 5438 --expected-lq 304`, checksum-gated on `LEGACY_TAXONOMIC_SCREEN_SUMMARY_COMPLETE` + `taxonomic_screen_summary.sha256`, fails closed. Downstream consumers read `M/05_refine/_gtdbtk_r232/` — a verbatim subset of the same R232 run to our reps (**V**: 4,847 bac120 + 591 ar53 = 5,438 exactly).

**First-empty-rank census — the method.** Walk `RANKS = (domain, phylum, class, order, family, genus, species)` per representative; assign it to the **first rank carrying no named assignment**; report `(order, family, genus, species, fully_assigned)`. Classes are **mutually exclusive** — a genome counted at "genus" is not also counted at "species", so the counts do **not** accumulate and "11 novel orders" would misread it.

**The 10 adjudicated order gaps (V, read row by row).** Of the 11 empty-order reps, 10 passed the staging gate; the 11th, `48Sr2__single_SemiBin_167_sub` (MQ, Gammaproteobacteria, RED 0.75507), was excluded on the recorded rule **`msa_percent_below_50`** (46.49 %) and routed to 48_S coassembly recovery. The 10 went through GTDB-Tk 2.7.2 `de_novo_wf` trees + per-phylum R232 RED; 10/10 present in trees.
| disposition | n | detail |
|---|---|---|
| `unresolved_single_genome_order_gap` | 6 | `no_multi_genome_candidate_only_clade_for_rank_test` (22PRr1, 28Dr3, 54PRr2, F41Dr3, S22PRr1, V51PRr2) |
| `unresolved_candidate_clade_red_incompatible` | 4 = 2 pairs | 59PRr1+V45Dr1\_\_self (Rhodothermia) lineage RED **0.989308** → **genus**-compatible · F41PRr1+V45Dr1\_\_single (Thermoplasmatota SW-10-69-26) RED **0.736935** → **family**-compatible |
| automatic novel-rank assignments | **0** | `phylogeny_adjudication_summary.json` |
| priority manual taxonomic review | **0** | same |
Source: `M/rebuild_20260729/coassembly_novelty/taxonomic_novelty_evidence/legacy_phylogeny_adjudication_r232_v2/{candidate_phylogeny_adjudication.tsv,phylogeny_adjudication_summary.json,method.txt}`; mirrored in `reproducibility/r232_phylogeny_adjudication_summary.json`.

**Chimerism.** ⚠ **Correct a standing project claim.** "GUNC never ran" is true of the 14,426 refined bins and the 5,438 representatives, and M.2 is right to say no dedicated chimerism screen was run on them. But GUNC **1.1.1** *did* run on the 10 order-gap candidates against ProGenomes 3 (`gunc run --input_dir … --detailed_output --contig_taxonomy_output`), at `…/taxonomic_novelty_evidence/legacy_candidates_gunc/` (**V**, script `M/rebuild_20260729/scripts/118_gunc_legacy_high_rank.sbatch`). All 10 pass: clade-separation score **0.01–0.32**, all below GUNC's 0.45 threshold. The earlier "no artifacts" finding came from a `find -maxdepth 2`; the artifacts sit at depth 3–4 inside the PI's rebuild tree.

## Controls
There is **no control for catalogue construction itself** — no mock community, no spike-in, no simulated assembly, no held-out genome set. Say so plainly rather than implying one. What plays a control-like role, and who built it:
- **DAS_Tool consensus of two non-redundant binnings** (PI, `M/scripts/01_semibin2_single.sh` + Alejandra's `EQ_WGS` Snakemake). This is the only cross-check on bin composition, and it is a consistency check, not a control: the two tracks differ in model, mapper and min-contig length, so agreement is informative about bin boundaries but calibrates nothing.
- **GTDB R232 itself** is the reference against which "conventional" is measured. Its coverage is the yardstick, so the claim is release-relative by construction.
- **r226 → R232 re-run** (1,898 of 5,438 classifications changed) is a *stability* probe, not a control.
- **GUNC on the 10 candidates** (PI) is a genuine negative control for the only positive-sounding taxonomy claim in the section, and it comes back clean.
- Nothing in §1.1 has a blinded or matched-negative arm. Contrast §1.3 (3,684 blinded known-fold controls) and §1.4 (length- and size-matched annotated families) — reviewers who read the whole paper will notice §1.1 is the one section without one.

## Numbers
All rows below are **REGISTERED** in `notes/CATALOGUE_NUMBERS.md` at the stated line unless marked otherwise.
| value | quantity | source file | reg. line |
|---|---|---|---|
| 551 / 549 / 542 → 539 | cohort chain | metadata · assembly tree · `05_refine/sample_list.txt` | R:25, 26, 27 |
| 60 sites (61 labels); 117/166/268 compartments | design | `M/metadata/sample_metadata.tsv` | R:31, 32 |
| 1,043 km; 5 expeditions (17–20 Mar 2023 · 14 Jul 2023 · 17–21 Feb 2024 · 19–21 Aug 2024 · 12–15 Oct 2025) | design | M.1 | R:39 |
| 14,426 (4,225 HQ / 8,209 MQ / 1,992 LQ) | refined bins | `_checkm2_refined/quality_report.tsv` | R:59 |
| 12,433 of 12,434 | dRep input | `_drep/out/data_tables/{genomeInfo,Cdb}.csv` | R:60 |
| 5,438 (2,259 HQ + 3,179 MQ) | species reps | `_drep/out/data_tables/Wdb.csv` | R:61 |
| 51 named phyla (also 111 class, 280 order, 597 family, 1,303 genus, 184 species) | named-taxon census | `reproducibility/r232_named_taxon_census.tsv`; **V** — I re-derived 5,438 genomes / 51 phyla / 111 classes / 280 orders straight from `_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` | R:62 (phyla only) |
| 0/0/0 · 11 order (6 HQ) · 142 family (64 HQ) · 1,343 genus (530 HQ) · 3,755 species (1,555 HQ) · 187 fully assigned (104 HQ) | first empty rank | `reproducibility/r232_first_unassigned_census.tsv` | R:63 |
| 10 adjudicated → 0 proposals; 1 low-MSA in recovery | order gaps | adjudication v2 (path above; R lists the path as "to confirm" — **now confirmed**) | R:64 |
| 1,898 | r226→R232 reclassifications | `r232_taxonomic_novelty_summary.json` | R:65 |
| 37,949 = 16,333 + 21,616 | unrefined HQ+MQ union | see Basis table | R:272 |
| 304 | LQ coassembly-salvage bins (separate prospective cohort, 1 order gap) | `reproducibility/canonical_catalogue_summary.json` | not in R; in `verify.py` `scale.salvage` |
| 2 | bins at contamination exactly 10.0 | `_checkm2_refined/quality_report.tsv` (**V**, my awk) | R:66 footnote |

## Caveats — as written in the sources
- `R`: *"⚠ **60 is the site count; 61 is a count of field-notebook LABELS** … **Do not "correct" 60 to 61** — that was attempted and reverted."*
- `R` §1.1 footnote: *"HQ+MQ filter = completeness ≥50, contamination ≤10 *inclusive* — two bins sit at exactly 10.0 (pass ≤10, fail strict MIMAG <10); state the inclusive rule in Methods M.2/M.3."*
- `main.tex` M.2: *"No dedicated chimerism screen was run, so contamination control rests on CheckM2 and on DAS\_Tool's scoring of the two candidate binnings."* (carries a live `\todo`: *"optional: run GUNC on the refined bins"*).
- `PIPELINE_GUIDE.md` §4: *"⚠⚠ **TERMINOLOGY — locked decision, do not drift (2026-08-02).** **BOTH** binning runs are **single-sample** … **Always write "pretrained-model vs self-supervised", never "single-sample vs self-supervised."**"* And: *"the pretrained track binned only contigs ≥~2500 bp, dropping roughly a quarter of each assembly that the self-supervised track kept. That is why the two tracks recover different bins and why DAS_Tool merging them is not redundant."*
- `PIPELINE_GUIDE.md` §5: *"⚠ `checkM_method = lineage_wf` appears in the argument dump but is inert — quality came from the supplied `genomeInfo.csv`. Do not report CheckM1/lineage_wf in M.3."*
- `reproducibility/r232_claim_boundary.txt`: *"Empty-rank counts are candidate depths, not validated novel taxa; order candidates require GUNC, source-read support, dereplication context and explicit phylogeny."*
- GUNC `method.txt`: *"Reference representation below 0.3 is unscored, not a pass."* (54PRr2 sits at 0.32 — the nearest to unscored.)
- Adjudication `method.txt`: *"A singleton rank gap is never promoted automatically: it remains an unresolved single-genome lineage even when attachment context is strong."*

## Status
**REGISTERED and SHIPPED.** §1.1, M.2, M.3 and `tab:scale` are written on the current basis; `verify.py` carries ~25 machine-checked `scale.*` claims against the `reproducibility/r232_*.{tsv,json}` mirrors. A0a still reads `\status{RECOVERY RUNNING}` for the 11th candidate.

Dead bases this section must never be confused with (R:283–286): union-dRep **9,253** species / R226 / **56 phyla** / 311-2,830-8,978 novel-rank · older dRep **9,211** · single-sample GTDB on 16,333 at R226 with **44 / 38 phyla** · **15,469 / 15,477** single-sample bins (the strict `<10` variant; still the basis of the historical A1–A5 function tables) · everything from `_gtdbtk_refined/` (= GTDB-Tk 2.4.1 / r226).

⚠ **Three dRep runs sit on disk; only one is the paper's** (**V**, Wdb rows): `M/05_refine/_drep/out` → **5,438** (live; 5,438 `.fa` in `dereplicated_genomes/`) · `M/04_drep_merged/out` → **9,253** (dead union-dRep) · `M/04_drep_self_vs_single/out` → **9,211** (older dead run). The two dead runs have *empty* `dereplicated_genomes/`, so a genome-counting check against them returns 0 rather than erroring. Always dereference `05_refine/_drep/`.

## What can and cannot be claimed
**Can:** every one of 5,438 representatives is placed to class under GTDB-Tk 2.7.2/R232, across 51 named phyla; novelty is mutually-exclusively distributed 11/142/1,343/3,755 from order to species; the ten order gaps with adequate marker coverage are chimera-clean (GUNC ≤0.32) and, on de novo trees with R232 RED, yield zero automatic novel-rank assignments and zero manual proposals, the two multi-genome clades sitting at genus- and family-compatible depths; release provenance matters (1,898 reclassifications r226→R232); the 14,426 bins come from a consensus of two *non-redundant* SemiBin2 binnings, not one.

**Cannot:** that no novel high-rank lineage exists in the Rub' al-Khali — 11 of 5,438 remain order-gap candidates, one is unadjudicated, and 304 LQ salvage bins hold a further order gap; that the catalogue is chimera-screened (GUNC touched 10 genomes, not 14,426); that the HQ tier is MIMAG-high-quality in full (rRNA/tRNA never assessed); that "conventional" means anything beyond "well covered by GTDB R232"; that 37,949 is a MAG count; that the two binnings constitute a control; that the first-empty-rank counts are novel-taxon counts.

## Reviewer attack surface
- **"High-quality by MIMAG" is asserted but only half-tested.** MIMAG HQ requires 5S+16S+23S rRNA and ≥18 tRNAs; only completeness/contamination were applied (CheckM2's report has no rRNA/tRNA columns, and no barrnap/tRNAscan artifact exists for this catalogue). **The data do not answer this.** Either soften to "completeness and contamination tiers following MIMAG" or run the screen.
- **Contamination ≤10 inclusive is not MIMAG's `<10`.** Answered and immaterial: exactly 2 of 14,426 bins are affected (**V**) — but M.2/M.3 must state the inclusive rule, as the registry instructs.
- **No chimerism screen on the catalogue.** Partly answered: DAS_Tool consensus scoring plus clean GUNC on the 10 deepest candidates; not answered for the other 14,416 bins. M.2's own `\todo` concedes the point, so a referee will simply cash it.
- **"Taxonomically conventional" is reference-relative and close to circular.** Answered elsewhere, not here: §3.1/M.12 measures reference coverage directly. §1.1 should hand off explicitly rather than let the adjective stand alone.
- **DAS_Tool consensus could discard genuinely novel single-binner bins.** `--score_threshold 0.5` is the default and no sensitivity analysis exists. **Not answered.** The partial rebuttal on record is that the two tracks are non-redundant by construction (different min-contig length), so a novel bin recovered by one track is scored, not dropped unseen.
- **Four assembly-stage counts circulate: 549, 550, 552, 572.** Now fully reconciled (**V**): 572 = `*_contigs.fa` files in the assembly tree (the whole 573-library raw universe, not the study cohort); 552 = the files pyrodigal ran on; 550 = canonical samples with a contigs file (`has_contigs=True`); **549** = those minus the 0-byte `F5Dr1`. **Answered** — but M.2 should name `51Sr1` and `F5Dr1` rather than saying "549 of the 551 samples produced an assembly" and leaving it there.

## Open questions for Alejandra
1. **Register the two lost assemblies.** `51Sr1` (no contigs file) and `F5Dr1` (0-byte contigs) — verified this session, not currently in `R`. Add them to R:26 and name them in M.2? Note `51Sr1` and `F5Dr1` are also two of the three samples that R's footnote says lack genome profiles, so the same two failures propagate to the 548-sample CoverM cohort.
2. **MIMAG wording.** Soften the claim, or commission barrnap + tRNAscan-SE on the 5,438 (the PI already has both wired up, in an unrelated `enrichment_review_20260916/` workflow)?
3. **GUNC on the 14,426.** M.2 carries the `\todo`. Given GUNC 1.1.1 + the ProGenomes 3 reference are already staged in `rebuild_20260729/`, this is cheap. Run it, or fence the todo?
4. **Stale guidance to fix.** `PIPELINE_GUIDE.md` §5 and `DISCUSSION_NOTES.md` still present the r226/R232 label mix as unresolved and open; `_gtdbtk_r232/README.md` says it was fixed 2026-08-17 and I verified 4,847+591 = 5,438 on disk. Should the guide and notes be corrected, and does the Discussion still want its release-drift sentence?
5. **Same for "GUNC never ran"** — true of the catalogue, false of the 10 candidates. Correct the flat statement in `CORRECTIONS.md` and the guide?
6. **187 vs 184.** The A0 scaffold says 187 reps carry 184 named species labels; §1.1 prose gives only 187. Keep both, or is the scaffold's version the one to promote?
7. **Is the 48_S recovery still running?** A0a is tagged `RECOVERY RUNNING`; if it has landed, §1.1's "routed to targeted coassembly recovery" needs a result, and if it has not, the status tag needs a date.
8. **Registry fix:** R:64 gives "adjudication records (path to confirm)". The path is `M/rebuild_20260729/coassembly_novelty/taxonomic_novelty_evidence/legacy_phylogeny_adjudication_r232_v2/` — shall I register it?
