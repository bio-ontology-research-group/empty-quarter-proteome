# Revision review: dark-matter/MAG manuscript and supporting data

Date: 2026-07-29

Status: complete.

Scope: `paper/main.tex`, `reproducibility/claims.tsv`, `RQ/A3/`, `results/`,
the workbench linked at `workbench/`, and the analysis code that produced the
headline tables. I did not edit `paper/main.tex`.

This review uses the corrected raw/QC-read source:

```
/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads/<sample>_R1.fastq.gz
/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads/<sample>_R2.fastq.gz
```

It explicitly excludes the incomplete 150-sample directory
`/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_reads/processed/`.
The corrected source has 573 paired samples according to
`workbench/DATA_LOCATIONS.md`. The IBEX mount was not available in this
sandbox, and direct SSH failed at name resolution, so IBEX-only files and the
running co-assembly job could not be re-read directly. Findings below labelled
"recomputed" were reproduced from repository-local raw tables; IBEX-only
claims are identified as provenance checks still required.

## Executive assessment

The revision is not ready to act on unchanged. The taxonomy census, the
recomputed fold-threshold counts, and both negative group-1l hydrogenase
results match their supporting tables. The endemic-versus-EQ-wide dN/dS
contrast is also robust to family size and usable sequence count.

Four issues are blocking:

1. The replacement PICRUSt2 claim is invalid as implemented. The reported
   1.11-fold ratio thresholds PICRUSt2 but leaves the MAG-presence denominator
   unthresholded; the nominal 440 rows contain only 436 unique samples; and the
   MAG abundance matrix is empty for most samples.
2. The 0.39-versus-0.29 dN/dS control is selected by lexicographically taking
   the first non-dark families in each size bin. The resulting control
   representatives come from only 28 samples, versus 430 for the dark set.
   The dark-only evidence for purifying selection remains, but the
   dark-versus-control contrast is confounded.
3. The 21.4% "no homolog anywhere" intersection is supported only by prose.
   The query-level GEM and desert hit sets and the query FASTA are absent, and
   no-hit probability is strongly family-size dependent in the one raw hit
   table that is present.
4. The paper calls all 15,469 entries bacterial/archaeal MAGs even though 17
   are domain-unplaced, including six eukaryotic-like bins and four candidate
   large mobile elements. Marker-workflow counts (`bac120`/`ar53`) are being
   mistaken for final domain classifications.

The manuscript does not literally interchange the 51.9M MAG-protein catalogue
and the 753.3M assembly catalogue. It does, however, fail to define their
different scopes, and its novelty ladder skips the 23.416M representative
sequences that are the actual denominator for the 2.839M UniRef50 no-hit
count.

## Required sentence-level revisions

### 1. Define the catalogue and restore the missing denominator

**Quote** (`paper/main.tex:39`, `:164-166`, `:340-344`):

> "15,469 metagenome-assembled genomes; ~51.9M proteins"
>
> "51.9M proteins → 2.84M UniRef50 no-hit reps"

**Problem.** The 51,920,778 entries are Bakta proteins from the
`single_sample` 15,469-bin set. They are not the 753,305,533-protein
assembly-level catalogue, which covers all contigs, binned and unbinned, from
both binning strategies across 552 samples. More importantly, UniRef50 was
searched against 23,416,070 MMseqs2 representatives, not directly against all
51.9M proteins. Calling 2.84M "proteins" without the representative
denominator makes the dark fraction appear much larger than the measured
12.12%.

**Replace with:**

> Bakta predicted 51,920,778 proteins from the 15,469-bin
> `single_sample` catalogue. MMseqs2 linclust at 90% identity and 80%
> coverage reduced these to 23,416,070 representatives, of which 2,839,030
> (12.12%) had no UniRef50 hit. Of those sequence-dark representatives,
> 2,533,446 (89.2%) had no AFDB match in the initial
> ProstT5-predicted-3Di/FoldSeek screen.

Add one Methods sentence:

> The separate 753,305,533-protein assembly catalogue includes predicted
> genes on binned and unbinned contigs from both binning strategies and was
> not used for the analyses reported here.

**Evidence.**

- `docs/binning-and-annotation-scope.md`, "Two protein catalogs".
- `workbench/DATA_LOCATIONS.md`, "The protein / dark-matter chain".
- `workbench/results/soil_proteins/README.md`.
- IBEX: `$M/04_structure/pooled_proteins.faa`,
  `$M/04_structure/seq_vs_uniref50/{all_rep_ids,nohit_ids,hit_ids}.txt`.
- Separate assembly product:
  `$M/04_protein_catalog/assembly_proteins_pooled.faa.gz` and
  `assembly_proteins.mmseqsdb`.

### 2. Do not call the bac120/ar53 workflow counts domain classifications

**Quote** (`paper/main.tex:105-108`):

> "GTDB-Tk classified 15,469 MAGs (14,752 bacterial and 717 archaeal). No
> classified MAG defines a novel phylum or class..."

**Problem.** The 14,752 and 717 values are numbers routed through the bac120
and ar53 marker workflows. They include 14 and 3 domain-unplaced bins,
respectively. They therefore do not prove 14,752 bacteria plus 717 archaea.
Calling all 15,469 entries MAGs is also inappropriate after identifying six
eukaryotic-like bins and four candidate mobile elements.

**Replace with:**

> GTDB-Tk processed 15,469 HQ/MQ `single_sample` bins, routing 14,752
> through bac120 and 717 through ar53. Seventeen bins received no domain
> placement. Among the 15,452 domain-placed prokaryotic bins, none lacked a
> phylum or class assignment, whereas 12 lacked an order assignment, 297 a
> family assignment, 3,677 a genus assignment, and 10,815 a species
> assignment; 3,664 HQ bins met our operational SeqCode-candidate criteria.

Use "bins entering the MAG analysis" until the ten non-prokaryotic/noncellular
entries have been excluded and all affected denominators have been
recomputed. In the claims ledger, rename `scale.domains` from "classified
bacterial / archaeal MAGs" to "bac120 / ar53 workflow counts".

**Evidence.**

- `workbench/results/soil_gtdbtk/novelty_census.txt`.
- `workbench/results/soil_gtdbtk/novelty_summary.tsv`; its 17 `d__` rows
  comprise 14 bac120 and 3 ar53 records.
- `workbench/results/soil_gtdbtk/checkm2_all_samples.tsv`.

### 3. State what the 17 domain-unplaced bins are now known to be

**Quote** (`paper/main.tex:112-119`):

> "Follow-up composition and marker checks indicate that at least some may
> be eukaryotic sequence or large mobile elements rather than cellular
> genomes. ... the remaining domain-unplaced bins could include
> deep-branching lineages..."

**Problem.** This understates the completed triage and leaves the reader with
17 apparent novel-phylum candidates. Six are eukaryotic-like, four are
candidate large mobile elements with zero ribosomal proteins, and seven—not
17—remain uncharacterised. CheckM2 "MQ" is a prokaryotic-marker score and
should not be presented as biological validation of the six eukaryotic and
four mobile-element bins.

**Replace with:**

> Seventeen rhizosphere bins met CheckM2 MQ thresholds but received no
> GTDB-Tk domain placement. Six are eukaryotic-like (6.4–17.6 Mb; 11–29%
> coding density). Four small, gene-dense bins are candidate large mobile
> elements: all lack ribosomal proteins, are 88–98% hypothetical, and two
> carry unusually large tRNA complements; these four are excluded from the
> cellular MAG analysis and transferred to the mobilome study. Seven bins
> with intermediate coding density remain uncharacterised. None supports a
> novel-phylum claim.

The six eukaryotic-like IDs are the six largest `d__` rows in
`novelty_summary.tsv`; the four mobile-element candidates are
`33PRr2__SemiBin_18`, `48PRr2__SemiBin_146`,
`33PRr1__SemiBin_92`, and `4PRr2__SemiBin_192`. Deposit a small triage TSV
with category, genome size, coding density, marker counts, ribosomal proteins,
tRNAs, and final inclusion/exclusion disposition. Then recalculate every
"15,469 MAGs" and "51.9M MAG proteins" denominator or show that the exclusion
does not change the rounded values.

**Evidence.**

- Sizes, GTDB status, and marker workflow:
  `workbench/results/soil_gtdbtk/novelty_summary.tsv`.
- Coding density and CDS counts:
  `workbench/results/soil_gtdbtk/checkm2_all_samples.tsv`.
- Derived Bakta/marker triage: `workbench/CHANGELOG.md`, 2026-07-29 item 3.
- Primary Bakta files, to be cited in the deposited triage TSV:
  `$M/RQ/novel_phylum_check/bakta/<bin>/`.
- Refined marker outputs:
  `$M/02_taxonomy/gtdbtk_out_chunks/chunk_*/identify/`.

### 4. Reframe "no novel phyla" so it cannot mean "no taxonomic novelty"

**Quotes** (`paper/main.tex:41-42`, `:309`, `:458`):

> "GTDB-Tk recovers 44 phyla and no novel phylum or class, but leaves..."
>
> "Novel phyla / classes 0 / 0"
>
> "No novel phyla/classes, but substantial lower-rank novelty..."

**Assessment.** The sibling edit improved the abstract substantially by
placing the lower-rank census in the same sentence. The standalone table row
and the lead of Box 1 can still be detached from that context, and "novel
phyla" is not the operational measure: these are ranks left unassigned by
GTDB-Tk r226, not formally established new taxa.

**Use this framing consistently:**

> All domain-placed prokaryotic bins were assigned at phylum and class in
> GTDB-Tk r226; taxonomic novelty was nevertheless extensive below class
> (12 order-, 297 family-, 3,677 genus-, and 10,815 species-level
> unassigned bins).

Change the table label to:

> Domain-placed prokaryotic bins unassigned at phylum / class: 0 / 0

Change Box 1 to lead with:

> Extensive lower-rank novelty despite no phylum- or class-level gaps among
> domain-placed prokaryotic bins; 17 domain-unplaced bins resolved as six
> eukaryotic-like, four candidate mobile elements, and seven
> uncharacterised.

Also revise `plans/analysis-roadmap.md` and `plans/outline.md`, which still use
"no novel phyla" to imply that discovery is "not taxonomic".

**Evidence.** `workbench/results/soil_gtdbtk/novelty_census.txt`.

### 5. Keep the fold-threshold result, but qualify the join and model flow

**Quote** (`paper/main.tex:175-194`):

> "Explicit ESMFold structures were built for 68,527 ... We joined the
> 63,654 successfully scored models ... counting no hit as novel."
>
> "the novelty rate remains 16–24% ... argues against low-confidence models
> alone manufacturing the novelty signal."
>
> "TM-align scores are unreliable for short structures..."

**Assessment.** The threshold and length numbers match
`RQ/A3/fold_validation_summary.tsv` exactly. The conclusion is directionally
reasonable but is not yet a matched sensitivity analysis:

- 68,527 built models versus 63,935 pLDDT-shard records leaves 4,592 built
  models unaccounted for; only 63,654 are `ok`. The manuscript does not explain
  that flow.
- `best_afdb_tm.tsv` has 58,235 query rows and `best_pdb_tm.tsv` 60,055.
  `finalize_tm.py` assigns every missing row TM=0. This is valid only if a
  completed-query manifest proves that absence means a genuine no-hit rather
  than a missing conversion/search result.
- Normalised pLDDT is stored on a 0–1 scale. Call it "mean pLDDT/100" or report
  the conventional thresholds 40, 50, 60, 70, and 80.
- The code uses non-overlapping bins `100–<150`, `150–<250`, and
  `250–<400`, and the long count is length **≥100 aa**. Use those exact
  boundaries everywhere.
- Stability across confidence thresholds is evidence against confidence
  alone explaining the calls, but length composition also changes across
  those thresholds. It is not a length-controlled result.
- "TM-align scores are unreliable" is too categorical. The defensible
  statement is that fixed TM-score classification and search sensitivity are
  length dependent for short proteins.

**Replace the interpretation with:**

> Across thresholds retaining at least 5,709 models, the dual-no-match rate
> remained 16.4–24.0%, which is inconsistent with model confidence alone
> explaining all calls. This sensitivity series is not length matched,
> however. At mean pLDDT/100 ≥0.60, 2,536/9,660 proteins shorter than 100 aa
> were called dual-no-match (26.3%), compared with 300 total calls among
> proteins ≥100 aa (2.7–5.4% by length bin). These remain candidate
> structural novelties pending a known-fold false-novelty calibration.

Add a model-flow row to the ladder:

> 68,808 candidate representatives → 68,527 models reported built →
> 63,935 pLDDT-shard records → 63,654 scored `ok` models.

Do not call any of these a validated novel-fold count.

**Evidence.**

- `RQ/A3/fold_validation_summary.tsv`.
- `RQ/A3/finalize_tm.py`, especially the missing-row default in
  `summarise()`.
- `RQ/A3/best_afdb_tm.tsv` and `RQ/A3/best_pdb_tm.tsv`.
- `workbench/DATA_LOCATIONS.md`, protein/dark-matter chain.

#### Exact known-fold false-novelty control to run

1. Build an independent positive-control pool of single-domain proteins with
   experimentally established CATH or ECOD folds. Record sequence, domain
   boundaries, experimental structure, and fold label.
2. Remove the query itself and its close sequence-homology cluster from the
   AFDB/PDB search references, while retaining at least one independently
   established remote member of the same fold. Otherwise exact self-hits make
   the false-novelty estimate trivially zero.
3. Fold every control sequence with the identical ESMFold version, settings,
   length limit, and hardware path used for the EQ candidates. Run two arms:
   experimental control coordinates (search/alignment error alone) and
   ESMFold coordinates (search plus modelling error).
4. Match or reweight controls to the candidate distribution using 10-aa
   length bins and 0.02-wide mean-pLDDT/100 bins, with at least five controls
   per candidate where possible. Match additionally on single- versus
   multidomain status, predicted disorder/low complexity, and transmembrane
   status, or restrict the primary analysis to single-domain soluble
   proteins.
5. Run exactly the same FoldSeek/TMalign commands, database snapshots,
   `--alignment-type`, `--max-seqs`, score extraction, missing-hit rule, and
   dual-TM<0.5 classification used for the candidates.
6. Report the false-dual-no-match rate with binomial or bootstrap 95%
   confidence intervals in every length × pLDDT stratum. Reweight those rates
   to the candidate distribution to estimate the expected false count and a
   lower bound on the positive predictive value. Predefine the ≥100-aa,
   pLDDT/100≥0.60 primary stratum.
7. Deposit the complete query manifest, per-database best-hit rows including
   explicit `SEARCHED_NO_HIT` states, software/database versions, and all
   stratum assignments.

### 6. Retain dark-family purifying selection, but withdraw the current control interpretation

**Quote** (`paper/main.tex:198-202`, figure caption `:424-425`):

> "dark families have median dN/dS = 0.39 << 1 (size-matched annotated
> control 0.29) ... refutes the random/erroneous-ORF null"

**Problem.** The dark median is robust to removing small families, but the
control is not a valid general "annotated" control. In
`code/ibex/93_darkfam_member_cds.sh`, controls are sorted and the first `n`
identifiers in each log2 family-size bin are taken. Family identifiers begin
with sample IDs, so the final controls have representatives from only 28
samples, compared with 430 for dark families. They are also merely
`not in isdark_ids.txt`, not necessarily functionally annotated.

Reanalysis of `results/darkmatter_ecology/dnds.tsv`:

| minimum total family members | dark families | dark median dN/dS |
|---:|---:|---:|
| 3 | 22,767 | 0.391 |
| 5 | 8,769 | 0.418 |
| 10 | 2,479 | 0.449 |
| 20 | 671 | 0.489 |
| 30 | 302 | 0.499 |

Thus the conclusion that recurrent, alignable dark ORFs are under purifying
selection is robust to family size. The 0.39-versus-0.29 effect is not:
restricting dark representatives to the same 28 source samples gives 0.302
versus 0.286 (`p=0.0063`), so most of the raw gap tracks control composition.
Protein-length robustness cannot be assessed because `dnds.tsv` lacks
representative length and retained codon-alignment length, and the member
FASTAs are not archived locally.

**Replace with:**

> Among 22,767 recurrent dark families that yielded at least one usable
> NG86 comparison, median dN/dS was 0.391 and remained below 0.50 after
> requiring at least 30 family members, consistent with pervasive purifying
> selection. The present non-dark comparison is sample-confounded and is not
> interpreted pending rematching. These data argue against, but do not by
> themselves refute, an erroneous-ORF explanation.

Rerun the control with a fixed random seed and multiple matched draws, matching
on exact or narrow family size, protein/alignment length, mean pairwise
identity, source sample/compartment/trip, representative MAG taxonomy and
quality. Use a matched model or mixed model with sample/taxon random effects,
not only pooled medians. Preserve the existing ascertainment in Methods:
CDS ≥90 nt, pairwise retained alignment ≥30 codons, non-identical pairs,
`0<dS≤2`, and at most 80 pairs per family.

**Evidence.**

- `results/darkmatter_ecology/dnds.tsv`.
- `results/darkmatter_ecology/darkfam_meta.tsv`.
- `code/ibex/93_darkfam_member_cds.sh`, control selection.
- `code/ibex/94_darkfam_dnds.sh`, pair filters and NG86 calculation.

### 7. The endemic-versus-wide selection result is family-size robust, but not length tested

**Quote** (`paper/main.tex:199-201`):

> "endemic (≤3 sites) families are more constrained than EQ-wide (≥30
> sites) (median 0.332 vs 0.420, MWU p<10^-67)"

**Assessment.** This result reproduces: 1,809 endemic families have median
0.3324; 10,498 EQ-wide families have median 0.4201; two-sided MWU
`p=5.32e-69`. The contrast remains after requiring progressively larger
families:

| minimum usable sequences | endemic median | EQ-wide median | p |
|---:|---:|---:|---:|
| 3 | 0.332 | 0.420 | 5.3e-69 |
| 5 | 0.368 | 0.442 | 3.8e-22 |
| 10 | 0.384 | 0.476 | 2.3e-10 |
| 20 | 0.363 | 0.510 | 3.1e-7 |
| 30 | 0.432 | 0.523 | 5.8e-4 |

Within every total-family-size stratum from 3 through ≥41 members, the endemic
median is lower. Exact matching on total family size and usable sequence count
also retains the difference (0.333 versus 0.414), and partial rank correlation
of site occupancy with dN/dS controlling both counts is `r=0.169`,
`p=1.9e-145`.

The manuscript may retain the median contrast, adding:

> The association was robust to total family size and the number of usable
> aligned sequences; robustness to protein length could not be tested from
> the current output.

Fix the sign explanation in `code/phase_d_endemism.py` and
`results/darkmatter_ecology/ENDEMISM.md`: the observed correlation is
**positive** (`r=+0.176`), which means occupancy increases with dN/dS and is
consistent with widespread families being less constrained. The current
parenthesis incorrectly starts "negative =>".

**Evidence used for the reanalysis.**

- `results/darkmatter_ecology/dnds.tsv`.
- `results/darkmatter_ecology/darkfam_meta.tsv`.
- `results/darkmatter_ecology/darkfam_abundance.tsv`.
- `results/darkmatter_ecology/env_matrix.tsv`.

### 8. Replace "no homolog anywhere" and remove the unsupported plateau inference

**Quote** (`paper/main.tex:257-264`, table `:371-382`, figure caption
`:449-452`):

> "21.4% have no homolog anywhere sequenced"
>
> "Endemism plateaus at 32–47% at the deepest references (the signature of
> genuine endemism, not undersampling)"

**Problem.**

- Only MGnify90, GEM, and five selected desert catalogues were searched.
  "Anywhere sequenced" is not supported, and SPIRE is still absent.
- The 5,011-family intersection appears only in
  `results/global_specificity/README.md`. The query-level GEM and desert hit
  files needed to recompute the union are not deposited.
- The catalogues are neither nested nor comparable merely by protein count.
  The non-monotonic 31.9–47.4% rates do not constitute a saturation curve.
- Searches use DIAMOND `--sensitive --evalue 1e-3` without a query-coverage
  or identity threshold. E-values and hit probability are length dependent.
- Family size is already a strong confound in the deposited MGnify result:

| family members | MGnify no hit |
|---:|---:|
| 3 | 4,255/9,807 (43.4%) |
| 4 | 1,727/4,721 (36.6%) |
| 5 | 810/2,508 (32.3%) |
| 11–20 | 386/1,458 (26.5%) |
| 41–80 | 19/152 (12.5%) |
| ≥81 | 2/27 (7.4%) |

Protein-length robustness and family-size robustness of the 5,011 intersection
cannot be tested because `query.faa` and the per-query GEM/desert hits are
missing.

**Replace with:**

> Under DIAMOND `--sensitive` at E<10^-3, 8,552/23,382 representatives
> lacked an MGnify90 hit, 11,075 lacked a GEM hit, and 7,469 lacked a hit in
> the union of five desert catalogues. A summary reports that 5,011 lacked a
> hit in all three catalogue sets, but this intersection remains
> catalogue- and threshold-specific and awaits query-level archival and
> length/family-size sensitivity analysis.

Until that archival and analysis are complete, remove "no homolog anywhere",
"genuine endemism, not undersampling", and "plateau". Use "no detected hit in
the searched catalogue sets". Treat 21.4% as an upper bound that may decrease
as additional catalogues are searched, not a final global-endemism fraction.

**Evidence.**

- Raw MGnify hits:
  `results/darkmatter_ecology/endemism/mgnify_hits.m8`.
- Query metadata and family sizes:
  `results/darkmatter_ecology/darkfam_meta.tsv`.
- Aggregate-only GEM/desert counts:
  `results/global_specificity/{global,desert}_specificity_summary.tsv`.
- Unsupported intersection prose:
  `results/global_specificity/README.md`.
- Search definitions:
  `workbench/scripts/ibex/{98b_mgnify_search,131_global_specificity,131b_desert_specificity}.sh`.

### 9. Withdraw the current PICRUSt2 replacement result; it is not a shared-threshold analysis

**Quotes** (`paper/main.tex:50-53`, `:231-240`, table `:355-368`, Box 1
`:463-464`):

> "The apparent 2.4× excess ... falls to 1.11× at a shared
> relative-abundance threshold of 10^-5."
>
> "the rebuilt 440-sample comparison"
>
> "93% KO presence-recall"
>
> "56.7% of MAG proteins carry no KO"

**Problem.** The sibling edit reproduces numbers in the new TSV but the code
does not implement the stated comparison.

1. In `workbench/scripts/analysis/174_picrust2_rebuild.py`, PICRUSt2 and the
   MAG relative-abundance matrix are thresholded at `t`, but MAG presence is
   always counted with `v>0`. The `ratio_picrust_over_mag_presence` therefore
   divides a thresholded PICRUSt2 numerator by the same unthresholded 4,346.48
   MAG-presence denominator at every threshold. At `10^-5`, the actual
   same-threshold relative-abundance ratio written in the adjacent column is
   14.633, not 1.115.
2. The "440 samples" are 440 PICRUSt2 columns but only 436 unique sample IDs.
   `F41Dr3`, `F58Sr3`, `T1Dr1`, and `F1PRr2` each occur twice.
3. `mag_ko_relab_by_sample.tsv` has any non-zero abundance for only 84/532
   WGS samples. Only 310/15,469 MAG rows have finite CoverM abundance;
   15,159 are `nan`. Consequently, 381/440 rebuilt rows have
   `mag_kos=0` at `10^-5`.
4. The original 93% presence recall is also a zero-threshold result:
   PICRUSt2 has 10,543 non-zero KOs in every one of the 436 rows. Calling this
   recovery of a "housekeeping core" is not justified without a detection
   sensitivity analysis.
5. The reported signature-function "abundance" correlations compare PICRUSt2
   predicted abundance with unweighted MAG KO copy counts
   (`magcount_s`), not CoverM-weighted community abundance.
6. The 56.7% denominator mixes 55,380,378 CheckM2
   `Total_Coding_Sequences` with 23,967,112 KO-annotated EggNOG rows. These
   are not the paper's 51,920,778 Bakta proteins. Five bins even have more
   EggNOG rows than CheckM2 CDS, and two have more KO rows than CheckM2 CDS,
   producing impossible negative "no annotation" fractions.

**Replace the entire C1 numerical interpretation, provisionally, with:**

> PICRUSt2 assigned non-zero abundance to essentially all 10,543 predicted
> KOs in every matched sample; therefore the zero-threshold 2.40-fold
> richness ratio and 93% presence recall are detection-rule dependent. The
> current 1.11-fold sensitivity estimate thresholds only the PICRUSt2
> numerator and is not retained. A deduplicated comparison with complete MAG
> abundance and prespecified detection rules is required before reporting a
> corrected ratio, recall, or abundance concordance.

The valid conceptual point—that KO-based prediction cannot represent genes
without a KO—can remain, but recompute its fraction using one gene-calling
catalogue for both numerator and denominator. Do not call 56.7% a fraction of
the 51.9M Bakta catalogue until it is calculated from that catalogue.

The rerun should:

1. collapse duplicate PICRUSt2 columns to 436 unique biological samples before
   matching;
2. repair the MAG abundance-name join and require finite CoverM abundance for
   all included MAGs/samples;
3. predefine two separate estimands: binary encoded presence versus predicted
   detection, and relative community abundance versus relative predicted
   abundance;
4. apply an explicitly justified threshold to the relevant side(s), reporting
   asymmetric detection limits honestly rather than calling them shared;
5. report sample-level distributions and bootstrap confidence intervals, not
   only ratios of means;
6. recompute KO-negative fractions within the exact Bakta or EggNOG query
   universe used in the paper.

Set C1 status to `REQUIRES RERUN`, remove the current table/figure numbers, and
do not regenerate prose from `results/picrust2_bias/SUMMARY.md` until the code
is corrected.

**Evidence.**

- Threshold bug:
  `workbench/scripts/analysis/174_picrust2_rebuild.py`, lines 162–188.
- Contradictory columns:
  `results/picrust2_bias/rebuild/threshold_sweep.tsv`.
- Duplicate/zero rows:
  `results/picrust2_bias/rebuild/per_sample_rebuilt.tsv`.
- Missing CoverM values:
  `results/picrust2_bias/mag/mag_annotation_stats.tsv` and
  `mag_ko_relab_by_sample.tsv`.
- Zero-threshold recall and unweighted function comparison:
  `code/picrust2_bias_analysis.py`.
- Mixed denominator:
  `results/picrust2_bias/dark_fraction.tsv` and
  `workbench/results/soil_gtdbtk/checkm2_all_samples.tsv`.

### 10. Define all sample universes and exclusions

**Quote** (`paper/main.tex:39`, table `:304-307`):

> "60 sites × 5 expeditions"
>
> "Soil samples 547"

**Problem.** "547 soil samples" is the row count of
`results/darkmatter_ecology/env_matrix.tsv`, not the raw-read count. The
project currently contains several legitimate nested universes: 573 paired
read samples, 571 primary assemblies in one audit, 552 canonical
assembly-catalogue samples, 551 dark-family abundance samples, 547
environment-matrix samples, 535 samples contributing the single-sample
protein catalogue, and 527 samples in 126 corrected co-assembly groups. The
workbench changelog also records inconsistent shorthand exclusions
(`17 C + 3 controls` in one entry; `16 C + 4 controls` in another). A reader
cannot reconstruct the transitions.

**Change the table row to:**

> Samples in dark-family environmental matrix: 547

Add a supplementary sample-flow table with one row per stage, manifest path,
inclusion rule, exclusion counts/reasons, and final sample IDs. Begin Methods
with:

> Paired QC reads were available for 573 samples at
> `soil/wgs_assemblies/processed_reads/`; analyses used explicitly
> documented nested subsets for assembly, binning, MAG annotation,
> abundance, environment matching, and co-assembly.

Do not use partial output from co-assembly job 49568250 in this revision.

**Evidence.**

- Correct read source and 573 count: `workbench/DATA_LOCATIONS.md`.
- Ecology subsets: `docs/dark-matter-pipeline.md`.
- Protein/assembly scopes: `docs/binning-and-annotation-scope.md`.
- Historical counts/exclusions: `workbench/CHANGELOG.md`.

### 11. Keep both group-1l aridity negatives

**Quote** (`paper/main.tex:128-132`):

> "Two independent aridity tests were also negative..."

**Assessment.** This edit matches the data and should remain. The within-EQ
test is `rho=-0.1716`, permutation `p=0.19610`, `n=58`; the cross-biome
fractions are non-monotonic. Keep the distinction that the atmospheric
enzyme census is retained while the aridity-gradient hypothesis is not.

**Evidence.**

- `results/atmospheric/aridity/aridity_summary.tsv`.
- `results/atmospheric/aridity/cross_biome_hydrogenase.tsv`.

## Audit of the 2026-07-29 sibling edits

| Area | Verdict | Mismatch or required qualification |
|---|---|---|
| Fold threshold counts | Matches | All six denominators, counts, and rates match the corrected values. |
| Fold length strata | Matches numerically | Labels should show upper-exclusive bounds; 300 is for length ≥100 aa, not strictly >100 aa. |
| Known-fold control | Correctly marked `NOT_RUN` | A control protocol is specified above; no validated-fold count is available. |
| Fold analysis universe | Incomplete | 4,592-model gap and missing-hit-as-zero assumption are unaudited. |
| Taxonomic rank census | Matches | 12/297/3,677/10,815 and HQ counts, plus 3,664, reproduce. |
| 14,752 / 717 domains | Does not match the stated meaning | They are bac120/ar53 workflow counts, not final bacterial/archaeal classifications. |
| Domain-unplaced bins | Manuscript incomplete | Must say 6 eukaryotic-like + 4 mobile + 7 uncharacterised and remove the four from this paper. |
| PICRUSt2 2.40 → 1.11 | Number reproduced, interpretation false | Threshold is not shared; sample duplication and missing MAG abundance invalidate the replacement. |
| PICRUSt2 93% recall | Numerically reproduced, biologically threshold-dependent | PICRUSt2 has 10,543 non-zero KOs in every sample. |
| PICRUSt2 56.7% | Arithmetic reproduced, denominator inconsistent | CheckM2 CDS and EggNOG counts are mixed; it is not a measured fraction of the 51.9M Bakta catalogue. |
| Hydrogenase aridity negatives | Matches | Both independent negatives are supported. |
| Claims verifier | Misleading success | `46 pass` checks file arithmetic, not semantics; all new A3 claims are `TABLE`, and invalid PICRUSt2 arithmetic still passes. |

## Robustness summary for the three requested dark-ecology claims

| Claim | Family-size robustness | Protein-length robustness | Conclusion |
|---|---|---|---|
| Dark dN/dS 0.39 versus control 0.29 | Dark median remains <0.50 through ≥30-member families, but the control is sample-confounded | Not testable from deposited output | Retain dark dN/dS<1; withdraw the interpreted dark/control gap pending rematching |
| Endemic 0.332 versus EQ-wide 0.420 | Robust across total family size, usable sequence count, exact matching, and partial rank analysis | Not testable from deposited output | Retain with a stated length limitation |
| 21.4% no hit in three catalogue sets | Intersection not reproducible; MGnify no-hit rate is strongly family-size dependent | Not testable because query FASTA/lengths are absent | Withhold as a global-endemism headline; report catalogue-specific aggregate rates only |

To complete the length audit, add to every family-level analysis table:
representative amino-acid length, median member length, retained codon-alignment
length, fraction gapped/low-complexity, and domain count. Then repeat the
contrasts in prespecified length strata and with length-matched controls.

## Catalogue identity and provenance

No direct 51.9M/753M substitution occurs in `paper/main.tex`. The paper uses
the 51.9M `single_sample` MAG-bin proteome throughout. The following wording
still needs correction:

- "proteins (total)" must be "Bakta proteins in the `single_sample` bin
  catalogue"; 753M is the actual assembly-level total.
- "2.84M proteins" must be "2.839M of 23.416M representative sequences".
- "56.7% of MAG proteins" uses a third universe—55.380M CheckM2 CDS—and cannot
  be placed on the 51.9M ladder.
- The refined 14,426-bin catalogue (HQ 4,225 / MQ 8,209 / LQ 1,992) is distinct
  from the 15,469-bin `single_sample` catalogue. The manuscript must say which
  is used at every stage. Co-assembly of the 304 salvageable LQ bins is still
  running and is not current evidence.

The novelty ladder should explicitly show units:

```
51,920,778 Bakta protein records
  -> 23,416,070 90%-identity representatives
  -> 2,839,030 UniRef50-no-hit representatives
  -> 2,533,446 initial AFDB-no-hit representatives
  -> 334,810 multi-member sequence-novel families
  -> 68,808 candidate novel-fold family representatives
  -> confidence/length sensitivity on 63,654 scored explicit models
```

## Missing, moved, vague, or stale supporting evidence

### Claims ledger pointers that are not resolvable files

`reproducibility/claims.tsv` still contains free-text sources rather than
evidence paths:

| Claim | Current source | Required source |
|---|---|---|
| `scale.proteins` | `workbench (Bakta; IBEX_SOIL_PROGRESS)` | `$M/04_structure/pooled_proteins.faa` plus a count manifest |
| `scale.samples` | `sampling design metadata` | a committed stage-by-stage sample manifest |
| `b1.dark` | `workbench/results/darkmatter (mmseqs/UniRef50)` | `$M/04_structure/seq_vs_uniref50/nohit_ids.txt` |
| `b1.structnovel` | `workbench (ProstT5+FoldSeek vs AFDB)` | `$M/04_structure/prostt5_foldseek_orix/merged_afdb_hit_ids.txt` plus the complete query manifest |
| `b1.families` | `workbench novel-fold clustering` | `$M/04_structure/novelty_analysis/novel_family_summary.tsv` |
| `b1.folds` | `workbench is_novel_family` | `$M/06_esmfold_novelfold/novelfold_reps.txt` |
| `b1a.structures` | `ORIX ESMFold run` | an explicit model-success/failure manifest under `$M/06_esmfold_novelfold/` |

The local summary `workbench/results/soil_proteins/README.md` resolves the
intended provenance, but the ledger should cite the actual file rather than
requiring a reader to discover it.

### Wrong or insufficient sources

- `c3.mgnify` cites `darkfam_meta.tsv`; the hits are in
  `results/darkmatter_ecology/endemism/mgnify_hits.m8`. The verifier currently
  hard-codes 8,552 and only divides by the metadata row count.
- `c3.intersection` cites prose in `results/global_specificity/README.md`.
  Deposit and cite the query-level GEM and desert hit tables plus code that
  computes the intersection.
- `scale.domainmq` cites CheckM2 alone, but CheckM2 does not prove GTDB
  marker-workflow assignment or rhizosphere origin. Cite
  `novelty_summary.tsv` plus the sample metadata.
- `scale.mags` is "verified" by counting `dark_fraction.tsv`, a downstream
  PICRUSt2 table. Cite the canonical MAG/bin manifest.
- All A3 threshold claims are `TABLE`; `reproducibility/verify.py` does not
  rerun `finalize_tm.py` or validate its join.

### Required raw artifacts absent locally

- `results/darkmatter_ecology/endemism/query.faa`.
- Query-level `gem_hits.m8` and `deserts_combined_hits.m8`.
- A SPIRE result/hit table.
- The dN/dS family member FASTAs and per-family retained alignment lengths.
- `results/darkmatter_ecology/darkfam_presence.tsv`, named as an output in
  `docs/dark-matter-pipeline.md`, is absent; the abundance table is present.

Large files need not be committed to Git, but checksummed manifests and a
stable archive path are required.

### Stale documents that can reintroduce withdrawn claims

- `handover.md`, `RQ.md`, and `RQ/INDEX.md` still call 10,007 a validated
  novel-fold floor.
- `reproducibility/REPORT.md` is stale: it records the older 37-pass/14-doc
  ledger, whereas the current verifier prints 46 pass/35 documented. Regenerate
  only after correcting the semantic errors above.
- `results/picrust2_bias/{README,SUMMARY}.md` repeat the invalid "shared
  threshold", 440-sample, 93%-core, and 56.7%-of-all-proteins statements.
- `results/global_specificity/{README,MULTICATALOG}.md` repeat "anywhere on
  Earth" and the unsupported plateau inference.
- `results/darkmatter_ecology/ENDEMISM.md` contains the wrong sign explanation
  for its positive correlation.

## Additional manuscript-wide revisions

1. Change the abstract's "We argue—and test with explicit nulls" because the
   known-fold null is not run, the dN/dS control is confounded, and the
   multi-catalogue intersection is not reproducible. Use "we evaluate with
   sensitivity analyses and explicit negative results" until the controls are
   complete.
2. Change B2's heading "It is real, not spurious ORFs" and Box 2's
   "defeated by purifying selection" to conditional language: recurrent,
   alignable dark ORFs show purifying selection. The current pipeline excludes
   singletons, short/unusable alignments, identical pairs, and saturated pairs.
3. Change B3's "genuine desert-endemic reservoir" to
   "catalogue-specific no-hit reservoir" until length/family-size-matched global
   searches are complete.
4. Expand the pointer-level Methods. At minimum, give every catalogue
   definition and denominator, sample manifest, database snapshot, software
   version, search parameters, gene caller, matching rule, missing-data rule,
   and statistical unit.
5. Add an analysis-universe flow figure or table. It should keep separate the
   15,469 `single_sample` bins, 14,426 refined bins, 5,438 dereplicated refined
   MAGs, 51.9M MAG-bin proteins, and 753.3M assembly proteins.

## Recommended order for the revision pass

1. Remove/triage the ten non-prokaryotic/noncellular bins and freeze all
   catalogue/sample manifests.
2. Correct the claims-ledger paths and regenerate the full catalogue ladder.
3. Repair and rerun PICRUSt2 before retaining any C1 number or figure.
4. Rematch the dN/dS control and add protein/alignment lengths.
5. Archive query-level global-specificity results and run length/family-size
   sensitivity.
6. Run the known-fold false-novelty calibration.
7. Apply the sentence replacements above, regenerate figures and
   `reproducibility/REPORT.md`, run ontology-independent manuscript QC, and
   inspect the final diff.
