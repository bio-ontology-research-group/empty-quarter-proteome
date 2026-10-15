# Quantitative-claim and provenance audit of the two Empty Quarter manuscripts

Date: 2026-07-29

## Scope and verdict key

This report audits every quantitative claim in the abstract, Results, figures,
and tables of:

1. `empty-quarter-mags/paper/main.tex`, with
   `empty-quarter-mags/reproducibility/claims.tsv`;
2. `empty-quarter-phages/manuscript/main.tex`, with
   `empty-quarter-phages/paper/outline.md`.

The audit recomputes each value from its source file where possible, checks its
denominator and aggregation level, and checks repeated values for consistency.
The manuscripts were not edited.

Findings use the requested severity classes:

- **(a) WRONG** — the printed number does not agree with recomputation.
- **(b) UNFOLLOWABLE PROVENANCE** — the number is right, but the manuscript or
  ledger points only to a repository, directory, or otherwise non-resolving
  source rather than the file needed to reproduce it.
- **(c) MISSING CONTEXT** — the number is right, but its denominator or
  aggregation level is not stated.
- **(d) INTERNAL INCONSISTENCY** — the same quantity has different values in
  different places.

`PASS` means that the recomputed and printed values agree at the displayed
precision and that the claim states the denominator and aggregation level
needed to interpret it. A claim can have more than one finding class.

## Audit log and source controls

### Check 0 — canonical data-location map

**Status: completed.** Read
`/home/leechuck/Public/software/empty-quarter-workbench/DATA_LOCATIONS.md`
before resolving any data pointer. The audit therefore treats the following as
distinct products and aggregation levels:

- 51,920,778 Bakta proteins predicted from MAGs versus 753,305,533 proteins in
  the assembly-wide catalogue;
- 23,416,070 90%-identity representatives, of which 2,839,030 (12.1%) are
  UniRef50 no-hits;
- the 2.84-million-representative structural-novelty chain versus the separate
  tier1/tier2 protein analyses;
- 334,810 multi-member sequence-dark families versus 156,743 tier2
  novel-fold families.

The canonical read source is
`/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads/`
(573 samples), not the truncated
`.../soil/wgs_reads/processed/` tree (150 samples). A source recorded merely as
`workbench` or `workbench/results/darkmatter` is treated as a provenance defect,
not as evidence that the associated value is numerically wrong.

### Check 1 — manuscript and companion-file inventory

**Status: completed.** Inspected both complete LaTeX manuscripts and their
companion files. The MAG ledger contains 81 claim entries (46 designated for
local computation and 35 designated only as table/provenance checks), but it
does not by itself exhaust all numbers printed in the Results, tables, and
figure captions; this audit therefore also includes unregistered values such as
the photolyase test statistic, p-values, structural-model threshold strata,
table-only counts, and plotted values. The phage paper has no equivalent
machine-readable claim manifest: `paper/outline.md` mixes evidence notes,
provisional statements, and source pointers, so every printed value in
`manuscript/main.tex` was inventoried independently.

Repeated occurrences are grouped only when they denote the same object,
denominator, and analysis. Values from different aggregation levels remain
separate even when the manuscript uses the same noun for them.

### Check 2 — MAG repository verifier and exact BGC arithmetic

**Status: completed.** `python3 reproducibility/verify.py` reported 46/46
registered `COMPUTE` claims as passing and 35 claims as documented rather than
recomputed. That output is only a screening result: its tolerance accepts the
printed terpene value even though exact rounding does not.

Source `results/bgc_product_summary.tsv` contains 9,984 MAGs with a terpene BGC.
With the catalogue denominator of 15,469 MAGs from
`results/picrust2_bias/dark_fraction.tsv`,

`100 × 9,984 / 15,469 = 64.541987...%`,

which rounds to **64.5%**, not the printed **64.6%**. The NRPS and ectoine rows
recompute as `100 × 3,914 / 15,469 = 25.3022...%` (**25.3%**) and
`100 × 2,626 / 15,469 = 16.9759...%` (**17.0%**), respectively.

### Check 3 — defense carriage by plasmid mobility class

**Status: completed from the complete canonical IBEX files.** Retrieved and
checked:

- `/ibex/scratch/projects/c2014/rob/metagenomics/06_mobilome/mge_mobility/potu_mobility.tsv`
  (312,567 data rows);
- `/ibex/scratch/projects/c2014/rob/metagenomics/06_mobilome/defense_on_mobilome/defense_on_mobilome.tsv`.

Filtering the latter to `element_type == plasmid` gives 9,780 system calls on
8,221 distinct pOTUs. These counts agree exactly with the joined
`n_defense_systems` column in `potu_mobility.tsv`. The three mobility classes
partition all 312,567 pOTUs exactly:

| Class | pOTUs | Defense-positive pOTUs | Defense-system calls |
|---|---:|---:|---:|
| Conjugative | 33 | 11 (33.3%) | 17 |
| Mobilizable | 466 | 31 (6.7%) | 34 |
| No detected mobility marker | 312,068 | 8,179 (2.6%) | 9,729 |

I independently repeated the recorded length-matching procedure: for each
mobility-typed pOTU, sample a no-marker pOTU within ±20% of its length; 2,000
draws with seed 20260726; use the median matched-positive rate; one-sided
resampling p-value with the `(extreme + 1)/(draws + 1)` correction. Results:

| Class | Observed | Median length-matched no-marker control | Ratio | Resampling p |
|---|---:|---:|---:|---:|
| Conjugative | 11/33 = 33.3% | 3/33 = 9.1% | 3.6667 → **3.7×** | 1/2,001 = **0.0005** |
| Mobilizable | 31/466 = 6.7% | 16/466 = 3.4% | 1.9375 → **1.9×** | 1/2,001 = **0.0005** |

`manuscript/main.tex` correctly says that length matching was used and explicitly
states that `n=33` is small, but it omits both effect sizes and both p-values.
`paper/outline.md` goes further and says that the estimates await complete-source
re-verification; that status is now stale. Neither file cites the two concrete
input files above or a saved resampling-output file. The numerical result is
right, but its printed provenance is unfollowable.

### Check 4 — protein-catalogue and dark-matter ladder

**Status: completed.** Counts were recomputed from the canonical products rather
than inferred from rounded manuscript values:

| Product / aggregation level | Recomputed count |
|---|---:|
| Bakta proteins from MAGs | 51,920,778 |
| 90%-identity/80%-coverage representatives | 23,416,070 |
| UniRef50 no-hit representatives | 2,839,030 |
| UniRef50-hit representatives | 20,577,040 |
| Predicted-3Di no-AFDB-match representatives | 2,533,446 |
| Predicted-3Di AFDB-hit representatives | 305,584 |
| Sequence-dark multi-member families (30% clustering) | 334,810 |
| Sequence-dark singletons | 1,540,730 |

The representative partition closes exactly:
`2,839,030 + 20,577,040 = 23,416,070`. The sequence-dark rate is
`2,839,030 / 23,416,070 = 12.1243%`, and the initial structural-no-match rate is
`2,533,446 / 2,839,030 = 89.2363%` (89.2%). The latter partition also closes:
`2,533,446 + 305,584 = 2,839,030`.

The manuscript does not state the **23,416,070-representative denominator** in
the abstract, Results, or novelty-ladder table. In the abstract it calls the
2.84 million objects “proteins”; the Results and table more accurately call
them representatives. The displayed arrow from 51.9 million proteins to 2.84
million representatives therefore invites the incorrect inference that 2.84
million is the dark fraction of the 51.9-million-protein catalogue. It is
12.1% of the representatives, not 5.5% of the original protein calls.

There are two further cohort switches inside the displayed ladder:

- `334,810` is the number of multi-member families obtained by clustering the
  **full 2,839,030 sequence-dark representative set**. It is not a family count
  obtained by clustering only the 2,533,446 predicted-3Di structural no-hits,
  even though the prose and table place it after that stage with an arrow.
- `68,808` is numerically correct, but it comes from
  `/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/tier1_cluster/novel_fold_analysis/family_summary.tsv`
  (`is_novel_family == 1`; materialised as 68,808 IDs in
  `/ibex/scratch/projects/c2014/rob/metagenomics/06_esmfold_novelfold/novelfold_reps.txt`).
  That clustering starts from the separate **11,418,079 MAG-associated
  sequence-novel-protein tier**, not from the 334,810-family cohort immediately
  above it. `CHANGELOG.md` records that those 68,808 families contain 169,651
  proteins.

Thus the manuscript never directly substitutes the 753,305,533-protein
assembly catalogue for the 51,920,778 Bakta-from-MAG catalogue, but it does
silently import a downstream result from the larger assembly-catalogue/tiering
workflow into a ladder introduced as the 51.9-million MAG-protein analysis.
That is an aggregation/provenance defect, not evidence that 68,808 was
miscounted.

The separate tier analyses were also checked. Tier1 reports 11,418,079 input
proteins, 883,688 AFDB hits and 10,534,391 no-matches:
`10,534,391 / 11,418,079 = 92.26%`. Its own
`results/tier1_structural_novelty/SUMMARY.md` explicitly warns that this raw
rate is fragment-inflated (median no-hit length 56 aa; 77.8% shorter than 100
aa) and gives the defensible length-gated result as 2,339,154 no-matches at
at least 100 aa, about 74%. Tier2 independently reports
`10,113,496 / 13,088,646 = 77.27%` for non-MAG proteins already restricted to
at least 100 aa. The manuscript does not print either the fragment-inflated
92.26% or the tier2 77.27%, so it does not misstate those percentages; it does
discuss length sensitivity for the final ESMFold subset, but it never tells the
reader that its 68,808-family input came from the fragment-inflated tier1
cohort.

Finally, tier2's 156,743 novel-fold families are not the manuscript's 334,810
sequence-dark families. The manuscript does not print 156,743 or claim that
those tier2 families are generally large. The source shows only 710 tier2
families with at least 10 members and 11 with at least 50, so the known
large-family trap is avoided in the text as written.

### Check 5 — sample universe and GTDB-Tk census

**Status: completed.** The canonical committed sample metadata,
`workbench/data/metadata/sample_metadata.tsv`, contains 551 soil-sample rows:
all 551 have reads and 550 have contigs. The concrete assembly-input list,
`workbench/data/metadata/binnable_samples.txt`, contains exactly 550 samples.
The metadata have five expedition labels, 61 distinct site identifiers and 217
observed site-by-expedition combinations; they are not a complete 60-by-5
factorial.

The MAG paper's value of **547 “Soil samples”** does not describe that source.
It is instead the row count of
`results/darkmatter_ecology/env_matrix.tsv`, a downstream environmental
analysis subset. That table contains 547 samples: it excludes four samples
with contigs (`62PRr2`, `64Dr3`, `64PRr2`, and `64Sr2`) but includes `51Sr1`,
the one metadata sample without contigs. The label “Soil samples” in the scale
table is therefore wrong. The phage paper's 550 metagenomes agrees with the
assembled/binnable input set. The two manuscripts present 547 and 550 as
unqualified scale counts for the same survey and never explain the subset,
creating an internal cross-manuscript inconsistency.

The even broader canonical read tree named in `DATA_LOCATIONS.md` contains 573
paired-read samples. None of 547, 550, or 551 should be described as that
archive-wide read count. These levels are:

| Level | Count | Concrete source |
|---|---:|---|
| Canonical paired-read archive | 573 | `$EQD/soil/wgs_assemblies/processed_reads/` |
| Rows in committed soil-sample metadata | 551 | `workbench/data/metadata/sample_metadata.tsv` |
| Binnable/read+contig metagenomes | 550 | `workbench/data/metadata/binnable_samples.txt` |
| Environmental-matrix subset | 547 | `results/darkmatter_ecology/env_matrix.tsv` |

GTDB-Tk values were independently aggregated from the 14,752 data rows in
`gtdbtk.bac120.summary.tsv`, 717 rows in `gtdbtk.ar53.summary.tsv`, and the
per-MAG quality join in `novelty_summary.tsv`. All printed census values agree:

| Quantity | Recomputed | Printed |
|---|---:|---:|
| Classified MAGs | 14,752 bacterial + 717 archaeal = 15,469 | 15,469 |
| Phyla | 37 bacterial + 7 archaeal = 44 | 44 |
| Actinomycetota / Pseudomonadota | 6,788 / 3,740 MAGs | 6,788 / 3,740 |
| Novel phyla / classes | 0 / 0 | 0 / 0 |
| Candidate novel orders | 12, including 3 HQ | 12 (3 HQ) |
| Candidate novel families | 297, including 80 HQ | 297 (80 HQ) |
| Candidate novel genera | 3,677, including 996 HQ | 3,677 (996 HQ) |
| Candidate novel species | 10,815, including 2,585 HQ | 10,815 (2,585 HQ) |
| HQ novel-rank candidates | 3 + 80 + 996 + 2,585 = 3,664 | 3,664 |

Filtering `novelty_summary.tsv` at `novel_rank == d__` also exactly supports the
domain-unplaced description: 17/17 are MQ, 0 are HQ, all identifiers are
rhizosphere (`PRr`) samples, and 3 used the archaeal marker set. Completeness is
51.38–65.96% (printed 51.4–66.0%) and contamination is 2.71–9.74% (printed
2.7–9.7%).

### Check 6 — MAG-paper A1–A5 arithmetic and plotted values

**Status: completed.**

For A1, direct counting of `hydB_class.tsv`, `coxL_class.tsv`, and
`guild_core_mags.txt` gives:

- group 1l: 5,292 proteins on 2,609 MAGs;
- group 1h: 986 proteins on 893 MAGs;
- form-I CO dehydrogenase: 5,645 proteins on 3,828 MAGs;
- H2+CO guild core: 1,204 MAGs, or
  `1,204/15,469 = 7.7833%` (correctly printed as about 8%).

The 16 bars in `fig_lifeonair_by_phylum.pdf` agree with
`results/atmospheric/integrate/atmo_guild_by_phylum.tsv`: CO/H2 MAG counts are
Actinomycetota 2,175/2,198; Pseudomonadota 1,007/24; Chloroflexota 332/305;
Gemmatimonadota 132/257; Halobacteriota 64/40; Acidobacteriota 50/118;
Myxococcota 13/131; and Bacteroidota 32/16.

The abstract's “atmospheric H2/CO oxidation in about 3,400 genomes” needs an
aggregation label. In `atmo_guild_per_mag.tsv`, 3,386 is the number of MAGs
with either group 1l or 1h **H2** hydrogenase. The CO set has 3,828 MAGs; 1,204
are in the H2+CO intersection; and the H2-or-CO union is 6,010. Thus about
3,400 is numerically defensible only as the H2 set, not as either the joint
guild or the H2/CO union.

The aridity statistic was rerun from all 58 rows in
`per_site_hydrogenase.tsv` using the recorded script, seed 20260726, and 9,999
permutations. It returns Spearman `rho = -0.171588986...` and
`p = 0.1961`, exactly supporting the printed `-0.172`, `0.196`, and `n=58`.

For A2, source counts in `group_prevalence.tsv` give NER
`15,382/15,469 = 99.4376%`, BER `15,398/15,469 = 99.5410%`, and NHEJ
`12,122/15,469 = 78.3632%`, matching 99.4%, 99.5%, and 78.4%.
`compartment.tsv` records 2,241 Surface and 3,889 Deep MAGs with rounded
photolyase prevalences 22.4% and 13.1%. Every integer contingency table
compatible with those one-decimal percentages gives a two-sided Fisher p
between `6.7e-21` and `2.9e-20`, so the manuscript's `p<1e-13` is true. The
source file retains only rounded percentages, however, so the exact cell count
and exact photolyase p-value cannot be recovered from it.

For A3, the printed osmolyte-category null recomputes from the contingency
counts in `results/darkmatter_ecology/DECODE_ECOLOGY_BRIDGE.md`: 113 osmolyte
families among 9,425 structured families versus 145 among 13,342 cosmopolitan
families gives OR `1.10444` and two-sided Fisher `p=0.44587`, matching OR 1.10
and `p=0.45`. The claims ledger points instead to `OSMOLYTE.md`, which reports
a different plant-root test (`p=0.11`, `n=212`) and does not contain the
printed OR. The number is right but its recorded source is wrong.

For A4, all 12 numbers annotated in `fig_bgc_product.pdf` exactly match
`bgc_product_summary.tsv`: 12,589; 9,984; 6,878; 4,900; 4,785; 3,914; 3,505;
3,151; 2,868; 2,626; 2,042; and 2,041 MAGs, in the plotted order. The plot
states its 15,469-MAG denominator. The Results percentages recompute as
64.5420%, 25.3022%, and 16.9759%; only terpene is misrounded, as recorded in
Check 2.

For A5, beta-lactam resistance is
`3,816/15,469 = 24.6687%` (24.7%). `metal_census.tsv` contains 2,423
arsenic-positive and 1,077 mercury-positive MAGs. The named-pathogen file has
five *Escherichia coli* plus two *Acinetobacter baumannii* MAGs:
`7/15,469 = 0.04525%`, which rounds to the printed 0.05%. All 12 annotated
bars in `fig_amr_class.pdf` match `amr_class_summary.tsv`: 3,816; 2,002;
1,108; 1,086; 705; 622; 294; 127; 66; 65; 29; and 25 MAGs. The figure states
the 15,469-MAG denominator.

### Check 7 — MAG-paper B1a novel-fold sensitivity

**Status: completed.**

The canonical IBEX objects were checked directly rather than inferred from the
paper repository. `novelfold_reps.txt` has exactly 68,808 records and the
canonical `06_esmfold_novelfold` tree has 68,527 PDB model files. The 16
`plddt_shard*.tsv` files contain 63,935 data records: 63,654 `ok` plus 281
`too_long`. Thus the paper's 68,527/68,808 and 63,654 successfully scored
counts are exact.

Joining the pLDDT population to the AFDB and PDB best-TM tables is summarized
in `RQ/A3/fold_validation_summary.tsv`. Its integer numerators and
denominators recompute as follows:

| minimum pLDDT | dual-novel / retained | exact rate | printed |
|---:|---:|---:|---:|
| 0.00 | 15,248/63,654 | 23.9545% | 24.0% |
| 0.40 | 9,258/50,154 | 18.4591% | 18.5% |
| 0.50 | 5,514/33,524 | 16.4479% | 16.4% |
| 0.60 | 2,836/17,196 | 16.4922% | 16.5% |
| 0.70 | 1,164/5,709 | 20.3889% | 20.4% |
| 0.80 | 15/187 | 8.0214% | 8.0% |

The high-confidence population is `187/63,654 = 0.2938%`, correctly printed
as 0.3%. The same joined source gives dual-novel/matched median pLDDT
0.440/0.520, median lengths 182/153 aa, and the five pLDDT>=0.60 length
strata printed in the text: 26.3% (`n=9,660`), 5.4% (`n=1,839`), 2.7%
(`n=3,421`), 4.8% (`n=1,831`), and 4.9% (`n=445`). Exactly 300 dual-novel
models satisfy both pLDDT>=0.60 and length >=100 aa.

The current paper has propagated the 2026-07-29 correction: it prints 15/187
and explicitly says this is not a corrected novel-fold total. It does not
print either withdrawn value (“about 12” or “10,007 high-confidence”). The
matched known-fold false-novelty control is also correctly disclosed as not
run.

The arithmetic passes, but `claims.tsv` still sources the 68,527 count only to
“ORIX ESMFold run,” not to the canonical model directory or a manifest.

### Check 8 — MAG-paper B2/B3 selection, ordination, occupancy and abundance

**Status: completed.**

Direct aggregation of `dnds.tsv` gives 22,767 dark families with median
family dN/dS 0.391 and 21,579 controls with median 0.2859, matching the
printed 0.39 and 0.29. Of the dark families, 22,619/22,767 = 99.3499% have
dN/dS below 1, supporting the figure title's integer “99%.”

Rejoining `dnds.tsv` to site occupancy recomputed directly from
`darkfam_abundance.tsv` and `env_matrix.tsv` gives:

- endemic (<=3 sites): `n=1,809`, median 0.3324;
- EQ-wide (>=30 sites): `n=10,498`, median 0.4201;
- two-sided Mann–Whitney `p=5.31784e-69`;
- occupancy/dN/dS Spearman `rho=0.175681`, `p=3.18919e-157`.

These support 0.332 versus 0.420 and `p<1e-67`. The group sizes are not stated
in the manuscript.

The raw occupancy matrix contains 23,382 families and 551 abundance columns.
After mapping detections to the 60-site environmental scheme, 10,669/23,382 =
45.6291% occupy at least 30 sites and 535/23,382 = 2.2889% occupy exactly one
site, matching 45.6% and 2.3%. The generator's reported single-site count of
543 includes eight families with **zero mappable sites** (they are detected
only in abundance columns excluded from `env_matrix.tsv`); this does not
alter the one-decimal percentage but makes the source's integer label
technically wrong. Median occupancy is 26 sites and five trips.

Streaming all 23,382 rows of `darkfam_abundance.tsv` independently gives 307
families in at most one sample (`307/23,382 = 1.3130%`), 23,261 with peak
abundance >=10 (`99.4825%`), median peak 584.7518, and median prevalence 52
of 551 samples. These match 1.3%, 99.5%, about 585, and 52. The manuscript
does not state the 23,382-family or 551-sample denominators in B3.

The PCoA was recomputed from the raw abundance matrix and has axes 27.7216%
and 7.3196%, matching the PDF's 27.7% and 7.3%. The deterministic PERMANOVA
was rerun from the same files and exactly reproduced the figure:
compartment 2.0%, space 3.8%, climate 4.7%, climate conditional on space
3.3%, geochemistry 6.6%, geochemistry conditional on space 5.1%, and
vanadium 0.4%; all permutation `p=0.001` except vanadium `p=0.325`.

There is a silent analysis-universe mismatch inside the paired B3 figure.
The PCoA code filters prevalence across all 551 columns and retains 21,961
families before aligning to 547 environmental samples. The PERMANOVA code
aligns to 547 samples first and then retains 21,923 families. Neither
denominator appears in the figure or caption.

### Check 9 — MAG-paper B5 decoding

**Status: completed, with a provenance and cohort defect.**

The manuscript's cited source, `results/darkmatter_ecology/family_ec.tsv`,
has 23,382 rows and no review-verdict field, so **10 cannot be derived from
that file**. The concrete source is instead the workbench file
`results/darkmatter_function_assignment/scores_v3/fanout/top100_reviewed.tsv`.
Counting its verdict column gives exactly 10 `well_supported_EC`, 39
`plausible_EC`, and 51 `pathway_associated` calls. Thus 10 is right but the
claims ledger points to the wrong file.

The “3 high-priority structural targets” are also exact, but arise from a
different 100-row shortlist. Applying mean pLDDT >=0.60 and top-pocket
probability >=0.4 to
`dossier/reallydark_top100_secondary_pockets.tsv` gives exactly three:
periplasmic nitrite reductase (0.75/0.870), uroporphyrinogen-III synthase
(0.76/0.598), and a malate:quinone oxidoreductase/PEP-carboxylase call
(0.80/0.524). The denominator collapses to 73 distinct folded proteins, but
the manuscript does not state 3/73 or the thresholds.

The printed calibration range does **not** describe the current ten-call
cohort. In `top100_reviewed_calibrated.tsv`, the ten well-supported calls
have calibrated `P(correct)` 0.219–0.432 (median 0.262), i.e. about
0.22–0.43, not 0.34–0.48. The latter range belongs to older/different
shortlists. This is a numerical and aggregation error, not merely a vague
source citation.

### Check 10 — MAG-paper C1–C3 reference bias, cross-feeding and specificity

**Status: completed.**

For C1, direct aggregation of `per_sample_concordance.tsv` gives 436 rows,
median NSTI 0.1591, PICRUSt2 KO count 10,543, MAG KO count 4,482.5, and
presence recall 0.9327. The rounded 0.16, 10,543, 4,482 and 93% are
numerically reproducible (the MAG median is exactly half-integral). The
10,543 value is constant in all 436 rows, not merely a typical median.

From `dark_fraction.tsv`, the no-KO result is
`31,413,266/55,380,378 = 56.7227%`. It is therefore 56.7% of CheckM2
`Total_Coding_Sequences`, not of the 51,920,778 Bakta protein catalogue.
The manuscript calls both products “MAG proteins” without disclosing this
change of denominator. The per-MAG median is 56.05% (the figure prints 56%).
Figure group sizes are exactly 651 known species, 10,815 novel species,
3,677 novel genera and 326 novel-family-or-higher MAGs.

The workbench rebuild script was rerun read-only from the original PICRUSt2
and MAG matrices. It exactly reproduces the seven figure ratios:
2.39619, 1.97170, 1.78596, 1.48413, 1.11480, 0.573007 and 0.0165484 at
thresholds 0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4 and 1e-3 over 440 samples.
However, the manuscript's description of a **shared** threshold is false.
Those ratios divide thresholded PICRUSt2 richness by an unthresholded,
constant MAG-presence denominator (mean 4,346.4818 at every threshold).
When the same relative-abundance threshold is actually applied to both
PICRUSt2 and the MAG relative-abundance matrix, the ratio at 1e-5 is
`4,845.4523/331.125 = 14.6333`, not 1.11, and it remains 3.8846 at 1e-3
rather than reversing. The 1.11 and reversal are correct only for the
hybrid numerator/denominator encoded in `threshold_sweep.tsv`.

Filtering `marquee_functions.tsv` at MAG prevalence >=0.5 leaves exactly ten
functions. Their correlations range from -0.123 to 0.417 with median 0.086,
matching -0.12 to 0.42 and 0.09.

For C2, the 16 SMETANA global files give MRO range
0.549746–0.777843 and median 0.688896, supporting “about 0.7.” The FBA source
contains 316 MAGs from 16 pilot communities: 27 guild and 289 non-guild.
Donor medians are 232 versus 243, and the recorded one-sided Mann–Whitney
calculation gives `p=0.693779` (0.69). Donor-fraction medians are 0.503
versus 0.495 (`p=0.387733`). Completeness correlations recompute as 0.0344
for receiver and 0.0605 for donor. The manuscript states the 16-community
pilot size but not the 27/289 denominator for its `p=0.69` claim.

For C3, direct set operations on the canonical hit files give:

| comparison | hit | absent | absent fraction |
|---|---:|---:|---:|
| MGnify90 | 14,830 | 8,552 | 36.5751% -> 36.6% |
| GEM | 12,307 | 11,075 | 47.3655% -> 47.4% |
| union of five other deserts | 15,913 | 7,469 | 31.9434% -> 31.9% |
| hit in any of those three sources | 18,371 | 5,011 | 21.4310% -> 21.4% absent everywhere |

All use the 23,382-orthogroup ecology subset, not the 334,810 family
clustering. The abstract attaches 21.4% to the 334,810 families and is
therefore at the wrong aggregation level. The same 23,382 subset is also
mislabelled in the phage paper; see the phage checks below.

Remote job logs give 122,978,191 GEM proteins and five desert catalogue
sizes that sum to 147,531,178. The manuscript's 123M is correctly rounded,
but 147M is truncated; nearest-million rounding is 148M, which is what the
generated `MULTICATALOG.md` table prints. The hard-coded 718M MGnify depth
has no local count/manifest from which it can be recomputed.

### Check 11 — phage repository verifier and coverage boundary

I ran the phage repository's
`reproducibility/reproduce_numbers.py` without modification. It reports
**87/87 PASS**. I also inspected what those checks actually compute rather
than treating the exit code as an audit result.

The verifier independently recomputes, from committed `results/*.tsv` files,
the vOTU catalogue and lifestyle partition; IMG/VR novelty; CheckV tiers;
DefenseFinder totals; defense-island fractions; retrons; defense and
anti-defense carriage; bulk load correlations; three endemicity values;
rarefaction; SNTC detection; the top operon candidate; pangenome flux;
spacer-network statistics; pN/pS directions; AFDB novelty; AMGs; TerL tips;
the 30-family PDB100 screen; and vConTACT2 cluster novelty. Those checks are
useful and all of their printed/recomputed comparisons pass.

They do **not** establish complete provenance for the manuscripts. Six
quantities are merely printed as “DOCUMENTED constants,” not recomputed:
550 metagenomes, 15,469 MAGs, 185 isolates/118 cultures, 312,567 pOTUs,
23,382 proteins, and 455 candidate families. Several other manuscript
claims are not checked at all, including the flagged-control count, plasmid
length distribution, host-phylum counts, most defense-system table rows,
the network null denominators, candidate-table rows, Borg count, exact
matched-enrichment statistics, and all external positioning-table values.
The script also calls the 23,382 records “structurally-dark proteins”
without testing whether that cohort belongs to this paper. It does not:
that issue is audited below.

### Check 12 — phage scale, saturation, controls and CheckV

The committed phage tables recompute as follows:

| claim | numerator / source rows | denominator | recomputed | printed |
|---|---:|---:|---:|---:|
| vOTU catalogue | 18,604 rows in `votu_catalog.tsv` | — | 18,604 | 18,604 |
| ≥Medium CheckV | 962 + 955 + 3,037 | 18,604 | 4,954; 26.6287% | 4,954; 26.6% |
| MIUViG ≥90% | 1,917 | 18,604 | 1,917 | 1,917 |
| complete genomes | 1,474 | 18,604 | 1,474 | 1,474 |
| CheckV median | per-vOTU estimates | 17,646 with an estimate | 27.7% | 27.7% |
| vOTU Chao2 coverage | 18,604 | 18,620.8 | 99.9098% | 99.91% / 99.9% |
| pOTU Chao2 coverage | 312,567 | 312,589.3 | 99.9929% | 99.99% |
| SNTC signal | 3,182 | 331,171 total OTUs | 0.96083% | 0.96% |
| control ≥ mean real | 3,722 rows | 331,171 total OTUs | 1.12389% | 1.12% |

The CheckV quality tiers also partition 18,604 exactly. The contamination
denominator is itself 18,604 vOTUs + 312,567 pOTUs = 331,171 OTUs, and the
file contains exactly 3,722 data rows after its comments/header.

The 550 metagenomes and 15,469 MAGs agree with the canonical metadata and
taxonomy sources already checked above. The pOTU count agrees with the
canonical path/count in `DATA_LOCATIONS.md`; it is not recomputed from a
committed pOTU catalogue in this repository. No source file for the printed
**185 isolates / 118 enrichment cultures** is present in either paper
repository or the workbench metadata. The phage verifier labels both numbers
“DOCUMENTED” and cites only “sequencing manifest,” without naming or supplying
one. This is a **(b)** provenance defect.

The outline additionally says “7 show strong control dominance” and that
excluding flagged OTUs changes no result. Neither the per-OTU table, its
header, nor the named upstream script (which is absent at the cited workbench
path) defines/recomputes the seven or contains the exclusion sensitivity
results. These are also **(b)**.

### Check 13 — anti-defense arithmetic exposes a wrong contig count

`results/antidefense_systems.tsv` contains **483** system rows, so that total
is correct. The number of carrier vOTUs is not 480:

- using `sys_beg`, remove only the terminal protein ordinal (`_<integer>`);
- this yields **455 distinct contig IDs**;
- all 455 match IDs in `votu_catalog.tsv`;
- they comprise 432 contigs with one system, 19 with two, three with three,
  and one with four: 432 + 38 + 9 + 4 = 483 systems.

Thus the abstract/result/table/outline claim **“483 systems on 480 phage
contigs” is (a) WRONG; the source gives 483 systems on 455 vOTUs**.
The repository verifier checks only the 483 rows and never checks 480.
The upstream merge script attempts to infer contigs from `sys_id` with
`sed 's/_[^_]*_[0-9]*$//'`; anti-defense subtype names themselves contain
underscores, so `sys_id` is not a safe field for this count. Even that
command gives 464 on the committed table, not 480.

There is a second, visible arithmetic inconsistency in the manuscript's
anti-defense table. Its eight displayed types sum to
110+103+101+68+51+20+17+8 = **478**, not the stated total 483. The source's
five omitted systems are Anti_TA (2), Anti_RecBCD (1), Anti_Gabija (1), and
Other (1). Because the table does not label its body as a selected subset or
include an “other” row, this is **(d)** as printed.

### Check 14 — phage identity, host, gene-sharing, TerL and AMG claims

The vOTU lifestyle and species-novelty claims recompute exactly:

- 7,357 virulent + 5,713 temperate + 5,534 unclassified = **18,604**;
- 18,251 IMG/VR no-hits / 18,604 = **98.1026% -> 98.1%**, with 353
  hits completing the partition;
- 2,494/18,604 = **13.4057%** have a named genus, supporting the
  outline's 13.4%;
- 13,583/18,604 = **73.0112%** have a non-empty lineage, supporting
  the outline's 73.0%;
- host calls partition into 4,135 assigned + 14,469 unassigned = 18,604;
  4,135/18,604 = **22.2264%**. Every phylum count in the host figure,
  including 2,114 Actinomycetota and 1,347 Pseudomonadota, matches
  `votu_host_phylum_summary.tsv`.

The external comparator is also transcribed correctly from Ma et al.: 1,824
soil metagenomes, 80,750 partial viral genomes, and 96.7% taxonomically
unassigned. The cited Atacama PDF reports 84 viral populations; therefore
18,604/84 = **221.5-fold**, consistent with the non-exact “about two orders”
positioning claim.

The TerL tree file/summary gives 3,009 EQ tips + 4,570 reference tips =
**7,579**, all correct. The vConTACT2 overview independently gives:

- 22,476 inputs = 17,942 EQ + 4,534 references;
- 2,378 EQ-only + 138 shared = **2,516 EQ-containing VCs**;
- 2,378/2,516 = **94.514% -> 94.5%**;
- 8,954 EQ sequences are in numbered VCs, of which 7,970 are in EQ-only VCs.

The main manuscript's genus-scale sentence states its immediate
2,516-cluster denominator, but it does not say that only 8,954/18,604 =
**48.13%** of the full vOTU catalogue is represented in those numbered
clusters, nor that only 17,942 entered the network. This is a **(c)**
inclusion/aggregation omission when the result is generalized to the
“desert virome.”

The outline's DRAM-v claim contains a genuine counting bug. In
`dramv_amg_summary.tsv`:

- 1,965 rows have `potential_amg=True` (correct);
- those rows belong to **1,537 distinct scaffolds/vOTUs**, not 2,928;
- 191 of the 1,965 are CAZy (correct);
- category counts 1,040 information systems, 586 amino acid, 56
  corrin-ring, 43 central-carbon, and the smaller printed categories
  also match the filtered rows.

The verifier obtains 2,928 by counting distinct scaffold IDs over **all
5,958 summary rows**, without applying its own `potential_amg=True` filter.
Consequently, the outline and claims ledger statement **“1,965 candidate
AMGs across 2,928 vOTUs” is (a) WRONG; the source supports 1,965 across
1,537 vOTUs**. The current main manuscript leaves AMG numbers pending and
does not print this error.

### Check 15 — phage plasmidome and defensome arithmetic

The EQ pOTU arithmetic is internally consistent:

- 33 conjugative + 466 mobilizable + 312,068 without a detected marker =
  **312,567** pOTUs;
- 312,068/312,567 = **99.8404% -> 99.84%**;
- the catalogue audit records 584,654 input contigs, 312,567
  representatives, median length 4,635 bp (**4.6 kb**), 57,751 ≥10 kb,
  and total length 2.50 Gbp;
- the recorded maximum is **643,799 bp = 643.799 kb**, so the unqualified
  main-text “max 643 kb” is truncated; nearest-whole-kb reporting would be
  **644 kb** (**(a)**);
- phage carriage is 361/18,604 = **1.9404% -> 1.9%**;
- plasmid carriage is 8,221/312,567 = **2.6302% -> 2.6%**;
- ≥10-kb plasmid carriage is 4,199/57,751 = **7.2709% -> 7.3%**;
- the system counts 419, 9,780, and 5,608 all match
  `defense_on_mobilome.tsv` after the stated element/length filters.

The outline and RQ documentation, however, label 9,780/419 = **23.3×**
and 5,608/419 = **13.4×** as enrichment, “per-bp,” or a
plasmid-versus-chromosome contrast. The verifier simply divides raw plasmid
system counts by raw phage system counts; it has neither chromosome bases
nor plasmid bases in its calculation. From the printed element rates, the
actual plasmid/phage *carriage-rate* ratios are only **1.36×** (all pOTUs)
and **3.75×** (≥10 kb). The current main manuscript wisely does not print
23.3×/13.4×, but the outline's robustness/enrichment claim is **(a)
WRONG as labelled** and its supporting provenance document describes an
analysis that the source/verifier did not perform.

The external GSPR numbers themselves are present in the cited PDF: 6,860
soil data sets, 98,728 collected sequences, and 57,265 pOTUs. The claimed
reason for the count difference is not. GSPR filters metagenomic plasmid
contigs at ≥5 kb, combines them with isolate plasmids, dereplicates at 100%
nucleotide identity/alignment fraction to 98,281, then forms 57,265 groups
with a BLAST/Leiden procedure using a 70% alignment-fraction floor. It is
neither a complete-plasmid-only catalogue nor clustered by the same 95/85
rule. Therefore:

- main manuscript: “all plasmid contigs ≥2.5 kb vs complete plasmids” is
  **(a) unsupported/wrong**;
- outline: “pOTU-vs-pOTU comparison, same 95/85 rule” is **(a) WRONG**.

For the defensome, 264 rows in `defense_system_summary.tsv` sum to
**25,096 systems**; 7,959 carrier MAG rows give
7,959/15,469 = **51.4513% -> 51.5%**. All 12 system/MAG pairs in the
manuscript table match the source. The island file has 24,776 positioned
systems, of which the printed 13.1% island / 86.9% solitary split and the
Gabija 32.4%, pAgo 31.5%, CBASS 26.8%, and Cas 6.9% values recompute at
the stated type/positioned-system level. The 357 DefenseFinder retron rows
and 305 reconciled carrier MAG rows support “about 357 across about 300.”
No isolate/culture result or manifest is supplied from which the additional
claim “absent from isolates/cultures” can be recomputed (**(b)**).

### Check 16 — phage network, load, pangenome and selection

`spacer_host_virus_network.tsv` has 885 unique host-vOTU pairs, 361 unique
host MAGs and 637 unique vOTUs. Its `n_spacers` column sums to **1,632**,
so the outline's 1,632 raw links and the manuscript's 885 collapsed edges
are both right and represent different aggregation levels. The outline does
not explain that transition where the two values appear (**(c)**). No
committed source contains the separate outline-only **3,281 iPHoP
host-linked** count (**(b)**).

All saved network statistics match the manuscript:

| metric | source / recomputation | printed |
|---|---:|---:|
| modularity | 0.965 vs null 0.9132; z=9.98 | 0.965 vs 0.913; +9.98 |
| empirical modularity p | 1/1,001 = 0.0009990 | ≈0.001 |
| NODF | 0.733 vs null 1.065; z=-15.92 | z=-15.9 |
| co-host phylum | 0.862 vs 0.212; p=1/1,001 | 86.2% vs 21.2%; ≈0.001 |
| co-host site | 0.253 vs 0.027; ratio 9.370× | 25.3% vs 2.7%; ≈10× |
| iPHoP phylum agreement | 26/27 = 96.296% | 26/27 (96%) |

The two co-host tests are over **154 multi-host vOTUs**, a denominator saved
as `cohost_*_nVOTU=154` but absent from the main table/prose; this is
**(c)**. The source also records `NODF_p=1.0`, apparently the wrong
one-sided tail for a lower-than-null statistic. The z=-15.92 does support
an extreme anti-nested departure, but the repository does not preserve a
correct lower-tail empirical p-value; the manuscript appropriately prints
no NODF p-value.

The unblocked load table has exactly 512 samples and recomputes to Spearman
rho=0.186 (normal-approximation p=1.9e-5) for viral density and rho=0.030
(p=0.49) for richness. The current manuscript/outline instead make the
more cautious claim that the association does not survive site × trip ×
compartment blocking. No blocked-permutation output file is committed, and
the outline explicitly says its exact p-value is withheld pending source
re-verification. The 512 denominator is right, but the negative blocked
result is **(b) unfollowable** from the supplied sources.

For defense repertoire turnover, the 135 genus rows have median zero and
mean 0.296 core systems, mean accessory fraction 0.955319 (**95.5%**), and
mean pairwise Jaccard distance 0.679052 (**0.68**). The genome-wide
“about 100% accessory” comparison comes from only **12 genera**, with mean
0.9985 and median 1.0, not from all 135. The outline discloses 12; the main
manuscript does not, making its sentence a **(c)** aggregation omission.
The outline-only 132-genera within-genus analysis also checks: 64% positive
and median rho=0.077.

All outline H1/H2 selection numbers match `microdiv_pnps_summary.tsv`:
defense pN/pS 0.314199 (`n=149`) versus control 0.231662
(`p=0.010759`); defense FST 0.257559 (`n=78`) versus 0.183448
(`p=0.013131`); and anti-defense pN/pS 0.140750 (`n=177`,
`p=2.8097e-7`). The abandoned dN/dS diagnostic also matches its source:
the control has median omega 0.6131 and 88.0% BUSTED-significant families.

### Check 17 — phage structural novelty, candidate chain, endemicity and Borg

The independent mobilome/defense structural search is internally exact:
316,650 AFDB hits + 296,510 no-hits = **613,160 clust90
representatives**, and 296,510/613,160 = **48.3560% -> 48.4%**.
This is a representative-level count and the denominator/aggregation is
stated where it appears.

All seven displayed candidate-table rows match
`candidate_defense_families.tsv` after the manuscript's rounding:
Armada 6 MAGs/OR 1,645/`p=2.3e-7`; CBASS and Gabija
5/1,391/`4.8e-7`; RM_Type_IIG and McrBC 4/1,138/`1.2e-6`;
Cas 5/464/`1.3e-11`; and RloC 5/126/`5.1e-16`.
However, the committed file contains only the top **200**, not all 455
families. No supplied file or extant named workbench script contains the
upstream **541,569 unknown proteins -> 364,544 clusters -> 455 recurring
families** totals. The verifier explicitly treats 455 as a documented
constant while checking only that the shortlist has 200 rows. All three
chain counts therefore have **(b) unfollowable provenance**.

The structural shortlist itself is auditable. It has 30 rows: **29**
`structurally_novel` and **1** `known_fold`. For the nominated pair,
the source gives:

- OR 1,137.9806, Fisher `p=1.17738e-6`, and four island MAGs;
- pLDDT 0.72 and 0.71, lengths 170 and 62 aa;
- no PDB100 hit for either;
- five predicted pockets on the larger protein, best score 8.21.

Thus 29/30, the rounded OR/p, four-MAG recurrence, 0.72/0.71, two-gene
description, and outline pocket score 8.2 all match. The outline's excluded
Boltz-2 values (0.126/0.128; controls up to 0.371) are present only in a
prose summary, not a raw result file (**(b)**).

The endemicity result is attached to the wrong object. Direct canonical hit
set arithmetic, already shown in Check 10, verifies every percentage in the
phage table: Atacama 89.1%, autotroph desert 62.8%, Antarctic 55.9%,
Mojave 51.9%, arid 41.7%, combined deserts 31.9%, and GEM 47.4%.
But the 23,382 queries are 50%-identity orthogroups derived from Prodigal
calls over 5,438 dRep **soil MAGs**. They are neither mobilome proteins nor
a subset of the phage paper's 296,510 AFDB no-hits. Therefore the caption
**“23,382 structurally-dark mobilome proteins” and the accompanying
endemicity claims are (a) WRONG at the object/provenance level**, even
though the percentages correctly describe the MAG-paper cohort. This is
also **(d)** across manuscripts: the same seven values appear under two
different labels. The 147.5M combined-desert denominator and rounded 123M
GEM denominator are otherwise correct here.

Finally, `borg_borg_candidates.tsv` contains **9** rows in the specified
600 kb–1.2 Mb / rounded-GC 30–40% window (length 600,387–956,840 bp,
GC 32.7–40.0%). The GTDB census contains zero Methanoperedens MAGs.
The main manuscript's “0; 9” is correct. The outline and `results/README.md`
say **10**, so the outline is **(a) WRONG** and the pair is **(d)**
internally inconsistent.

### Check 18 — phage figures and external positioning table

Text extraction from the three included figure PDFs confirms that their
numeric labels are synchronized with their TSVs: the lifestyle figure prints
18,604 and 7,357/5,713/5,534; the host figure prints 4,135 assigned,
14,469 unassigned and the ten displayed phylum counts; and the defensome
figure prints the 15,469 denominator and all 15 displayed carrier-MAG counts.

The literature numbers in the manuscript positioning table also match their
specifically cited primary articles:

| study | paper prints | cited source |
|---|---:|---:|
| Hwang | 11 samples; 84 viral populations | 11; 84 |
| Ma | 1,824 metagenomes; 80,750 partial viral genomes; 96.7% unassigned | 1,824; 80,750; 96.7% |
| Li | 607 metagenomes; 29k MAGs | 607; 29,031 |
| Han | 191 metagenomes; 13,336 viral genomes | 191; 13,336 |
| GSPR | 6,860 data sets; 57,265 pOTUs | 6,860; 57,265 |
| Meaden | 1,075 metagenomes | 1,075 |

The “about two orders” comparison is conservative:
18,604/84 = **221.48×** (2.345 orders of magnitude). These background
values have article-level provenance rather than repository-only pointers;
no category (a)–(d) defect arises from their transcription. The GSPR
*interpretation* remains defective as described in Check 15.

### Check 19 — final canonical-source reconciliation

The last pass returned to the canonical remote products rather than relying
on their copied summaries.

For the plasmid catalogue,
`/ibex/scratch/projects/c2014/rob/metagenomics/06_mobilome/plasmidome/potu_reps.fna`
contains exactly **312,567** representatives. Direct sequence-length
aggregation gives minimum 2,500 bp, median 4,635 bp, 90th percentile
15,387.4 bp, 95th percentile 23,812.4 bp, 99th percentile 55,670.36 bp,
maximum **643,799 bp**, 57,751 sequences at least 10 kb, and
2,500,595,314 total bp. This independently confirms all of Check 15's
catalogue figures and confirms that “643 kb” is truncation rather than
rounding.

The candidate-family chain is recoverable from the canonical remote
directory
`/ibex/scratch/projects/c2014/rob/metagenomics/06_mobilome/operon_candidates/`:

- `unknowns.faa` and `unknown_meta.tsv` each contain **541,569** unknown
  proteins;
- `clu_rep_seq.fasta` contains **364,544** representatives, and
  `clu_cluster.tsv` maps all 541,569 proteins to exactly 364,544 families;
- joining `clu_cluster.tsv` to `unknown_meta.tsv` yields **455** families
  with at least two *member proteins* annotated as defence-island
  neighbours, but only **412** families in at least two distinct
  island-bearing MAGs.

Thus the first two printed chain counts are correct. The phrase “455
families recurring in at least two defence islands” is not: it silently
uses member-protein multiplicity as island multiplicity. Under the
manuscript's own `n_island_MAGs` aggregation, the recomputed count is
**412**, so this is **(a)** as worded and **(c)** as an unstated aggregation
switch. The committed `candidate_defense_families.tsv` contains only the
top 200 and the manuscripts do not identify the canonical upstream files,
so repository-level provenance remains **(b)** even though the audit could
locate the data through `DATA_LOCATIONS.md`.

Finally, the global-endemism selection direction was checked against
`results/darkmatter_ecology/GLOBAL_ENDEMISM.md`. MGnify-absent families
have median dN/dS **0.393**, versus **0.390** for MGnify-shared families
(`p=0.00709`). A higher dN/dS is marginally *weaker*, not stronger,
purifying selection. The MAG Results claim that globally endemic families
are under “stronger purifying selection” is therefore **(a) WRONG**. It
appears to conflate this global comparison with the separate within-EQ
occupancy result (site-restricted 0.332 versus EQ-wide 0.420), which does
show stronger constraint in locally restricted families.

### Check 20 — remaining thresholds, qualifiers and supplement-only numbers

Several quantitative statements do not look like headline counts but still
require an aggregation definition.

In MAG A1, form-I CO dehydrogenase is called the “dominant marker form.”
`coxL_form_summary.tsv` contains 13,675 CoxL-family proteins: 5,645 form I,
27 form II/other and **8,003 ambiguous**. Form I is dominant only among
*assigned* forms, not among all CoxL-family calls. The assignment denominator
is omitted, so the statement is **(c)** as written.

In MAG B1a, 68,527 PDB models were built but only 63,935 have a pLDDT record:
63,654 usable and 281 rejected as too long. The paper explains the 281 but not
the **4,592 models (6.70%)** absent between the model-directory and pLDDT
manifest. The downstream rates correctly use 63,654; the unreported attrition
is **(c)**.

In MAG B4, direct aggregation of `positive_families.tsv` gives **148/22,767
= 0.6501%** families with mean pairwise NG86 dN/dS greater than one. The
paper describes this tail parenthetically as “per-site FUBAR,” but
`STRUCT_FUNCTION.md` records HyPhy per-site output for only 65 families and
says only 12 have at least eight members. Thus the `>1` tail exists, but it is
a pairwise-NG86 screen, not a per-site-FUBAR cohort; the method and aggregation
are silently mixed (**(c)**). The lectin/GlcNAc vignette traces to
`boltz_affinity.tsv`; no numeric Boltz score is printed in the main paper.

For MAG C1, the printed 10,543 “median” is arithmetically true but is actually
constant in all 436 rows of `per_sample_concordance.tsv`; only the MAG-side KO
count varies (median 4,482.5). This is not a rounding error, but the
aggregation description conceals that PICRUSt2 assigns the same KO richness
to every sample (**(c)**).

For the phage controls, 3,182/331,171 = 0.96083% and
3,722/331,171 = 1.12389%, where 331,171 is 18,604 vOTUs plus 312,567 pOTUs.
That denominator is not shown with the two percentages (**(c)**). Likewise,
the CheckV median 27.7% is over the **17,646** vOTUs with an estimated
completeness value, not all 18,604; the table does not state that denominator
(**(c)**).

Finally, the four domain-unplaced-bin supplement was checked against
`RQ/BORG/candidates.tsv`. Total-length range 710,481–1,027,754 bp
(0.71–1.03 Mb), contigs 48/75/221/157, N50 4,679–24,157 bp, hypothetical
fractions 88.70–98.95%, tRNAs 6–73 and GC 29.1–50.1% all reproduce. The CDS
counts are 908, 1,046, 1,425 and **901**. The printed range
**908–1,425** therefore excludes the true minimum and is **(a) WRONG**; it
should be 901–1,425. These four fragmented bins are distinct from the
nine-contig 600-kb–1.2-Mb screen in `borg_borg_candidates.tsv`.

## Findings ranked by severity

### (a) Number is wrong

1. **MAG sampling scale:** “547 soil samples” is the environmental-matrix
   subset, not the survey total. The concrete universes are 551 metadata rows
   and 550 assembled/binnable samples. “60 sites × 5 expeditions” also does not
   follow from the named metadata, which contain 61 site IDs, five expedition
   labels and 217 observed site×expedition cells.

2. **MAG A4 terpene prevalence:** printed **64.6%**; recomputed
   **9,984/15,469 = 64.541987% -> 64.5%**. The verifier's loose tolerance
   masks this one-decimal rounding error.

3. **MAG abstract, 21.4% of 334,810 families:** the **21.4%** is
   5,011/23,382 ecology orthogroups, not a fraction of the 334,810
   sequence-dark-family clustering. The percentage is attached to the wrong
   object.

4. **MAG B5 calibration:** printed `P(correct) ≈ 0.34–0.48`; the current ten
   `well_supported_EC` calls span **0.219–0.432** (median 0.262).

5. **MAG C1 shared-threshold claim:** the printed 1.11× at `10^-5` and
   reversal above that threshold use a thresholded PICRUSt2 numerator divided
   by an unthresholded MAG denominator. A genuinely shared threshold gives
   **14.633×** at `10^-5` and **3.885×** at `10^-3`; it does not reverse.

6. **MAG C3 selection direction:** MGnify-absent families have median dN/dS
   **0.393**, versus **0.390** for shared families (`p=0.00709`). That is
   marginally weaker, not “stronger,” purifying selection.

7. **MAG C3 desert-catalogue depth:** the five source catalogues total
   **147,531,178** proteins. Nearest-million reporting is **148M**, not the
   unqualified 147M printed in the table.

8. **Phage anti-defence carriers:** printed **483 systems on 480 contigs**;
   `antidefense_systems.tsv` has 483 system rows on **455 distinct vOTUs**.

9. **Phage outline AMG carriers:** printed **1,965 candidate AMGs across
   2,928 vOTUs**; after applying `potential_amg=True`, the 1,965 rows occupy
   **1,537** distinct vOTUs. The 2,928 count was taken over all 5,958 DRAM-v
   rows before filtering.

10. **Phage plasmid maximum:** the canonical pOTU FASTA has maximum
    **643,799 bp = 643.799 kb**; nearest-whole-kilobase reporting is 644 kb,
    not the unqualified 643 kb in the main text and outline.

11. **Phage outline 23.3×/13.4× “enrichment”:** these are merely
    `9,780/419` and `5,608/419`, ratios of raw system counts. The calculation
    contains neither base-pair exposure nor chromosomes. Carrier-rate ratios
    are **1.36×** (all pOTUs versus vOTUs) and **3.75×** (>=10-kb pOTUs
    versus vOTUs).

12. **Phage GSPR comparison rationale:** GSPR did not analyse only complete
    plasmids and did not use the same 95% ANI/85% AF clustering. Its 57,265
    groups derive from >=5-kb metagenomic contigs plus isolates, 100%
    dereplication and a BLAST/Leiden grouping. Both the main-paper explanation
    and the outline's “same 95/85 rule” are wrong.

13. **Phage operon recurrence:** **455** is the number of families having at
    least two island-neighbour *member proteins*. Only **412** occur in at
    least two distinct island-bearing MAGs. “455 families recurring in >=2
    defence islands” uses the wrong aggregation.

14. **Phage endemicity table:** all percentages belong to the MAG paper's
    23,382 soil-MAG orthogroups, not to mobilome proteins or the 296,510
    AFDB-no-hit mobilome representatives. The table's quantitative object is
    wrong even though its arithmetic describes the other paper's cohort.

15. **Phage Borg outline:** the 600-kb–1.2-Mb/GC-screen file has **9** rows,
    not 10.

16. **Phage four-bin supplement:** printed CDS range **908–1,425**; the four
    source counts are 908, 1,046, 1,425 and **901**, so the correct range is
    **901–1,425**.

### (b) Number is right but provenance is unfollowable

1. **MAG headline chain:** seven `claims.tsv` entries cite only `workbench`,
   the nonexistent `workbench/results/darkmatter`, “novel-fold clustering,”
   “is_novel_family,” or “ORIX ESMFold run.” The numbers 51,920,778;
   2,839,030; 2,533,446; 334,810; 68,808; and 68,527 are recoverable only
   after locating concrete canonical files independently.

2. **MAG sample scale:** `claims.tsv` says only “sampling design metadata”;
   it names neither the 551-row metadata, 550-row binnable list, 547-row
   environmental matrix, nor a rule producing 60 sites.

3. **MAG A2 photolyase test:** 22.4%/13.1% and `p<10^-13` are supported, but
   the committed table keeps rounded rates rather than exact cells or the
   original test result.

4. **MAG A3 osmolyte null:** OR 1.10 and `p=0.45` are exact, but
   `claims.tsv` points to `OSMOLYTE.md`, which contains a different
   `p=0.11, n=212` test. The correct source is
   `DECODE_ECOLOGY_BRIDGE.md`.

5. **MAG B5 decode count:** the printed ten calls are in
   `top100_reviewed.tsv`; the ledger instead cites `family_ec.tsv`, whose
   23,382 rows have no review-verdict column.

6. **MAG MGnify depth:** 718M is hard-coded in prose with no local manifest
   or count file. The other catalogue depths were recomputed.

7. **Phage isolates/cultures:** no named manifest supports **185 isolates /
   118 enrichment cultures**.

8. **Phage controls:** no file defines the outline-only “7 strong
   control-dominance” OTUs or preserves the claimed exclusion sensitivity.

9. **Phage mobility matching:** the carrier data support **3.7×** and
   **1.9×**, both `p=0.0005`, and the paper correctly warns that conjugative
   `n=33` is small. Neither manuscript file identifies a saved matched-control
   resampling output; the outline still says re-verification is pending.

10. **Phage retrons:** 357 systems and 305 MAGs are reproducible, but no
    isolate/culture result supports “absent from isolates/cultures.”

11. **Phage host links:** no committed file supports the outline-only
    **3,281 iPHoP host-linked** count.

12. **Phage blocked load result:** the unblocked 512-sample correlations are
    reproducible, but no blocked-permutation output supports the retained
    qualitative negative result.

13. **Phage structural pair:** the outline's Boltz-2 values
    0.126/0.128 and controls up to 0.371 occur only in narrative, not a raw
    result file.

14. **Phage operon upstream chain:** the canonical unknown-protein and
    clustering files exist remotely, but the paper repository supplies only
    a top-200 table and points to a non-resolving workbench script/log rather
    than those files.

### (c) Number is right but denominator or aggregation is unstated

1. **MAG dark fraction:** 2,839,030 is **12.1243% of 23,416,070
   representatives**, not of the 51.9M Bakta catalogue. The denominator is
   absent and the abstract calls representatives “proteins.”

2. **MAG novelty ladder:** the arrows silently cross proteins,
   representatives and families, then cross catalogues. The 334,810 cluster
   all 2,839,030 sequence-dark reps; 68,808 comes from the separate
   11,418,079-protein tier1 workflow rooted in the 753,305,533 assembly
   catalogue. Tier1's source warns that its 92.26% no-hit rate is
   fragment-inflated; the paper does not propagate that input-cohort caveat.

3. **MAG atmospheric shorthand:** about 3,400 is the H2 union (3,386), while
   H2+CO is 1,204 and H2-or-CO is 6,010. “H2/CO in ~3,400 genomes” does not
   specify the operation.

4. **MAG CoxL dominance:** 5,645 form-I calls are dominant only among
   assigned forms; 8,003 of 13,675 CoxL-family proteins are ambiguous.

5. **MAG B1a attrition:** 68,527 structures become 63,935 pLDDT records
   (63,654 usable + 281 too long); 4,592 models disappear without a stated
   reason.

6. **MAG B2/B3 ecology:** group sizes 1,809/10,498 and the
   23,382-family/551-sample abundance denominators are omitted. The paired
   PCoA/PERMANOVA panels also use 21,961 versus 21,923 filtered families after
   applying alignment and prevalence filters in opposite orders.

7. **MAG B4 selection:** the 148-family `dN/dS>1` tail is a mean-pairwise
   NG86 screen; per-site HyPhy exists for 65 families and is adequately
   populated for about 12. The parenthetical “per-site FUBAR” hides that
   aggregation switch.

8. **MAG B5 structural targets:** three targets are selected from 73 folded
   proteins using pLDDT>=0.60 and pocket probability>=0.4; neither 3/73 nor
   the gates are stated.

9. **MAG C1 denominator:** 56.7% is
   31,413,266/55,380,378 CheckM2 total coding sequences, not a fraction of
   the 51,920,778 Bakta proteins. Separately, “median 10,543 predicted KOs”
   conceals that 10,543 is constant in all 436 samples.

10. **MAG C2 donor test:** the `p=0.69` comparison contains 27 guild and 289
    non-guild MAGs; the text states 16 communities but not those group sizes.

11. **Phage control/quality denominators:** 0.96% and 1.12% use 331,171
    combined vOTU+pOTU rows; CheckV median 27.7% uses the 17,646 vOTUs with
    estimates. Neither denominator is printed with its statistic.

12. **Phage vConTACT2:** 94.5% is a cluster-level result over 2,516
    EQ-containing clusters, but only 17,942 EQ vOTUs entered and only
    8,954/18,604 (**48.13%**) reached a numbered cluster.

13. **Phage spacer network:** 1,632 raw spacer links collapse to 885 unique
    host-vOTU edges; the outline does not explain the switch. The co-host
    concordance tests use only 154 multi-host vOTUs, an unstated denominator.

14. **Phage pangenome comparison:** the “~100% accessory overall” baseline
    comes from only 12 genera, not the 135 genera in the defence-repertoire
    analysis.

15. **Phage candidate recurrence:** the same printed 455 also has the
    unstated member-protein versus distinct-island aggregation described
    under (a).

### (d) Internal inconsistency between two places

1. **Sample count:** MAG Table 1 says 547; the phage paper repeatedly says
   550; canonical metadata contain 551. These are defensible subsets only if
   labelled, which they are not at the claims.

2. **“MAG proteins” denominator:** the MAG paper uses 51,920,778 Bakta
   proteins for its catalogue and 55,380,378 CheckM2 coding sequences for its
   no-KO fraction without distinguishing them.

3. **Endemicity identity:** the same 23,382-query percentages are labelled
   soil-MAG dark families in the MAG paper and “structurally-dark mobilome
   proteins” in the phage paper.

4. **Anti-defence table:** its eight rows sum to 478 while its total row says
   483; five omitted systems are not represented by an “other” row or subset
   label.

5. **Borg screen:** the phage main text/source has 9 qualifying long contigs,
   while the outline and `results/README.md` say 10.

## Complete claim ledger

Repeated occurrences are grouped below only when their object, denominator
and source are identical. All figure bars and table rows are included as
groups; the exact component values are recorded in Checks 5–20 above.

### Empty Quarter MAG/proteome paper

| Claim bundle | Recomputed from source | Verdict |
|---|---|---|
| Survey scale: 15,469 MAGs; ~51.9M proteins; 60×5; 547 samples | 15,469 and 51,920,778 exact; metadata give 61 site IDs, 5 trips, 217 cells, 551 samples, 550 binnable; `sample_metadata.tsv`, `binnable_samples.txt`, `env_matrix.tsv` | Counts partly PASS; **(a), (b), (d)** for design/sample labels |
| A0 taxonomy: 14,752 bacteria + 717 archaea; 44 phyla (37/7); 0/0 novel phyla/classes; 12/297/3,677/10,815 lower ranks; HQ 3/80/996/2,585; 3,664 eligible; phyla 6,788/3,740 | All exact from GTDB-Tk summaries, `novelty_summary.tsv`, `seqcode_candidates.tsv` | PASS |
| A0a: 17 domain-unplaced, all MQ rhizosphere; 0 HQ; 3 ar53; 51.4–66.0% complete, 2.7–9.7% contaminated | 17; all MQ/PRr; 3 ar53; ranges 51.38–65.96 and 2.71–9.74 | PASS |
| A1 marker table: 5,292/2,609 group-1l; 986/893 group-1h; 5,645/3,828 form-I; guild core 1,204 (~8%) | Exact; 1,204/15,469=7.7833%; `hydB_class_summary.tsv`, `coxL_form_summary.tsv`, `guild_core_mags.txt` | PASS |
| Abstract ~3,400 H2/CO; “dominant” form-I | H2=3,386, CO=3,828, intersection=1,204, union=6,010; CoxL total 13,675 includes 8,003 ambiguous | **(c)** |
| A1 aridity and figure | rho=-0.171589, p=0.1961, n=58; all 16 phylum-marker bars exact; `aridity_summary.tsv`, `atmo_guild_by_phylum.tsv` | PASS; source summaries separately contain 2,602/885 versus main table's 2,609/893 marker-MAG counts |
| A2 NER/BER/NHEJ | 15,382/15,469=99.4376%; 15,398=99.5410%; 12,122=78.3632%; `repair_summary.tsv` | PASS |
| A2 photolyase | 22.4% surface vs 13.1% deep; every integer table compatible with rounded rates has Fisher p<2.9e-20 | Numerically PASS; **(b)** exact cells/test absent |
| A3 osmolyte null | 113/9,425 vs 145/13,342; OR=1.10444, p=0.44587 | Numerically PASS; **(b)** wrong ledger source |
| A4 BGC prevalence and figure | Terpene 9,984/15,469=64.541987% -> 64.5%; NRPS 3,914=25.3022%; ectoine 2,626=16.9759%; all 12 bars exact; `bgc_product_summary.tsv` | Terpene **(a)**; remainder PASS |
| A5 resistance/pathogens and figure | Beta-lactam 3,816/15,469=24.6687%; As=2,423; Hg=1,077; pathogens (5+2)/15,469=0.04525%; all 12 bars exact | PASS |
| B1 catalogue and representative partition | 51,920,778 Bakta proteins -> 23,416,070 reps = 2,839,030 no-hit + 20,577,040 hit; structural partition 2,533,446 + 305,584; 89.2363% no-hit | Arithmetic PASS; **(b), (c)** |
| B1 family ladder | 334,810 multi-member families over all 2,839,030 sequence-dark reps +1,540,730 singletons; 68,808 tier1-derived families over another 11,418,079-protein cohort | Counts PASS; **(b), (c)** |
| Separate structural-tier trap control | Tier1 10,534,391/11,418,079=92.26%, explicitly fragment-inflated in source; tier2 10,113,496/13,088,646=77.27%; tier2 has 156,743 families, only 710 >=10 and 11 >=50 | Not printed as competing manuscript values; paper avoids the “typically large” error, but omits tier1 input caveat for 68,808 |
| B1a model inventory | 68,808 reps; 68,527 PDBs; 63,935 pLDDT records =63,654 usable+281 too long | Counts PASS; **(b), (c)** for pointer and 4,592-record attrition |
| B1a six thresholds | 15,248/63,654=24.0%; 9,258/50,154=18.5%; 5,514/33,524=16.4%; 2,836/17,196=16.5%; 1,164/5,709=20.4%; 15/187=8.0%; `fold_validation_summary.tsv` | PASS |
| B1a medians/length strata | pLDDT .440/.520; length 182/153; 26.3% n9,660; 5.4% n1,839; 2.7% n3,421; 4.8% n1,831; 4.9% n445; exactly 300 confident >=100-aa models | PASS; corrected paper does not print withdrawn “10,007” or “~12” |
| B2 constraint | Dark n22,767 median .391, control n21,579 .2859; 22,619/22,767=99.3499% below 1 | PASS |
| B2/B3 occupancy-selection | <=3 sites n1,809 median .3324; >=30 n10,498 median .4201; MWU 5.31784e-69; rho=.175681, p=3.18919e-157 | Printed values PASS; **(c)** group sizes absent |
| B3 occupancy/abundance | Raw 535/23,382=2.2889% single-site and 10,669=45.6291% >=30; medians 26 sites/5 trips; 307=1.3130% <=1 sample; 23,261=99.4825% peak>=10; median peak584.75, prevalence52/551 | Printed rounded values PASS; **(c)** denominators absent; generator's integer 543 includes 8 zero-mappable families |
| B3 ordination/PERMANOVA figures | PCoA 27.7216/7.3196%; variance 2.0/3.8/4.7/3.3/6.6/5.1/0.4%, p=.001 except V=.325 | PASS; **(c)** 21,961 vs 21,923 family universes |
| B4 dN/dS>1/FUBAR/Boltz | 148/22,767 pairwise-NG86 tail; HyPhy 65 families, ~12 powered; lectin values trace to `boltz_affinity.tsv` | Threshold exists; **(c)** method/aggregation mixed |
| B5 ten EC calls / three targets | 10/39/51 review verdicts; three targets from 73 folded proteins at pLDDT>=.60 and pocket>=.4 | Counts PASS; **(b), (c)** |
| B5 calibration | Current ten span .219–.432, median .262 | **(a)** |
| C1 436-sample core | NSTI .1591; PICRUSt KOs 10,543 constant; MAG median 4,482.5; recall .9327 | Printed rounding PASS; **(c)** constant described as median |
| C1 no-KO | 31,413,266/55,380,378=56.7227%; per-MAG median56.05%; figure groups 651/10,815/3,677/326 | Arithmetic PASS; **(c), (d)** denominator differs from Bakta catalogue |
| C1 threshold figure | Hybrid ratios 2.39619,1.97170,1.78596,1.48413,1.11480,.573007,.0165484; true shared ratios include 14.6333 at 1e-5 and 3.8846 at 1e-3 | Source figure reproduced, but manuscript interpretation **(a)** |
| C1 ten functions | rho range -.123 to .417, median .086 | PASS |
| C2 cross-feeding | 16 MRO files, range .549746–.777843, median .688896; 27 guild/289 non-guild, donor medians232/243, p=.693779; fractions .503/.495, p=.387733; completeness rho .0344/.0605 | Printed values PASS; **(c)** group sizes absent |
| C3 homology sets | Of 23,382: MGnify absent8,552=36.5751%; GEM11,075=47.3655%; deserts7,469=31.9434%; absent all5,011=21.4310% | Percentages PASS for this cohort; abstract object **(a), (c)** |
| C3 depths/selection | GEM122,978,191 ->123M; deserts147,531,178 ->148M; MGnify718M unmanifested; endemic/shared dN/dS .393/.390, p=.00709 | Desert/selection **(a)**; MGnify **(b)** |

### Empty Quarter phage/mobilome paper

| Claim bundle | Recomputed from source | Verdict |
|---|---|---|
| Survey scale: 550 metagenomes; 15,469 MAGs; 185 isolates; 118 cultures | 550 `binnable_samples.txt` rows and 15,469 census rows; no isolate/culture manifest | First two PASS, cross-paper **(d)**; latter **(b)** |
| Core catalogues: 18,604 vOTUs; 312,567 pOTUs; Chao2 99.91%/99.99% | 18,604 rows; 18,604/18,620.8=99.9098%; 312,567/312,589.3=99.9929% | PASS |
| Controls | 3,182/331,171=.96083%; 3,722/331,171=1.12389% | Printed values PASS; **(b), (c)** for “7”/sensitivity and denominator |
| CheckV table | 962 complete+955 high+3,037 medium=4,954 (26.6287%); low12,692; undetermined958; MIUViG1,917; complete1,474; median27.7 over17,646 | Counts PASS; median **(c)** |
| Lifestyle/IMG novelty | 7,357+5,713+5,534=18,604; 18,251/18,604=98.1026%, with353 hits | PASS |
| Taxonomy/host | Named genus2,494/18,604=13.4057%; lineage13,583=73.0112%; assigned4,135+unassigned14,469; all host-phylum figure bars exact | PASS |
| vConTACT2 | Input22,476=17,942 EQ+4,534 refs; EQ VCs2,378+138=2,516; 94.514%; 8,954 EQ in numbered VCs, 7,970 novel-only | Printed cluster result PASS; **(c)** catalogue coverage |
| TerL outline | 3,009 EQ+4,570 refs=7,579 tips | PASS |
| DRAM-v outline | Potential AMGs1,965 on1,537 vOTUs; CAZy191; categories1,040/586/56/43 and smaller rows exact | **(a)** for printed 2,928 carrier count; other counts PASS |
| pOTU catalogue | Direct FASTA: n312,567; min2,500; median4,635; p95 23,812.4; max643,799; >=10kb57,751; total2,500,595,314 bp; input584,654 | PASS except max rounding **(a)** |
| External GSPR | Article: 6,860 data sets, 98,728 collected sequences, 57,265 groups | Counts PASS; comparison method/rationale **(a)** |
| Defensome/table/figure | 264 types, 25,096 systems; 7,959/15,469=51.4513%; all 12 table rows and 15 figure labels exact | PASS |
| Defence islands | 24,776 positioned; 3,242=13.1% island,21,534=86.9%; Gabija55/170=32.4%, pAgo57/181=31.5%, CBASS83/310=26.8%, Cas161/2,348=6.9% | PASS |
| Retrons | 357 systems; 305 carrier MAGs | Counts PASS; isolate/culture absence **(b)** |
| Spacer network size | 885 unique pairs;361 MAGs;637 vOTUs; raw `n_spacers` sum1,632 | PASS; raw/collapsed switch **(c)**; 3,281 iPHoP **(b)** |
| Network statistics | Q .965 vs .9132,z9.98,p1/1001; NODF .733 vs1.065,z-15.92; phylum .862/.212; site .253/.027=9.37×; iPHoP26/27=96.296% | PASS; co-host n=154 **(c)** |
| Load | Unblocked n512: rho=.186,p=1.9e-5 density; rho=.030,p=.49 richness | Source values PASS; retained blocked negative **(b)** |
| Pangenome | n135; core median0/mean.296; accessory .955319; Jaccard .679052; genome-wide baseline n12, mean.9985; outline n132,64% positive, median rho=.077 | Printed values PASS; main baseline **(c)** |
| Microdiversity outline | Defence pN/pS .314199 n149 vs .231662,p=.010759; FST .257559 n78 vs .183448,p=.013131; anti-defence .140750 n177,p=2.8097e-7; diagnostic omega .6131/88% | PASS |
| Mobilome carriage table | Phages419 systems/361 of18,604=1.9404%; pOTUs9,780/8,221 of312,567=2.6302%; >=10kb5,608/4,199 of57,751=7.2709% | PASS |
| Mobility partition | 33+466+312,068=312,567; no-marker99.8404%; 11/33=33.333%; p95 length23.8kb | PASS; small n=33 stated |
| Length-matched mobility effect | Conjugative 11/33 vs matched3/33=3.6667×; mobilizable31/466 vs16/466=1.9375×; both p=1/2,001=.0005 | Qualitative claim supported; **(b)** saved output absent |
| Outline system-count “enrichment” | 9,780/419=23.3 and5,608/419=13.4 are raw count ratios | **(a)** as per-bp/enrichment labels |
| Anti-defence | 483 rows on455 vOTUs; multiplicity432×1,19×2,3×3,1×4; displayed types total478+5 omitted=483 | Carrier **(a)**; table **(d)** |
| Mobilome structural novelty | 316,650 hits+296,510 no-hits=613,160; no-hit48.3560% | PASS |
| Operon chain | 541,569 unknown proteins;364,544 families;455 with >=2 island-member proteins;412 with >=2 distinct island MAGs | First two PASS; recurrence **(a), (b), (c)** |
| Seven-row candidate table | Armada6/1,645/2.3e-7; CBASS5/1,391/4.8e-7; Gabija same; IIG4/1,138/1.2e-6; McrBC same; Cas5/464/1.3e-11; RloC5/126/5.1e-16 | All exact from `candidate_defense_families.tsv` |
| Structural shortlist/locus | 30 rows=29 novel+1 known; OR1,137.98,p1.177e-6,4 MAGs; pLDDT.72/.71; lengths170/62; five pockets,best8.21 | PASS; outline Boltz raw **(b)** |
| “Mobilome” endemicity | Percentages 89.1/62.8/55.9/51.9/41.7/31.9/47.4 and depths147.5M/123M exactly reproduce the MAG 23,382-query source | Arithmetic PASS for another cohort; object **(a), (d)** |
| Long-contig/Borg screen | 9 rows, length600,387–956,840, GC32.7–40.0; zero Methanoperedens in GTDB census | Main “0;9” PASS; outline10 **(a), (d)** |
| Four domain-unplaced bins | Lengths710,481–1,027,754; contigs48/75/221/157; N50 4,679–24,157; CDS901–1,425; hypothetical88.70–98.95%; tRNA6–73; GC29.1–50.1 | All PASS except printed CDS minimum **(a)** |
| Three figure PDFs | Lifestyle 18,604/7,357/5,713/5,534; host4,135/14,469 and phyla; defensome15,469 and all carrier counts | PASS |
| External positioning table | Hwang11/84; Ma1,824/80,750/96.7%; Li607/29,031; Han191/13,336; GSPR6,860/57,265; Meaden1,075; 18,604/84=221.48× | PASS; GSPR interpretation separately **(a)** |
