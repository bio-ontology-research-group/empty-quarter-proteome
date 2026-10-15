# Paper numbers — the registry

**Purpose.** The single human-readable record of every number the paper uses, in the paper's own order.
One format throughout: tables of `number | value | source`; caveats are footnotes under the table they
qualify. Convention: each section's basis line names its workbench folder once; table sources are relative
to it. No dates and no job IDs here — time and change history live in `CORRECTIONS.md`, forensics live in
the workbench README the source names (exception: a still-RUNNING job keeps its ID until the result lands).
Numbers computed but not used in the paper stay in their workbench README and are not registered.
Machine-checked mirror = `reproducibility/claims.tsv` (PI-owned, rebound from this file at end-game).

Paths: `M = /ibex/scratch/projects/c2014/rob/metagenomics` ·
`P = /ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites`

**The pipeline in one line:** two binnings (single_sample + self_supervised) → HQ+MQ union (37,949), kept
for protein-mapping → DAS_Tool refinement (14,426) → dRep ANI95 → final catalogue (5,438 species reps,
R232). Assembly proteins are catalogued separately (753M → 449M reps) and the dark set is intersected with
the union to get the genome-resolved numbers.

---

## Samples and cohorts  (the denominators everything below divides by)

| cohort | value | source |
|---|---|---|
| canonical soil samples (headline N) | **551** (5 expeditions) | `M/metadata/sample_metadata.tsv` |
| assemblies | **549** (1 failed + 1 empty): `51Sr1` has no contigs file, `F5Dr1` has a zero-byte one | `P` assembly tree |
| sequencing depth per sample | median **21.58 Gbp** / 156.6 M reads after filtering, IQR 16.53–25.08, range 0–475.81; 3 samples below 2 Gbp (`F5Dr1` 1,232 reads, `F45PRr3`, `F13Sr1`); median **48.5%** of reads map to the 5,438 | fastp `*_report.json` + coverm logs |
| DAS_Tool inputs → samples producing refined bins | **542** in → **539** out (3 yielded none) | `M/05_refine/sample_list.txt`; `*_DASTool_contig2bin.tsv` |
| genome-mapping cohort (§2.1–2.4 abundance) | **548** usable of **569** profiles produced | `M/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/` |
| family-member mapping cohort (§1.4/§2.5) | **551** (complete) | `M/06_darkmatter/darkfam_basis_clean/results/coverm/` |
| compartment split (canonical / profiled) | 117 / 166 / 268 · 116 / 165 / 267 (Surface / Deep / root-adjacent) | `M/metadata/sample_metadata.tsv` |
| **sites**; site-paired test set | **60** sites (61 numeric labels); 45 Surface+Deep sampled; **44** both profiled | `M/metadata/sample_metadata.tsv` |
| XRF chemistry subset / climate axes | **248** / **453** profiled samples | `M/metadata/` |

Footnotes: gene prediction ran on 552 assembly files (a Methods-level basis only). Three canonical samples
lack genome profiles (51Sr1, F13PRr1, F5Dr1); the family-member mapping covers all 551 (F13PRr1 landed
2026-08-14, matrix refreshed — every §1.4 digit unchanged).
⚠ **60 is the site count; 61 is a count of field-notebook LABELS** (settled 2026-08-17 against the two
companion 16S papers — both say 60 primary sites on a 1,043 km transect, plus 10 named locations; the
data descriptor shows identifiers 61–64 are trip1-only aliases of *named* locations). Of our 61 labels,
59 are numbered transect sites (56 absent) and 2 are aliases (62, 64; trip1 only; 4 samples).
**Do not "correct" 60 to 61** — that was attempted and reverted (CORRECTIONS 2026-08-17).
Expeditions: 17–20 Mar 2023 · 14 Jul 2023 (heat-truncated to 8 sites → only 37 samples) ·
17–21 Feb 2024 · 19–21 Aug 2024 · 12–15 Oct 2025 (31 months).
Median unmapped read fraction ~51%, balanced between Surface and Deep (§2.3/M.10) → paired compartment
contrasts are unbiased; root-adjacent comparisons use the mapped denominator.
"Detected" = relative abundance > 0, i.e. ≥10% genome breadth — CoverM's never-passed
`--min-covered-fraction` default, defined once in **M.5** (inherited, not chosen).

## §1.1 Catalogue and taxonomy

Basis: `M/05_refine/` (DAS_Tool 1.1.7-refined → dRep 3.6.2 ANI95). ⚠ **Taxonomy has a different
basis:** the R232 census comes from `M/rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/`
(GTDB-Tk **2.7.2** / **r232**, 2026-07-31), **not** from `_gtdbtk_refined/`, which is the superseded
GTDB-Tk 2.4.1 / **r226** run. Rows below marked `_gtdbtk_refined/` for taxonomy are read from the R232
screen; path corrected 2026-08-17 (`CORRECTIONS.md`).

| number | value | source |
|---|---|---|
| refined bins | **14,426** (4,225 HQ / 8,209 MQ / 1,992 LQ) | `_checkm2_refined/quality_report.tsv` |
| dereplication input | **12,433** clustered (of 12,434 passing HQ+MQ filter) | `_drep/out/data_tables/{genomeInfo,Cdb}.csv` |
| species representatives (ANI 95%) | **5,438** (2,259 HQ + 3,179 MQ) | `_drep/` |
| named phyla | **51** (Actinomycetota ~32%, Pseudomonadota ~24%) | `_gtdbtk_refined/` |
| first-empty-rank census (mutually exclusive) | **0/0/0** domain–phylum–class · **11** order (6 HQ) · **142** family (64 HQ) · **1,343** genus (530 HQ) · **3,755** species (1,555 HQ) · 187 assigned through species | `_gtdbtk_refined/` |
| order-gap adjudication | **10 adjudicated → 0** novel-rank proposals; 1 low-MSA candidate in recovery | adjudication records (path to confirm) |
| R232-vs-R226 reclassifications | 1,898 | `_gtdbtk_refined/` |
| Thermoplasmatota representatives by class | **237** = 185 `c__SW-10-69-26` + 51 `c__E2` + 1 `c__Thermoplasmata`; 0 named genera | `M/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv` |
| novelty of the 237 (GTDB r232) | **237 of 237 new species** (no reference within the species radius); **135 unplaced at genus** (126 in a placeholder family, **9 in no family**); 102 in placeholder (MAG-derived) genera | same, columns `classification`, `closest_genome_reference` |
| GTDB r232 species representatives in those classes | **21** (`c__SW-10-69-26`) · **197** (`c__E2`) — write "species representatives in GTDB r232", never "genomes worldwide" | `M/rebuild_20260729/references/gtdbtk_r232/taxonomy/ar53_taxonomy_r232_reps.tsv` |

Footnotes: HQ+MQ filter = completeness ≥50, contamination ≤10 *inclusive* — two bins sit at exactly 10.0
(pass ≤10, fail strict MIMAG <10); state the inclusive rule in Methods M.2/M.3. One filter-passing bin was
dropped inside dRep (`25PRr2__self_SemiBin_2_sub`), so 12,433 is the count `Cdb.csv` membership analyses
see. The low-MSA candidate is `48Sr2__single_SemiBin_167_sub` (46.49% MSA).

## §1.2 Sequence-dark proteome

Basis: `M/04_protein_catalog/` (assembly universe) and its `darkmatter_clean/` (dark calls).

| number | value | source |
|---|---|---|
| assembly proteins | **753,305,533** (pyrodigal, 552 assemblies) | `assembly_proteins/` |
| non-redundant reps (90% id / 80% cov) | **449,265,306** | `all_linclust90/` |
| UniRef50-dark | **69,938,628 = 15.6%** (hit 379,326,678; uncovered 0 — 100% classified) | `darkmatter_clean/04_catalog_dark/` |
| confident dark (≥100 aa) | **16,971,287 / 315,606,983 = 5.4%** | `darkmatter_clean/06_dark_lengths/` |
| MAG-mapped proteins / reps | 88,247,222 / 43,423,322 | `darkmatter_clean/05_mag_dark/` |
| MAG-mapped dark reps | **4,292,446 = 9.89%** (34% ≥100 aa, 44% complete ORFs) | `darkmatter_clean/05_mag_dark/` |

Footnotes: the headline is the **rep** count (4,292,446); the two protein-count variants measure different
things (`CORRECTIONS.md`). Length, not completeness, discriminates darkness: dark 76% <100 aa vs all 30%;
incompleteness 86.7% vs 87.5% — not discriminating.

Dark definition (for Methods): a rep is sequence-dark iff `mmseqs search` vs **UniRef50** at
**`-s 7.5 -e 1e-3 --max-seqs 10`** returns **no hit** — pure E-value, no coverage/identity filter. Phrase as
*"no UniRef50 homolog"* (UniRef50 is cultured-biased → metagenomic novelty is carried by the Part 3
comparison; MMseqs misses remote homology → the structure layer is that check).

## §1.3 Structural dark

Basis: `M/04_structure/00_structure_clean/` — 3Di triage of the ≥100 aa genome-resolved dark
(`03_search/`), structural clustering (`04_cluster/`), ESMFold + full-coordinate validation
(`05_validate/`). References: AFDB v6 / UniProt 2025_03 (241,070,489) + PDB release 250101 (324,204),
both stages. Report both TM normalisations with the control floor — never a raw fold count.

| number | value | source |
|---|---|---|
| queries (≥100 aa genome-resolved dark) | **1,466,943** | `03_search/target_ids.txt` |
| triage retrieval: AFDB-any / PDB-any / both | 463,316 (31.6%) / 232,793 (15.9%) / 226,812 (15.5%) | `03_search/SUMMARY.txt` |
| triage: no hit in either (candidates) | **997,646 = 68.0%** — a triage, not a result | `03_search/SUMMARY.txt` |
| structural clusters / admitted fold families (≥3 members) | 901,638 / **17,221** (80,205 members) | `04_cluster/SUMMARY.txt` |
| excluded, no claim | 851,393 singletons + 33,024 pairs | `04_cluster/cluster_sizes.txt` |
| models (one per family) | **17,221** = 13,770 folded + 3,451 reused | `05_validate/SUMMARY_05a.txt` |
| confidence-qualified models **Q** (≥100 aa, pLDDT ≥0.70) | **2,838** of 17,221 (16.5%) | `05_validate/per_rep_structures.tsv` |
| **K** dual no-hit, query-normalised TM <0.5 | **1,031 = 36.3%** (95% CI 34.6–38.1) | `05_validate/summary_05c.json` |
| control floor, query-normalised | **101/3,013 = 3.35%** (95% CI 2.77–4.06) | `05_validate/per_control_structures.tsv` |
| candidate rate vs its floor | **6.6×** (CI 5.0–9.1) standardised on model confidence, **8.1×** (CI 6.0–11.5) on confidence × length. Crude 10.8× is superseded: controls average pLDDT 0.867 vs 0.775 and their error rate falls 5.49→1.53 % across that range | `05_validate/stratified_floor.tsv` (2026-09-17) |
| dual no-hit, alignment-normalised TM (sensitivity) | **95 = 3.35%**; floor 68/3,013 = 2.26%. ⚠ **NOT a positive result**: standardised on confidence × length it is 1.38× with CI 0.87–2.41, crossing 1. The shipped "1.5×, p = 0.011" is superseded and must not be quoted | `05_validate/stratified_floor.tsv` (2026-09-17) |
| locally familiar vs fully covered (of Q) | 96.7% have a matching aligned region; 63.7% covered end-to-end | `05_validate/per_rep_structures.tsv` |
| member proteins in the K families | 4,880 | `05_validate/fam_reps.tsv` |
| **Pfam-A 38.2 on the K members** (basis `M/04_structure/pfam_candidates_clean/`) | **427 of 4,880 (8.75%)** with ≥1 domain; **144 of 1,031 families (14.0%)**; hits cover median **0.25** of the chain, **89.5%** leave ≥80 aa unmatched (gap median **159 aa**). Q-not-K (1,807): 5.59% / 0.56 / 41.6%; controls (3,013): 90.67% / 0.885 / 17.5% | `results/set_summary.tsv` |

Footnotes: the length gate is nearly free (17,218/17,221 ≥100 aa) — **confidence** is what the models must
survive. Candidate models are longer than other qualified models (mean 301 vs 215 aa) at indistinguishable
confidence (pLDDT 0.763 vs 0.782), so whole-protein novelty tracks multi-domain architecture rather than model
quality. The 3,684 blinded known-fold controls are the PI's set, co-searched in our own run under the identical
invocation; our floor reproduces his exactly (same controls, references and parameters). The historical B1a
calibration (68,527 models, 608/1,875) is on the old May tier1 cohort and is **superseded** — it survives only
in the explicitly historical Table and the fenced B1a scaffold.

## §1.4 Dark families

Basis: `M/06_darkmatter/darkfam_basis_clean/` — 4,292,446 genome-resolved dark reps → mmseqs 50%/80% →
admission ≥3 diverged members, ≥3 samples, ≥2 sites. Mapping cohort: 551 (Samples and cohorts).

| number | value | source |
|---|---|---|
| dark families / member CDS | **106,042 / 536,437** | `README.md` |
| excluded classes | 2,784 paralog-only · 782,616 young (1–2 variants, no claims) | `results/fam_recurrence.tsv` |
| recurrence | median **12 samples / 9 sites** | `results/fam_recurrence.tsv` |
| dN/dS dark vs control | **0.425** (96.7%<1, n=97,200) vs **0.388** (99.6%<1, n=95,647); MWU p<1e-300 | `results/selection/dnds.tsv` |
| fast tail (>1) | 3.3% dark vs 0.4% control — hypothesis-generating only | `results/selection/dnds.tsv` |
| **frame-shifted null arm** (same dark CDS read in the +1 frame, stop codons masked; seeded 10,000-family subset) | median **1.087** (IQR 0.901–1.326; 37.7% <1; n 9,148) vs the same families in the true frame **0.427** (IQR 0.295–0.580; 97.0% <1; n 9,182; full arm 0.425); paired null − true **+0.655**, null higher in **92.4%** of 9,029 families; stop codons masked 2.0% (+1 frame) vs 0.7% (true) | `results/selection/dnds_null_frame1.tsv`, `dnds_dark_subset.tsv` |
| prevalence / abundance | median **98/551** · 4 single-sample · median peak **99 TPM** | `results/fam_{presence,abundance}.tsv` |
| persistence | **79.4%** all-5 expeditions (94.8% ≥4) · 28.9% re-assembled all-5 | `results/fam_presence.tsv` + `fam_recurrence.tsv` |
| persistence, length-stratified (presence rule favours short genes) | families ≥250 aa: median 55/551, **67.1%** all-5 · <60 aa: median 201/551, 89.7% all-5 | `results/fam_presence.tsv` × `darkmatter_clean/06_dark_lengths/` (CORRECTIONS 2026-08-17) |
| phylum attribution (**R232**) | **97.9%** single-phylum of 78,727 attributable (74%); Actinomycetota 32% / Thermoplasmatota 19% (denominator = all attributable) | `results/fam_phylum.tsv` |
| threshold robustness | ≥3-member clusters 108,826 (50%) → 130,843 (30%) | robustness log |

## R2 cultivation control (cultured genomes of the same basin)

Basis: `M/06_darkmatter/cultured_comparison_clean/` (promoted 2026-09-20 from `alelopezv/eqmags_poc/05_rakpr_overlap`;
rerun on all **722** RAK-PR genomes 2026-09-21, jobs 52196463–52196465). Cultured genomes = 722 (all Bacteria), all
with Bakta proteomes (651 in `eqmgc/results/annotation/bakta`, 71 in `EmptyQuarter_Data/cultures/annotation/thermophilic/bakta`);
3,024,784 proteins; median 2 contigs (49% single-contig). Darkness on the cultured side is per protein (Bakta), on the
catalogue side per 90/80 cluster representative; both by the identical UniRef50 2026_01 rule, all hit tables from the
rerun (`tmp/uniref_hits/`).

| number | value | source |
|---|---|---|
| species shared with a catalogue representative (skani 0.3.1, ANI ≥95, AF ≥30) | **189 of 722** cultured genomes (26.2%); **25 of 5,438** representatives, all Bacteria; 203 pairs; the 189 resolve to **23** representatives on best ANI | `results/summary.tsv`, `shared_species.tsv` |
| cultured dark rate | all 722: 48,741 / 3,024,784 = **1.61%** · 189 shared-species genomes: 6,713 / 777,916 = **0.86%** · unshared 533: 1.87% · by batch: soil 0.77, enrichment 1.62, isolate 2.05 | `results/step3_cultured_dark_rate_summary.tsv` |
| matched representatives (23 distinct) | 1,513 / 87,381 = 1.73% pooled, median 1.29% | same |
| **paired by species, n = 23** | cultured lower in **18 of 23**; median difference **−0.33 pp**; Wilcoxon **p 0.021**, sign test p 0.011; ≥10 strains (n=5) median −0.19 pp, single strain (n=9) −0.62 pp | `results/step3_paired_by_species.tsv` |
| cultured dark proteins re-found as catalogue dark reps | 2,746 / 6,713 = **40.9%** | `results/step3_dark_with_cultured_homolog.tsv` |
| catalogue dark reps with a cultured homolog (reciprocal cov ≥50%) | 4,226 / 4,292,446 = 0.10% (any hit 25,403 = 0.59%; 679 at ≥90% id in a shared-species genome) | same |
| control in the same DIAMOND run | annotated family reps 49,170 / 106,042 = **46.4%** vs dark family reps 187 / 106,042 = **0.18%** (263-fold) | same |
| lengths | catalogue dark reps median **70** aa (65.8% <100) · cultured dark proteins median **59** aa (78.5% <100) | `results/step3_length_summary.tsv` |

| **pair anatomy** (step3_06; reciprocal best hits ≥50% id, ≥50% cov both ways; 23 pairs, best-ANI culture) | shared proteins median **2,935** per pair (pooled 69,207 of 87,381 MAG / 85,949 cultured); shared-pair identity median **99.7%**; **annotated MAG proteins with a partner 80.3%** (median 84.7) vs **dark 15.7%** (median 14.5); from the culture side 81.2% vs 23.6%; dark < annotated in **23 of 23** pairs, both directions; on complete calls ≥100 aa 88.9% vs 15.7% (n 65,450 vs 166) | `results/pair_anatomy/` |
| pair anatomy, unique remainders | MAG-unique 18,174: median **137 aa** (shared 268), 45% partial, 51% within 200 bp of a contig end, 7% dark · culture-unique 16,742: median **210 aa**, 4.5% dark, 24% Bakta hypothetical; top products lipoprotein, phage protein, transposase, IS3, histidine kinase | same |

Superseded (651-proteome run, 2026-09-20): 1.71% / 0.99% / 16 of 23, p 0.042 / 37.8% / 3,988 / 46.5% vs 0.16% / 60 aa.

## R3 lineage darkness (per-genome dark fraction)

Basis: `M/06_darkmatter/lineage_darkness_clean/` (promoted 2026-09-20 from `eqmags_poc/01_dark_by_lineage` round two +
`03_genome_size_function`; reproduced byte-identical). Native key: protein instance → winner contig → 90/80 cluster →
dark call; n = 5,438. Dark fraction here is per **protein instance on the representative's contigs**, a third basis
distinct from 9.89% (cluster reps, 37,949-bin union) and 15.6% (cluster reps, whole catalogue).

| number | value | source |
|---|---|---|
| proteins / dark on winner contigs | 20,379,930 / 1,600,067 (global 7.85%); per-rep median **6.3%** (IQR 3.7–10.4) | `results/v2_scan_report.txt`, `v2_rep_dark.tsv` |
| R²: quality (completeness, contamination, log N50, log length) | **0.115** | `results/v2_model_report.txt` |
| R²: + phylum / + genus | **0.768** / **0.835** | same |
| variance partition | unique lineage **0.720** · unique quality **0.020** · shared 0.094 · unexplained 0.165 | same |
| novelty gradient (median dark by first unnamed rank) | phylum 0.217 (n 10) · class 0.175 (560) · order 0.108 (272) · family 0.083 (1,006) · genus 0.065 (1,613) · species 0.033 (1,855) · fully named **0.017** (122) | same |
| named vs placeholder genus, medians (MWU p) | Actinomycetota 0.022 vs 0.075 (2e-138) · Pseudomonadota 0.030 vs 0.051 (1e-56) · Bacteroidota 0.053 vs 0.068 (1e-12) | same |
| phylum increment after novelty | +0.653 → **+0.340** (52% retained) | same |
| phylum residuals after quality model (n ≥20) | Thermoplasmatota **+0.21** · Cyanobacteriota +0.14 · Bdellovibrionota +0.10 · Halobacteriota −0.05 · Bacillota −0.05 · Pseudomonadota −0.04 · Actinomycetota −0.02 | `results/v2_phylum_dark.tsv` |
| dark vs corrected genome size | slope **−0.096** per tenfold, R² **0.067**, Spearman −0.117 (bigger = less dark) | `results/model_report.txt` |
| KO richness vs size (Extended Data) | R² 0.516 (size) → 0.760 (+quality) → 0.910 (+phylum); +413 KOs per doubling; Halobacteriota −298, Thermoplasmatota −253 KOs below the line | same |

Footnote: 50.8% of representatives below 90% complete; genera <5 reps pooled; no phylogenetic mixed model (descriptive).

## R4 archaeal classes (Thermoplasmatota)

Basis: same folder. Class counts in §1.1 table above.

| number | value | source |
|---|---|---|
| median dark fraction | Thermoplasmatota **0.310** (n 237) · other archaea 0.046 (354) · bacteria 0.062 (4,847); 236/237 above the bacterial 95th percentile | `results/v2_phylum_dark.tsv`, `v2_matched_null_summary.txt` |
| novelty-matched contrast (class stratum) | 0.310 (n 236) vs 0.125 (n 324, all Bacteria), MWU p 7e-90 — superseded as the main contrast by the matched null | `results/v2_model_report.txt` |
| **matched null vs archaea** (caliper 1.0) | **218 of 237** paired; **0.308 vs 0.048**, **6.5x**; darker in **216/218**; worst residual SMD 0.28 | `results/v2_matched_null_summary.txt` |
| **matched null vs bacteria** (caliper 0.5) | **237 of 237**; **0.310 vs 0.069**, **4.5x**; darker in **237/237**; worst residual SMD 0.05 | same |
| matching covariates | completeness, contamination, log corrected size, log N50, log mean abundance, prevalence; 1:1 nearest neighbour without replacement | same |
| composition | species **4.4%** · abundance **0.8%** (median of per-sample shares, 549 samples) · dark families **19.8%** of 77,069 single-phylum · dark proteins **15.2%** · enrichment **4.6** (families) / **4.0** (proteins); Actinomycetota 30.0 / 52.8 / 32.8 / 21.5, enr 1.1 / 0.7 | `results/v2_dark_composition_by_phylum.tsv` |
| class tree | FastTree LG+gamma on the **r232** GTDB-Tk ar53 alignment, **1,914** taxa (237 desert + 1,647 r232 Thermoplasmatota reps + 30 outgroup); display only. The earlier tree on the r226 alignment (1,525 taxa) is kept in `results/thermo_tree_r226/` | `results/thermo_tree.nwk` (job 52210624) |
| occurrence (549 usable samples) | detected in **481 of 549**; median 21 of 237 reps per sample; summed share median **0.39%** of reads (IQR 0.06–1.01); Halobacteriota 81 of 549 (patchy, median 31 reps where present) | `results/v2_thermo_ecology_summary.tsv` |
| depth | median share Deep **0.55%** (155/165 samples) · PlantRoot 0.38% (224/267) · Surface **0.17%** (102/117); site-paired Surface vs Deep: **Deep higher at 40 of 45 sites**, median −0.38 pp, Wilcoxon p 4e-7 (38/45, p 8e-7 after dividing by total mapped); PlantRoot vs bulk **null** (26/53, p 0.66) | same |
| clean absences in the 237 (median completeness 85.4) | mcrABG **0** (also 0 in all 5,438: no methanogen recovered) · amoABC **0** (archaea overall 55/591) · rhodopsin bop/hop/sop **0** (Halobacteriota 51/269) · nifHDK 0 (all archaea 0) · photolyase **5** · ectoine ectC **1** · hydrogenase K06281 17 · coxL K03520 12 · trehalose otsAB 30 | `results/v2_thermo_ecology_by_rep.tsv` |
| aerobic respiration in the 237 | cytochrome c oxidase coxA K02274 **168 of 237** (SW-10-69-26 124/185, E2 43/51) · cytochrome bd 38/237 | same |
| encoded metabolism of the 237 (share of genomes; DRAM step-form modules ≥2/3 steps; comparators 1:1 completeness-matched: Halobacteriota 181 pairs, other archaea 216, bacteria 237) | EMP glycolysis **82.3%** · gluconeogenesis 81.0 · archaeal non-oxidative PPP **82.7** (Halobacteriota 9.9) · oxidative PPP 0.4 (75.1) · Entner–Doudoroff **0** · TCA 83.1 via 2-oxoglutarate:ferredoxin oxidoreductase 96.2 · glyoxylate cycle 79.3 (34.8) · complex I 70.9 · complex II 77.6 · coxA+coxB 58.6 · V/A-ATPase 91.1 · branched-chain amino-acid ABC **78.9** · peptide/Ni ABC 54.9 · multiple-sugar ABC **2.1** (72.9) · extracellular serine protease AprX 95.8 · endoglucanase 79.3 · chitinase 44.3 · GDH 94.5 + GS 82.3, GOGAT 19.4 (82.3) · nitrate assimilation 2.5 (71.8) · sulfate assimilation 10.1 · sulfide:quinone oxidoreductase 34.2 (1.7) · archaellum 99.6 · CRISPR-Cas 79.3 · catalase 28.7 (84.5) · trehalose otsAB 7.2 (53.6) | `results/v2_thermo_metabolism_summary.tsv` |
| auxotrophy proxy | ilvBN/C/D carbon skeleton **1.7%** (Halobacteriota 91.2) with the terminal transaminase 93.2 → valine/isoleucine imported; median 11 of 16 testable amino acids and 5 of 17 cofactors encoded (matched bacteria 12 / 9); cobalamin 0.4, F420 0, pantothenate 0, thiamine 0.4, biotin 7.6 | same |
| the two classes | E2 (51) vs SW-10-69-26 (185): same completeness (82.5 vs 85.6, p 0.18), smaller genomes (2.89 vs 3.84 Mb, p 3e-19), 819 vs 952 KOs at matched completeness; E2 sheds importers (peptide ABC 3.9 vs 68.6%), Trp/Pro/His/Met synthesis, SOD, cas3; keeps GOGAT (68.6 vs 5.4); dark 0.283 vs 0.320 | `results/v2_thermo_metabolism_report.txt` |

## §2.1 Reference prediction bias

Basis: `M/06_reference_bias/refbias_clean/` on 5,438 reps / `rep_ko_v2`; PICRUSt2 predictions from
`M/rebuild_20260729/picrust2_current_20260803/inputs/` (checksum-pinned, byte-identical to the companion
amplicon study's own input). Comparison cohort = 448 samples with both a CoverM profile and a 16S profile.

| number | value | source |
|---|---|---|
| **no-KO ceiling** | **55.6%** of 19,186,853 protein-coding genes (8,516,190 with a KO) | `results/ceiling_summary.tsv` |
| ceiling, CheckM2 denominator | **58.3%** of 20,445,320 (0.583465; 58.4 was a rounding slip) — the like-for-like successor to the historical 56.7% | `results/ceiling_summary.tsv` |
| ceiling, HQ only | 53.9% (2,259 reps) — not a genome-quality artefact | `results/ceiling_summary.tsv` |
| comparison cohort | 448 = canonical ∩ CoverM-profiled ∩ amplicon-profiled | `results/per_sample_comparison.tsv` |
| weighted NSTI | median 0.159 (IQR 0.115–0.209) | `results/comparison_summary.tsv` |
| presence-recall of MAG KOs | median 86.0% (IQR 85.0–87.1) | `results/comparison_summary.tsv` |
| Spearman ρ, shared KOs | median 0.839 (IQR 0.813–0.858) | `results/comparison_summary.tsv` |
| MAG KO richness / sample | median 8,666. PICRUSt2 returns 10,543 in **419 of 448** samples, not every one: 29 return fewer, minimum 1,948. Medians are unaffected (recall 0.8597 vs 0.8596), but "a constant" is wrong | `results/comparison_summary.tsv` |
| KOs outside the PICRUSt2 reference | 2,517 of 11,374 (22%); median 14.0% per sample | `results/ko_universe.tsv` |
| KOs inside the reference but unpredicted | **0 in 420 of 448** samples (median 0) — the recall shortfall is the universe gap | `results/per_sample_comparison.tsv` (`mag_in_ref_not_predicted`) |

Headline the no-KO ceiling; the over-prediction ratio stays withdrawn (asymmetric denominator — see the
fenced C1 scaffold, retained as historical). Ceiling denominators are not interchangeable: the KO numerator
comes from Bakta gene calls, so 55.6% is the proteome-consistent value; CheckM2's own calls run ~6% higher.

## §2.2 Trace gases

Basis: `M/06_atmospheric/lifeonair_clean/` on 5,438 reps / `rep_ko_v2` (`M/04_functional/repko_clean/`);
abundance over the 548-sample genome-mapping cohort.

| number | value | source |
|---|---|---|
| high-affinity hydB (HydDB 1h/1l) | **968** species (18%) | `README.md` |
| form-I coxL | **1,170** species (22%) | `README.md` |
| guild (either) / core (both) | **1,860 (34%) / 278 (5%)** | `README.md` |
| 1l vs 1h proteins | 1,124 vs 445 | `classify/` |
| coxL form calls (of 4,727 proteins) | form I **1,901** · **ambiguous 2,815 (60%)** · form II/other 11 — so 1,170 species is a floor | `classify/coxL_class.tsv` |
| guild abundance | median **18.4%** whole (IQR 12.8–25.6) / **43.9%** of mapped; CO 11% vs H₂ 9% | `README.md` |
| energetics / sink | ΔG ≈ −200 / −230 kJ mol⁻¹; ~0.01 Tg H₂ yr⁻¹ (<1% global soil sink) | Supp Note 1 derivation |

## §2.3 Radiation, desiccation, salt

Basis: `M/06_metabolism/uv_survival_clean/` + `M/06_metabolism/osmolyte_clean/` on 5,438 reps / `rep_ko_v2`;
compartment tests on the 44 site-paired subset. Full grids = Supplementary Tables S3–S4.

| number | value | source |
|---|---|---|
| within-phylum photolyase control (**R232**) | replicates in **3** phyla: Actinomycetota p=0.0014 · Pseudomonadota p=0.020 · Bacillota p=0.021. **Bacteroidota p=0.30 (not significant)** — it gained 107 reps in the r226→R232 `Bacteroidota_A` merge | `uv_survival_clean/results/within_phylum_tests.tsv` |
| repair floor (% of species) | BER 99.6 · NER 99.5 · oxidative 93.5 · SOS 90.1 · NHEJ 74.5 · MMR 57.7 · carotenoid 54.4 · photolyase 17.9 | `uv_survival_clean/README.md` |
| compartment (mapped, 44 sites) | photolyase 24.3 vs 12.1 · carotenoid 70.3 vs 61.0 · NHEJ 80.0 vs 85.7 | Supp Table S4 |
| compartment (whole-community bound) | photolyase 11.1 vs 4.9 · carotenoid 31.1 vs 28.8 | Supp Table S4 |
| ecotype test | null across 324/327 spanning species (all p≥0.18) → species sorting | `uv_survival_clean/README.md` |
| osmolyte census (% of species) | any 79.5 · trehalose-any 75.4 (otsAB 66.3) · ectoine 16.5 · GG 3.6 · MG 2.8 | `osmolyte_clean/README.md` |
| salt = turnover | Halobacteriota vs S ρ +0.451 (p 7.4e-14) · Cl +0.385 · proxy +0.311 (n=248, geography-partialled; CLR/detection-robust) | `osmolyte_clean/README.md` |
| haloarchaea pathways | ectoine 5/269 · GG+MG 0/269 · trehalose 172/269 | `osmolyte_clean/README.md` |

## §2.4 Biosynthesis and resistome

Basis: `M/06_metabolism/bgc_clean/` + `M/06_metabolism/resistome_clean/` (+ `pathway_check_clean/`
verification layer) on 5,438 reps; community carriage on the 548-sample CoverM cohort.

| number | value | source |
|---|---|---|
| BGC census | 97.5% of species ≥1 cluster · 44,002 regions · contig-edge 85.3% | `bgc_clean/results/bgc_census_summary.tsv` |
| BGC products (% of species) | terpene 62.1 · terpene-precursor 83.8 (separate) · T3PKS 41.5 · RiPP-like 36.5 · NRPS 25.1 · T1PKS 22.4 · ectoine 16.5 (899) | `bgc_clean/results/bgc_product_summary.tsv` |
| BGC carriage | any-BGC median 48% whole-community / 99% mapped | `bgc_clean/results/community_summary.tsv` |
| AMR census (% of species) | any 37.4 · β-lactam 24.6 (1,335) · aminoglycoside 11.9 · rifamycin 7.2 | `resistome_clean/results/amr_class_summary.tsv` |
| AMR carriage | median 15% whole-community / 33.5% mapped | `resistome_clean/results/community_summary.tsv` |
| metal layer (STRESS) | arsenic 17.3 (943) · mercury 7.5 (409) | `resistome_clean/results/stress_class_summary.tsv` |
| arsenic gradients (geography-partialled) | Ca ρ +0.224 (q 0.011, n=248) · precipitation −0.243 (q 2.4e-6, n=453) | `resistome_clean/results/gradient_tests.tsv` |
| pathogen census | 2 named pathogen species of 5,438 (0.04%) · pathogen-containing genera 0.9% | `resistome_clean/results/pathogen_census_summary.tsv` |
| VFDB tiers (PI compute, reused) | 5,435 → 3,769 → 202 species at id ≥40 → ≥70 → ≥90 | `M/06_resistome_safety/vfdb/VFDB.md` |
| ectoine triangulation | antiSMASH×KO Jaccard 0.921 (860/895/899) · trehalose KO×gapsmith 0.795 | `pathway_check_clean/results/triangulation_summary.tsv` |

Footnotes: BGC novelty never computed (deferred list); per-species prevalence conservative vs
pangenome content (one rep per species); salinity-proxy association attenuates geography-partialled
(0.24→0.12) — cite Ca/Na + aridity instead.

## Contamination control (blanks and the two reagent contaminants)

Basis: `M/06_darkmatter/contamination_clean/` (promoted 2026-09-20 from `eqmags_poc/07_contamination_check`, byte-identical
rerun; blanks profiled with the paper's own CoverM 0.7.0 invocation and the prebuilt 5,438 index, 10% breadth gate).

| number | value | source |
|---|---|---|
| shotgun negative control SNTC | 335,190 pairs; **8.56% assigned** to exactly **2** representatives: *E. coli* `V42PRr2__single_SemiBin_33` **4.30%**, *C. acnes* `51Dr3__single_SemiBin_390` **4.26%** | `results/sntc_profile.tsv` |
| shotgun positive control SPositive | 119,683,914 pairs; **0 reads assigned, 0 representatives** > 0 (ungated 4.1% of reads align without reaching breadth) | `results/spositive_profile.tsv` |
| what SPositive is | one organism, *Hafnia paralvei* (GTDB-Tk R232, ANI 98.99, AF 0.87 to GCF_001655005.1); absent from the 5,438 by skani, reads and contigs | `results/spositive_gtdbtk/` |
| what the SNTC assembly is | 11 contigs, 11 kb, laboratory plasmid: aminoglycoside 3'-phosphotransferase, TEM β-lactamase, Rop, chloramphenicol acetyltransferase, T4 SSB at 94–100% identity | `results/sntc_uniref50.tsv` |
| *E. coli* representative in soil | 316/550 samples (57.5%); trips 72 / 57 / 24 / 26 / **100**%; DNA rho **−0.49**; replicate-concordance excess 0.004 (null band 0.30–0.44) | `results/species_table.tsv` |
| *C. acnes* representative in soil | **290/550** (52.7%); trips 67 / 35 / 47 / 60 / 29%; DNA rho **−0.47**; not on the pathogen panel | same |
| pathogen census | strict 2 of 5,438 (0.04%); suffix-blind **5** (0.09%) | `results/pathogen_census.tsv` |
| ⚠ vector sequence inside the catalogue | **10** representatives carry contigs matching SNTC plasmid contigs; strongest `T5PRr1__self_SemiBin_614` (*Mixta alhagi*) 3,173 bp at 99.94% | `results/control_contig_hits.tsv`, `results/identify/SNTC_contigs_vs_reps.paf` |
| effect on the resistome census (checked 2026-09-20) | AMR hits on SNTC-matching contigs: **2** (both aph(3')-Ia, aminoglycoside; T5PRr1 coordinate-confirmed, V40PRr2 same contig only) → aminoglycoside 649 → 647 species, **11.9% unchanged at one decimal**; β-lactam 24.6, rifamycin 7.2, phenicol 5.02, any-AMR 37.4 unchanged. Wider: aph(3')-Ia occurs in **7** species, 5 on an identical 7,838 bp contig at 100% identity (plasmid-like); blaTEM-116 in 2 species at 100%; removing all 7 aph(3')-Ia gives **11.8%**. catA hits (10 species, 53–85% identity) are native | `resistome_clean/results/amrfinder/*.tsv` × the PAF |

## R7 differential function by compartment

Basis: `M/06_metabolism/differential_function_clean/` (promoted 2026-09-20 from `eqmags_poc/04_differential_function`,
byte-identical rerun; method review in its README). Carriage = share of the MAPPED community encoding a feature
(549 samples, 5,438 reps; mapped fraction mean 47.8%, root-adjacent 50.9 vs subsurface 45.3 vs surface 44.3).
Modules = 282 DRAM KEGG step-form modules present at ≥2/3 of steps; KOs = 4,127 at 5–95% prevalence. Contrasts
paired within a site visit, medians per site. Pass = BH q<0.05 within (basis × contrast) + floor 5 pp + same sign
recomputed inside Actinomycetota alone and Pseudomonadota alone (sign control; magnitudes reported alongside).
Only the compartment axis is paper-facing; chemistry (confounded with trip: XRF 100% of trip 4, 31–53% elsewhere)
and trip (8 pairs) are not reported.

| number | value | source |
|---|---|---|
| root-adjacent − bulk, 52 paired sites | modules: 145 q<0.05 → **26** at floor → **14** pass phylum sign (3 partial, 9 fail); KOs: 2,496 → **716** → **441** (68 partial, 207 fail) | `results/compartment_tests.tsv` |
| surface − subsurface, 37 paired sites | modules: 119 → **21** → **9** (1 partial, 11 fail); KOs: 2,458 → **734** → **372** | same |
| deduplicated over both contrasts | 32 modules / 1,070 KOs clear floor + FDR; 20 / 626 also pass the phylum control on their stronger contrast | same |
| root − bulk, passing modules (pp; within-Actinomycetota / within-Pseudomonadota) | cytochrome bd oxidase **+13.0** (+17.0 / +1.9) · serine→cysteine **+11.0** (+8.4 / +3.4) · NDH-1-type NAD(P)H:quinone oxidoreductase **−10.3** (−8.8 / −1.0) · F420 biosynthesis **−10.1** (−3.8 / −0.4) · Entner–Doudoroff **+9.4** (+14.6 / +5.8) · reductive pentose phosphate G3P→Ru5P +8.5 · UDP-glucose +7.9 · complex I (nuo) **−7.8** (−11.8 / −0.03) · pyrimidine deoxyribonucleotide +7.6 (+10.9 / +5.9) · serine +7.3 · heme (glycine route) +7.2 · heme (glutamate route) +6.9 (+11.0 / +2.5) · glutathione +6.6 (+0.1 / +2.5) · PEP-carboxykinase/malate set (KEGG "CAM, dark") −5.5 (−4.8 / −3.3) | same |
| surface − subsurface, passing modules | F420 biosynthesis **−10.1** (−1.7 / −3.1) · putrescine synthesis, single KO K00797 **+8.7** (+3.5 / +2.5) · formaldehyde assimilation RuMP **+8.7** (+9.1 / +4.6) · phosphotransacetylase–acetate kinase +7.0 (+5.7 / +0.8) · dTDP-rhamnose −6.2 (−5.5 / −6.8) · tyrosine +6.0 (+0.005 / +10.2) · propanoyl-CoA→succinyl-CoA −6.0 (−2.9 / −6.5) · phenylalanine +5.7 (+0.005 / +9.7) · methionine +5.1 (+5.5 / +1.8) | same |
| largest KO effects, root − bulk | MFS DHA1 transporter K07552 +21.9 · tRNA m7G methyltransferase +20.8 · starvation-inducible Dps K04047 +19.4 · methylmalonyl-CoA mutase −18.9 · SAM hydrolase −17.6 · cytochrome bd subunit II K00426 +15.9 · **hydrogenase large subunit K06281 −15.6** (−19.8 / −3.6) · sulfate permease +15.5 · Zn/Mn transport +15.4 · iron complex permease +15.0 | same |
| known-trait check (all phylum pass) | photolyase K01669 **+8.8** surface (+12.8 / +12.8) · betaine transport K02168/K05020 +10.4/+10.3 surface · hydrogenase K06281/K06282 −15.6/−13.9 root · Hyp maturation −12 to −13 root · spore photoproduct lyase −5.5 root · trehalose transport −5.3 root · CO dehydrogenase small subunit K03518 −5.4 surface | `results/known_traits_sanity.tsv` |
| within-phylum magnitude | smaller within-phylum effect is a median **0.16** of the whole-community effect (13 of 23 module rows below 0.2); only 3 rows have both phyla clearing 5 pp alone | `results/review_diagnostics.txt` |
| robustness | size normalisation reverses no sign (47/47 modules, 1,219/1,219 KOs); reps ≥90% complete (2,675): **no sign reverses**, 21 of 26 root and 14 of 21 surface modules still clear 5 pp; two-sided-bulk visits only (64): all signs hold, 19 of 26 pass | same |
| limits | mapped fraction root − bulk **+5.6 pp** (p 1e-5); completeness +0.66 root-adjacent; KO effects offset +0.52 pp, modules centred at 0 | same |

Footnotes: one-step "modules" (cytochrome bd, nuo, K00797) are single-KO or any-of sets and are named by enzyme in prose,
not as modules; the three KEGG modules titled for eukaryotic contexts (chloroplast NDH, CAM, animal heme) are
prokaryotic enzyme sets and are named by content. All methanogenesis modules fall to 0 prevalence despite 269
Halobacteriota reps (completeness), so module absence is not evidence of absent biology.

## §2.5 Dark function and environment

Basis: `M/06_darkmatter/darkfam_ecology_clean/` — 106,042-family × 551-sample member-aware matrix
(from darkfam_basis_clean) × the shipped §2.3 env table (`env_matrix.tsv`, 547 samples; XRF subset 251,
climate 456).

| number | value | source |
|---|---|---|
| community structure (PERMANOVA, 547 samples / 105,977 families at prev≥5) | geochem R² **0.132** · geochem beyond space **0.089** · beyond space + compartment 0.092 · climate beyond space 0.030 · compartment 0.044 · space 0.079 (all p=0.001). Corrected 2026-09-17: the block now includes S, Cl and Na, which carry the family-level signal; the old 0.123 / 0.080 omitted them and is superseded | `results/permanova_blocks.tsv` |
| variance partitioning (XRF cohort, n=251) | geochem unique R² **0.106** · space 0.055 · compartment 0.057 (p=0.002). ⚠ The two-block `varpart.tsv` cited by M.7/M.13 never partitioned the environmental blocks; that claim is withdrawn | `results/varpart.tsv` (2026-09-17) |
| vanadium block | null (R² 0.003, p=0.78; within-site p 0.88) | `results/permanova_blocks.tsv` |
| permutations restricted within site (59 sites) | every environmental block keeps p ≤ 0.002 (geochem 0.002, all others 0.001); space itself is not testable within site (p 1.0) | `results/permanova_blocks.tsv` col `p_within_site` |
| sequencing depth as a block / covariate | log10 bases R² **0.011**; geochem beyond space+depth **0.089** (vs 0.089 beyond space); climate beyond space+depth 0.030; aridity beyond space 0.006 → 0.006 with depth | same; depth table `results/sample_depth.tsv` (573 fastp reports) |
| distance-decay (Mantel vs geography) | r +0.341 (p=0.001) | `results/community_summary.txt` |
| families testable (≥25-present rule) | 94,568 (XRF axes) · 100,293 (aridity) | `results/indicators_summary.tsv` |
| families tracking an axis (q<0.05 · % negative) | S 54,227 (91.7) · Ca 47,174 (92.5) · Cl 46,972 (90.4) · salinity proxy 37,280 (86.9) · Fe 24,110 (95.0) · aridity 23,031 (95.9) · K 11,935 · Na 4,603 | `results/indicators_summary.tsv` |
| trace-metal nulls | V, Mo, Ni, Co, Mn, Zn = 0 families (Cu 1) | `results/indicators_summary.tsv` |
| space-partial retention (control lat+lon, q<0.05 surviving) | **plant-root 94.1%** (n 547; 61,182 → 59,820) · aridity 85.6% · Cl 82.6% · S 82.0% · Ca 71.1% · salinity proxy 41.4% | `results/partial_summary.tsv` |
| plant-root contrast, family sets (q<0.05 of 102,005 tested) | **19,818 root-enriched** (median ρ +0.154, IQR 0.121–0.200) · **40,495 root-depleted** (−0.186, IQR −0.235 to −0.140) | `results/compartment_families_summary.tsv` |
| phylum shares, root-enriched vs all attributable (single-phylum base 77,069) | Pseudomonadota **19.2 vs 9.6%** · Actinomycetota **16.8 vs 32.8%** · Myxococcota_A 11.3 vs 4.9 · Thermoplasmatota 11.2 vs 19.8 · Planctomycetota 8.3 vs 2.4 · Bacteroidota 5.5 vs 3.0; root-depleted: Actinomycetota **52.7%**, Chloroflexota 13.2, Thermoplasmatota 11.3 | `results/compartment_families_by_phylum.tsv` |
| over-representation among root-enriched (share of set / share of all; ≥20 families) | Bacillota 3.6 · Planctomycetota 3.4 · Bdellovibrionota 3.3 · Cyanobacteriota 3.0 · Verrucomicrobiota 2.9 · Myxococcota_A 2.3 · Pseudomonadota 2.0 | same |
| all-five-expedition persistence of the two sets | root-enriched **86.8%** · root-depleted **89.8%** · all families 79.4% | `results/compartment_families_summary.tsv` |

Footnotes: counts ≈10× the June object's (universe 4.5× larger; median family prevalence 98/551 clears
the ≥25-present power gate); June's Mo signal (177 families, 100% positive) does not reproduce (0) —
never a paper claim; plant-host genus axes not rerun (amplicon cache off-IBEX; plant-root compartment
contrast retained); EC-decode layers deferred to the dark-matter-annotator sibling paper.

## §3.1 Specificity

Basis: `M/06_darkmatter/specificity_clean/` — the 106,042-family universe, queried in ONE DIAMOND pass
together with the 106,042 length/size-matched annotated control families (212,084 seqs), so every
catalogue carries a known-protein absence floor. Rule: `diamond blastp --sensitive --evalue 1e-3 -k 1`,
v2.0.15. **All numbers are dark/control pairs.** Jobs 50594992 · 50594993 · 50613399 · 50613400 · 50618838.
**June's 21.4% (5,011/23,382) is DEAD** — superseded object; the abstract's "pending SPIRE" hedge is CLOSED.

| number | value (dark / matched control) | source |
|---|---|---|
| **HEADLINE — absent from every named reference** (MGnify90 ∪ GEM ∪ deserts ∪ SPIRE, 1.258 G seqs) | **28.74% (30,479) / 0.41% (435)** — 70-fold | `results/union_specificity_summary.tsv` |
| per catalogue, % absent | mgnify90 (717,738,164) 52.68/1.36 · gem (122,978,191) 64.38/2.97 · spire (269,496,580) 39.90/0.98 · deserts_combined (147,531,178) 42.60/2.07 | `results/{specificity,desert_specificity}_summary.tsv` |
| per desert study, % absent | atacama (6,193,359) 95.61/43.04 · desert_autotroph (7,200,701) 72.81/10.48 · arid (30,569,562) 52.59/3.91 · antarctic_colddesert (49,206,550) 67.55/8.67 · mojave (54,361,006) 62.54/6.16 | `results/desert_specificity_summary.tsv` |
| known **genome** space only (GEM ∪ SPIRE, 392.5M) | 39.10 / 0.89 | `results/union_specificity_summary.tsv` |
| **content beyond depth** — gap vs a size-matched random MGnify90 slice (expected − observed, pts) | deserts_combined **+24.33** (42.60 obs / 66.93 exp) · spire **+21.30** · arid **+28.38** · gem +4.27; control gaps +0.93 / +1.23 / +2.35 / +0.32 | `results/content_vs_depth.tsv` |
| rarefaction **self-validation** (licenses the above) | mgnify90 at p=1: 52.68 observed vs 52.53 expected = **0.15 pt** | `results/rarefaction_curve.tsv` |
| **depth robustness** — vs complete MGnify90 **2026_07** (1,661,142,385 seqs, MMseqs2 `-s 7.5 -e 1e-3`) | 1,760 of 30,479 dark gain a hit (5.8%) → **27.08% still unmatched**; control 149 of 435 (34.3%) → **0.27%** | `results/depth_robustness_summary.tsv` |
| e-value/DB-size confound, measured | 2,544 families (2.40%) hit only because the desert DBs were searched split rather than concatenated | README |
| occupancy of the reference-absent set (length-stratified) | median **37 of 61 sites** vs 38 for families with a hit; holds in all five length strata (44/50 · 38/40 · 33/36 · 31/32 · 26/27) | `results/b6_absence_occupancy_by_length.tsv` |

⚠ **Get the gap's sign right.** Gap = expected − observed. A **positive** gap means observed absence is
*below* the size-matched expectation, so the catalogue leaves **less** of the dark proteome unexplained
than its size alone predicts (deserts_combined 42.60 observed vs 66.93 expected → +24.33). A negative gap
is the reverse (atacama 95.61 vs 90.94 → −4.67). **Never write "a positive gap explains more absence"** —
that states the opposite of the result.
⚠ **Do NOT union the depth leg with the DIAMOND absences** — different tool (MMseqs2) and release. It is a
robustness check; **28.74% is the quotable headline**, 27.08% is its robustness figure.
⚠ **atacama's gap is not interpretable**: its *control* gap is −29.97 because a single 35-sample study is a
narrow, redundant slice that loses to a random draw even for conserved proteins. Quote the gap for the
large, broad catalogues only. See CORRECTIONS 2026-08-17.

---

## Shared inputs  (intermediate artifact still in use)

- HQ+MQ union = **37,949** bins: **16,333** single_sample (pretrained; 4,030 HQ + 12,303 MQ) + **21,616**
  self_supervised (SemiBin2; 4,507 HQ + 17,109 MQ). Sources: `M/results/checkm2_single_sample_full.tsv` ·
  `M/03_annotation/mag_list_self_supervised.txt`.
- **Role:** the genome-resolved protein census (§1.2) maps assembly proteins to contigs in these 37,949 bins,
  across both binnings; `contig_to_bin.tsv.gz` covers 100% of them. This is the *only* live use of the union;
  taxonomy/species counts come from the refined **5,438** catalogue (§1.1), not from here.

## Superseded → replacement  (do not resurrect the left column)

| Superseded | Replaced by | For |
|---|---|---|
| union-dRep **9,253** species · R226 · **56 phyla** · 311/2,830/8,978 novel-rank | refined **5,438** · R232 · **51 phyla** (§1.1) | taxonomy / species |
| single-sample GTDB (16,333, R226); older dRep 9,211; **44 / 38** phyla | 5,438 / R232 / 51 (§1.1) | taxonomy |
| MAG-Bakta proteome **51.9M** → 23.4M reps → **2.84M** dark (12.1%) | 449M reps → 69.9M dark; 4.29M genome-resolved (§1.2) | proteome |
| single_sample **15,469 / 15,477** (strict `<10`) | **16,333** (MIMAG `≤10`) (Shared inputs) | binning tiers |
| structural "**89% novel**" / "10,007 folds" (predicted-3Di) | calibrated structural search (§1.3) | structural |
| June dark orthogroups **23,382** · dN/dS **0.391** · 52/551 · peak ~585 · 1.3% single-sample | **106,042** families · 0.425-vs-0.388 · 98/551 · peak 99 · 4 single-sample (§1.4) | dark families |
| June B3 gene–env on 23,382: Ca **4,066** · S 3,964 · Cl 3,594 · salinity 3,179 · Mo 177 positive · PERMANOVA geochem 0.066 / distance-decay 0.161 | 106,042 basis (§2.5 table): S 54,227 · Ca 47,174 · Cl 46,972 · salinity 37,280 · Mo 0 · geochem 0.123 / decay 0.341 | dark-family ecology |
| BGC on 15,469 bins: terpene **64.5%** · NRPS 25.3% · ectoine 17.0% · edge "~90%" | 5,438 species: 62.1 · 25.1 · 16.5 · edge 85.3 measured (§2.4) | BGC census |
| AMR on 15,469 bins: β-lactam **24.67%** · As 2,423 · Hg 1,077 · pathogens 7 MAGs / 0.05% · "As tracks salinity proxy r=+0.26" | 5,438 species: 24.55 · As 943 · Hg 409 · 2 species / 0.04% · Ca/Na + aridity, geography-partialled (§2.4) | resistome |
| C1 on 15,469 bins / 436 matched samples: no-KO **56.7%** · recall **93%** · MAG richness **4,482** · ratios 2.40×/1.11× and 31.41×/14.63× | 5,438 / 448 samples: **55.6%** (58.4% CheckM2 denominator) · 86.0% · 8,666 · ratios stay withdrawn (§2.1) | reference bias |

_Last updated 2026-08-14 — uniform `number | value | source` tables throughout, cohorts restored to the
top, dates/status words removed (time lives in CORRECTIONS), pending-section placeholders added so every
paper section is listed. Values unchanged._

## 2026-10-13 current review overlay

This overlay supersedes conflicting older prose. The frozen source bindings and independently rechecked values are in `reproducibility/darkmatter/current_claims.tsv` (49revisionclaims); the full change rationale is in `notes/CORRECTIONS.md` and the three review/response pairs.

| Quantity | Current basis/value |
|---|---|
| MIMAG high/medium representatives | 126 / 5,312; CheckM2-only near-complete is2,259 |
| Recurrent coding-analysis scope | 106,042families,536,437members,12.5%of genome-resolveddarkrepresentatives |
| Evaluable dark/control family medians | 97,200/.4253 and95,647/.3884; divergence/target-population dependent |
| Flagged-family calibration | 91evaluableunionfamilies,median.7622,70below1 |
| Exact base-alignment culture support | Dark-unshared516/1,105;annotated-unshared4,877/16,056;dark-shared219/236;annotated-shared62,925/65,477 |
| Culture loci with no major Bakta CDS | 344/735DNA-supporteddark;607/67,802DNA-supportedannotated |
| Culture overlap with admitted families | 65of1,341eligible darkcalls,47distinctfamilies |
| Structural no-row counts | 22dark/66control; conditional local failures73/2,816vs2/2,947 |
| Gated local structural matches | 797of1,031query-thresholdfailures have>=100alignedresiduesand>=50%querycoverage |
| Thermoplasmatota family enrichment | 4.551874,reported4.6-fold |
| Locations | 61sampledlocations:59core+2auxiliary |

Historical pooled structural-floor equivalence, accessory-content inference from protein no-hits, and the assertion that rRNA/tRNA was unassessed are superseded. Environmental results are exploratory supplementary context. Submission remains blocked on documented provenance/access/deposition gaps.
