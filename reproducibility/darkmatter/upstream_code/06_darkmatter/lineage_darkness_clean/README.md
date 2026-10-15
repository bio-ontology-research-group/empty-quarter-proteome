# lineage_darkness_clean: dark fraction per species representative by lineage, novelty and genome size

Owner: velazqal. Promoted 2026-09-20 from the exploratory folders
`alelopezv/eqmags_poc/01_dark_by_lineage` (round two) and `alelopezv/eqmags_poc/03_genome_size_function`.
This folder is authoritative; the exploratory copies are superseded.

## STATUS (2026-09-20): all reruns COMPLETE, every headline number reproduced exactly

| step | script | how it ran | wall |
|---|---|---|---|
| native-key scan | `v2_scan_cluster90.py 8` | Slurm **52131414**, 8 CPUs, 200G | 00:02:08 |
| boundary audit | `v2_boundary_check.py 8` | login | seconds |
| models + novelty | `v2_dark_by_lineage.py` | login | ~2 min |
| composition | `v2_dark_composition.py` | login | ~1 min |
| matched null | `v2_matched_null.py` | login | ~1 min |
| size / function | `size_function.py` | login | ~1 min |
| Thermoplasmatota tree | `v2_thermo_tree.py` + `v2_thermo_tree.sbatch` | login + Slurm **52210624** (2026-09-21, on the r232 ar53 MSA; 1,914 taxa), 16 CPUs, 64G. The 2026-09-20 tree on the r226 run's MSA (job 52131597, 1,525 taxa) is kept in `results/thermo_tree_r226/` | 00:05:53 |
| Thermoplasmatota known relatives | `v2_thermo_relatives_fetch.py` (login, 2 min download) + `v2_thermo_relatives.sbatch prep/blast/aai` + `v2_thermo_relatives_gtdbtk.sbatch` + `v2_thermo_relatives.py` (login, 2 s) | Slurm 2026-09-24: prep **52468152** (32 CPUs) 00:00:52 · DIAMOND array **52468153_[0-19]** (16 CPUs, 48G) 00:02:50 to 00:07:37 each · AAI **52468154** 00:00:44 · GTDB-Tk **52468155** (32 CPUs, 128G) 00:11:27 | ~20 min end to end |

Scan submit line (the script takes the worker count as argv[1]; there is no sbatch file for it):

```
sbatch --job-name=ale_v2_scan --account=pi-hohndor --partition=batch \
  --cpus-per-task=8 --mem=200G --time=03:00:00 \
  --output=logs/v2_scan_%j.out --error=logs/v2_scan_%j.err \
  --wrap="/home/velazqal/miniconda3/envs/python-utils/bin/python scripts/v2_scan_cluster90.py 8"
```

**Rerun vs copied.** Everything in `results/` is a rerun of this folder's scripts except
`reps_all.tsv` (5,438 rows) and `reps_dark.tsv` (5,359 rows). Those two are round-one tables copied
from the exploratory folder: they are the output of `step0_join_check.py` / `dark_by_lineage.py`,
the superseded rename-join scripts, which were deliberately not promoted. `v2_dark_by_lineage.py`
reads them only for the round-one continuity block (Spearman, `form` labels, the 79 dropped reps).
No paper number depends on them.

Python `/home/velazqal/miniconda3/envs/python-utils/bin/python` (pandas 1.5.0, numpy 1.26.1,
scipy 1.11.3, statsmodels 0.14.2). FastTreeMP from `/home/velazqal/miniconda3/envs/phabox2/bin`.

## Inputs (all read-only, under `W = /ibex/scratch/projects/c2014/rob/metagenomics`)

| input | rows |
|---|---|
| `04_protein_catalog/all_linclust90/all_cluster90.tsv` | 753,305,533 rep/member lines (35,934,050,788 bytes) |
| `04_protein_catalog/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz` | 69,938,628 |
| `04_protein_catalog/contig_to_bin_refined.tsv.gz`, cut to the Wdb winners | 2,485,945 contigs |
| `05_refine/_drep/out/data_tables/Wdb.csv` | 5,438 species reps |
| `05_refine/_drep/out/data_tables/genomeInformation.csv` | 14,426 refined bins |
| `05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` | 4,847 + 591 = 5,438 |
| `04_functional/repko_clean/results/rep_ko_v2.tsv` | 5,438 |
| `06_darkmatter/darkfam_basis_clean/results/fam_phylum.tsv` | 78,727 attributable, 77,069 single-phylum |
| `10_thermal/abundance/rep_relative_abundance.tsv` | 5,438 x 550 samples (549 used, F5Dr1 all zero) |
| `05_refine/_gtdbtk_refined/align/gtdbtk.ar53.msa.fasta.gz` + r232 ar53 rep taxonomy | 1,525 of 1,914 requested taxa recovered |

## Headline numbers

**Scan** (`v2_scan_report.txt`, `v2_rep_dark.tsv`). 5,438 / 5,438 reps computed, none empty.
20,379,930 proteins, 1,600,067 dark, global 0.078512. Median dark fraction 0.0631, IQR 0.0365 to
0.1038. Boundary audit: 7 boundaries, 0 lines lost.

**Model ladder** (n = 5,438, `v2_model_report.txt`). A quality 0.1147 · B = A + phylum 0.7680
(+0.6533) · C = B + genus 0.8347 (+0.0668) · N = A + novelty 0.4955 (+0.3808) · N + phylum 0.8356
(+0.3401) · N + phylum + genus 0.8749 (+0.0393). Variance partition: unique lineage 0.7201,
unique quality 0.0203, shared 0.0944, unexplained 0.1653. The phylum increment survives novelty at
52.1 %.

**Novelty gradient** (median dark by first unnamed rank): phylum 0.2170 (n=10), class 0.1745 (560),
order 0.1084 (272), family 0.0828 (1,006), genus 0.0652 (1,613), species 0.0327 (1,855), fully
named 0.0170 (122).

**Named vs placeholder genus**, three testable phyla: Actinomycetota 0.0215 vs 0.0751,
p = 2.13e-138 · Pseudomonadota 0.0303 vs 0.0514, p = 1.1e-56 · Bacteroidota 0.0526 vs 0.0683,
p = 1.24e-12. Thermoplasmatota and Chloroflexota are not testable (0 and 3 named genera).

**Thermoplasmatota**: n = 237, 236 of them first unnamed at class and 1 at order; median dark
0.309565, median Model A residual +0.214914, median completeness 85.42. At matched novelty (class)
0.3101 against 0.1249 for 324 non-Thermoplasmatota reps, MWU p = 7.17e-90.

**Matched null** (`v2_matched_null_summary.txt`). Pools: other archaea 354, bacteria 4,847.

| pool | caliper | pairs | Thermo | control | ratio | darker in | worst residual SMD |
|---|---|---|---|---|---|---|---|
| archaea | 0.5 | 128/237 | 0.3020 | 0.0459 | 6.57x | 127/128 | 0.14 |
| archaea | 1.0 | 218/237 | 0.3078 | 0.0476 | 6.47x | 216/218 | 0.28 |
| archaea | 2.0 | 236/237 | 0.3101 | 0.0476 | 6.52x | 234/236 | 0.44 |
| bacteria | 0.25 | 195/237 | 0.3110 | 0.0673 | 4.62x | 195/195 | 0.04 |
| bacteria | 0.5 | 237/237 | 0.3096 | 0.0688 | 4.50x | 237/237 | 0.05 |
| bacteria | 1.0 | 237/237 | 0.3096 | 0.0688 | 4.50x | 237/237 | 0.05 |
| bacteria | 2.0 | 237/237 | 0.3096 | 0.0688 | 4.50x | 237/237 | 0.05 |

Caliper 0.25 against archaea yields 24 pairs and is skipped.

**Composition** (`v2_dark_composition_by_phylum.tsv`). Denominators: 77,069 single-phylum dark
families of 78,727 attributable; 20,379,930 protein instances, 1,600,067 of them dark.

| phylum | n_sp | species % | abundance % | dark family % | dark protein % | enr_fam | enr_prot |
|---|---|---|---|---|---|---|---|
| Thermoplasmatota | 237 | 4.36 | 0.81 | 19.84 | 15.21 | 4.55 | 4.01 |
| Actinomycetota | 1,633 | 30.03 | 52.79 | 32.76 | 21.46 | 1.09 | 0.72 |

**Dark versus genome size** (`model_report.txt`, n = 5,438): slope -0.09565 per 10-fold corrected
size, R2 0.0671, p = 4.83e-84, Spearman rho -0.1169. Bigger genomes are less dark.

**KO richness** (n = 5,438): M1 size 0.5163 (slope +1,370.9 KOs per 10-fold, +413 per doubling) ·
M1raw on observed length 0.7049 · M2 + completeness + contamination 0.7596 (+0.2433) · M3 + phylum
0.9099 (+0.1503). Phylum residual extremes (size-only model): Bacillota +251, Pseudomonadota +238,
Deinococcota +157 against Halobacteriota -298, Thermoplasmatota -253, Planctomycetota -227.
Rank agreement size-only vs size + quality, Spearman 0.917.

**Tree** (`thermo_tree.nwk`): 1,525 taxa (237 desert Thermoplasmatota, 1,647 GTDB r232
Thermoplasmatota reps requested, 30-genome Halobacteriota/Thermoproteota outgroup), FastTreeMP
-lg -gamma, approximate-ML.

## Reconciliation against the exploratory READMEs

Every one of the 18 result files this folder regenerates is **byte-identical** to the exploratory
copy (`cmp` on `v2_rep_dark.tsv`, `v2_scan_report.txt`, `v2_reps.tsv`, `v2_model_report.txt`,
`v2_phylum_dark.tsv`, `v2_genus_dark.tsv`, `v2_unjoined_reps.tsv`,
`v2_dark_composition_by_phylum.tsv`, `v2_composition_report.txt`, `v2_matched_null.tsv`,
`v2_matched_null_summary.txt`, `model_report.txt`, `phylum_ko_residual.tsv`, `thermo_tree.fasta`,
`thermo_tree_taxa.tsv`, `thermo_tree.nwk`). No headline number differs. The FastTreeMP newick came
out byte-identical as well (77,829 bytes, 1,525 taxa, 0/1,522 bad splits).

Two presentation notes, not differences in this rerun:
- The exploratory `01` README quotes round-**one** values in its first tables (n = 5,359; e.g.
  Thermoplasmatota n = 230, median residual +0.215, JAJPHT01 n = 33). The round-two numbers above
  are the ones this folder computes and are the ones to register.
- The exploratory `03` README rounds the dark-vs-size slope to -0.0957 and R2 to 0.067; the report
  carries -0.09565 and 0.0671.

## Caveats (carried forward, unchanged)

Completeness spans 50.05 to 100 and 2,763 reps (50.8 %) sit below 90 %, so MQ genomes carry the
analysis; contamination is capped at 10. Genera under 5 reps are pooled as `other`, so the genus
increment is a floor. Dark means no UniRef50 homolog, so the signal partly tracks reference
coverage. The abundance matrix carries literal `nan` cells from failed sample columns; both
consumers coerce them to absent. Novelty is near constant inside Thermoplasmatota (236 of 237
unnamed at class), so the ordinal is partly that phylum's own indicator; the matched-stratum and
matched-null tests are the evidence. `dark family %` uses a family denominator against a genome
denominator, but enr_fam 4.55 and enr_prot 4.01 agree. KO richness is counted on the observed
genome while corrected size estimates the whole one; only the residual extremes are safe to name.

## STATUS (2026-09-20): Thermoplasmatota ecology and clean functional absences

`scripts/v2_thermo_ecology.py`, login node, ~40 s. Adds three tables:
`results/v2_thermo_ecology_summary.tsv` (531 rows, long format, numerator and denominator on every
number), `results/v2_thermo_ecology_by_sample.tsv` (550 sample columns), and
`results/v2_thermo_ecology_by_rep.tsv` (5,438 reps, marker presence 0/1).
Inputs: `10_thermal/abundance/rep_relative_abundance.tsv`, `metadata/sample_metadata.tsv`,
`04_functional/repko_clean/results/rep_ko_v2.tsv`, `05_refine/_gtdbtk_r232/` for class,
and this folder's `results/v2_reps.tsv`. Denominators: 550 of the 551 canonical samples have a
column (F13PRr1 has none); F5Dr1 is all `nan`, so **549 usable samples**; 61 sites in the metadata.

**Occurrence.** Thermoplasmatota (237 of 5,438) is detected in **481 of 549** samples, median summed
share **0.389 %** (IQR 0.065 to 1.008), median 21 of 237 reps per sample. It is the cosmopolitan
archaeal lineage here: Halobacteriota (269) reaches only 81 of 549 samples, but in those 81 it
appears in bulk (median 31 of 269 reps, max 250), so it is patchy and site restricted, not rare.
All 591 archaea 534 of 549, median 1.006 %; all 4,847 bacteria 549 of 549, median 46.110 %.

| group | samples detected / 549 | median summed share % | IQR | median completeness |
|---|---|---|---|---|
| Thermoplasmatota 237 | 481 | 0.3892 | 0.0649 to 1.0078 | 85.42 |
| Halobacteriota 269 | 81 | 0.0000 | 0.0000 to 0.0000 | 79.73 |
| Archaea 591 | 534 | 1.0062 | 0.4261 to 1.9365 | 82.45 |
| Bacteria 4,847 | 549 | 46.1101 | 36.8170 to 53.4058 | 90.58 |

**Depth.** Thermoplasmatota is a subsurface lineage. Per compartment the median summed share is
Deep 0.5469 % (155 of 165 samples), PlantRoot 0.3799 % (224 of 267), Surface 0.1743 % (102 of 117).
Site-paired Surface against Deep on per site by compartment medians, 45 of 61 sites have both:
Deep higher at **40 of 45** sites, median paired difference -0.3846 percentage points,
Wilcoxon signed-rank p = **4.18e-07**. The result is not a mappability artefact: dividing by the
total mapped per cent per sample leaves 38 of 45 sites Deep higher, p = 7.94e-07. PlantRoot against
bulk (Surface plus Deep) over 53 sites is **null**: 26 of 53 sites PlantRoot higher, p = 0.66
(normalised 25 of 53, p = 0.37). Halobacteriota cannot be tested this way: 37 of 45 and 41 of 53
site pairs have a zero difference, so its p values are uninformative. Per trip the Thermoplasmatota
median share runs 0.4725 (trip1, 190 of 210 samples), 0.3077 (trip2, 34 of 37), 0.2057 (trip3,
105 of 126), 0.4541 (trip4, 72 of 85), 0.4994 (trip5, 80 of 91).

**Clean absences** (KO identities verified against the KEGG REST API on 2026-09-20; rhodopsins are
K04641 bop, K04642 hop and K04643 sop, the only microbial entries of KEGG's rhodopsin set present
in `rep_ko_v2.tsv`; K04646 was rejected because KEGG calls it clathrin heavy chain, not a rhodopsin).

| marker | Thermoplasmatota /237 | Halobacteriota /269 | Archaea /591 | Bacteria /4,847 |
|---|---|---|---|---|
| mcrA/B/G K00399+K00401+K00402 | 0 | 0 | 0 | 0 |
| amoA/B/C K10944+K10945+K10946 | 0 | 0 | 55 | 70 |
| rhodopsin K04641+K04642+K04643 | 0 | 51 | 54 | 9 |
| hydrogenase large subunit K06281 | 17 | 39 | 69 | 1,058 |
| form-I coxL K03520 | 12 | 151 | 163 | 1,701 |
| nifH/D/K K02588+K02586+K02591 | 0 | 0 | 0 | 17 |
| photolyase K01669 | 5 | 108 | 125 | 850 |
| trehalose otsA/otsB K00697+K01087 | 30 | 175 | 206 | 3,411 |
| ectoine ectC K06720 | 1 | 5 | 14 | 881 |

Counts are "carries at least one KO of the set"; per KO counts and the "all of the set" variants are
in the summary table. The four Thermoplasmatota zeros are not equally informative. `mcr` is absent
from all 5,438 representatives and `nif` from all 591 archaea, so those two absences are properties
of the catalogue, not of this phylum. Only the `amo` and rhodopsin zeros discriminate, and they are
read against a **higher** median completeness than the comparators (85.42 against 79.73 for
Halobacteriota and 82.45 for all archaea), so they are not a recovery artefact. Photolyase 5 of 237
and ectoine 1 of 237 are near absences on the same footing.

**Aerobic respiration is present.** In the 237: cytochrome c oxidase subunit I K02274 in **168**
(SW-10-69-26 124 of 185, E2 43 of 51, Thermoplasmata 1 of 1); cytochrome bd K00425 in 31 and K00426
in 35, either in **38** (SW-10-69-26 36 of 185, E2 2 of 51). So the lineage is an aerobe or
facultative aerobe that is dark, deep and light independent, not an anaerobe.

**Caveats.** Group membership is selected by genome name, never by row position: `v2_reps.tsv` and
the abundance matrix are not in the same order, and the script asserts their index sets match.
The summed share is per cent of reads, so it moves with total mappability; the `frac_of_mapped`
columns and the `_fracmapped` paired tests are the control, and the Deep result survives both.
PlantRoot samples carry a higher total mapped per cent, which is why the raw and normalised
PlantRoot comparisons for other groups disagree in sign. Absences are absences in the KO table, so
they inherit its annotation sensitivity; for a deeply novel lineage a missing KO can mean a missing
model rather than a missing gene, which is the same caveat that makes this phylum dark.

## STATUS (2026-09-22, NOT REGISTERED): Thermoplasmatota encoded metabolism, positive content

`scripts/v2_thermo_metabolism.py`, login node, ~2 min. Reads what the 237 CAN encode, where the
ecology script read what they cannot. Module definitions and the presence rule are reused verbatim
from `06_metabolism/differential_function_clean` (DRAM KEGG step form, 399 modules, present at
>= 2/3 of top-level steps); the raw steps fraction and the fraction of module KOs are carried on
every row, because the single-step respiratory complexes reduce the 2/3 rule to "any one subunit".
Outputs: `results/v2_thermo_metabolism_summary.tsv` (6,526 rows, feature x group),
`results/v2_thermo_metabolism_by_rep.tsv` (5,438 reps x 403 columns),
`results/v2_thermo_metabolism_report.txt`.

Comparators are 1:1 nearest-neighbour matched on completeness without replacement, caliper 1.0
points: Halobacteriota 181 of 237 pairs (81.87 vs 81.87), other archaea (the 354 = all archaea
minus the 237) 216 pairs (84.14 vs 84.14), bacteria 237 of 237 (85.42 vs 85.42).

**Encoded lifestyle.** An aerobic, archaeal, peptide-and-polymer heterotroph. EMP glycolysis 82.3 %,
gluconeogenesis 81.0 %, archaeal non-oxidative PPP 82.7 % (the oxidative phase is absent, 0.4 %),
no Entner-Doudoroff of any variant (0 %), TCA 83.1 % with 2-oxoglutarate:ferredoxin oxidoreductase
96.2 % rather than the E1/E2 dehydrogenase (3.0 %), glyoxylate cycle 79.3 %. Respiration: complex I
70.9 % at >= 7 of 14 nuo subunits, complex II 77.6 %, coxA+coxB 58.6 % (coxA alone 70.9 %), bd
oxidase 16.0 %, V/A-type ATPase 91.1 % and F-type 0 %; a canonical bc1 is absent (7.2 % at two
subunits) but a Rieske subunit K03886 sits in 38.8 % against 0.6 % of matched Halobacteriota.
Uptake is peptides and amino acids, not sugars: branched-chain amino-acid ABC 78.9 %, peptide/nickel
ABC 54.9 %, multiple-sugar ABC 2.1 % against 72.9 % in matched Halobacteriota. Median 28 peptidase
KOs and 5 glycoside-hydrolase KOs per genome; endoglucanase 79.3 %, alpha-amylase 56.5 %, chitinase
44.3 %, extracellular serine protease AprX 95.8 %. Nitrogen runs through glutamate dehydrogenase
94.5 % plus glutamine synthetase 82.3 %, not GOGAT (19.4 %); no urease, no nitrate reduction, no
nitrogen fixation. Sulfur: assimilatory sulfate reduction 10.1 %, PAPS reductase 15.6 %, but
sulfide:quinone oxidoreductase 34.2 % against 1.7 % in matched Halobacteriota. Osmolytes are thin
(trehalose otsAB 7.2 %, ectoine 0.4 %), archaellum 99.6 %, no bacterial flagellum or type IV pilus,
CRISPR-Cas 79.3 % with type I cas3 33.8 %.

**Auxotrophies.** Median 11 of 16 module-testable amino acids and 5 of 17 cofactors. The
branched-chain block is the sharp one: the ilvBN/ilvC/ilvD carbon skeleton is in 1.7 % against
91.2 % of matched Halobacteriota, while the terminal transaminase K00826 is in 93.2 %, so valine and
isoleucine are taken up, not made. Pantothenate, pyridoxal, thiamine (0.4 %), cobalamin (0.4 %) and
F420 (0 %) are absent; folate 33.3 %, riboflavin 81.0 %, heme 61.6 %.

**Classes.** SW-10-69-26 (185) and E2 (51) are at the same completeness (85.64 vs 82.45, MWU
p = 0.18) but not the same corrected size (3.84 vs 2.89 Mb, p = 2.85e-19), so the gaps are reduction,
not recovery: at 51 of 51 completeness-matched pairs E2 carries 819 KOs against 952. E2 drops the
importers (peptide/nickel 3.9 vs 68.6 %, branched-chain 35.3 vs 90.8 %), pyruvate kinase
(35.3 vs 88.1 %) and tryptophan biosynthesis (45.1 vs 94.1 %), and is the class that keeps GOGAT
(68.6 vs 5.4 %).

**Caveats.** Everything here is KO-model sensitivity read forwards instead of backwards: in a phylum
that is 31 % dark a missing KO can be a missing model, so absences are floors and presences are
safe. KEGG has no biosynthesis module for Ala, Asp, Asn, Glu, Gln or Gly, and its Phe and Tyr
modules score 0 in every archaeal group here, so the prototrophy count is out of the 16 testable
amino acids and is not an auxotrophy call. Single-step modules (complexes I to IV, ATPases,
ribosome) are present on one annotated subunit, which is why the strict subunit-count features are
reported beside them. Draft genome-scale reconstructions were looked for and not used: gapseq 1.4
output exists at `05_metabolism/gapseq_out/` (3,835 per-MAG directories) but on the ambiguous
`<sample>__SemiBin_<n>` basis, and it covers only 1,037 of the 5,438 representatives, 24 of the 237
and 11 of the 269 Halobacteriota, so no comparison can be made from it. No Bakta annotation exists
for the representative set, so transporter and CAZyme-like counts are KO-based (EC 3.2.1 and EC 3.4
classes, transport-system definitions) and are floors.

## STATUS (2026-09-24, NOT REGISTERED): Thermoplasmatota known relatives (ANI, AAI)

Question: how many of the 237 have any known relative, and how close. Comparators: all 218 GTDB r232
species reps in c__SW-10-69-26 (21) and c__E2 (197), accessions from the r232 ar53 reps taxonomy,
FASTAs from the NCBI FTP, 218/218 md5-verified; all 83 MGnify genomes in those classes (genomes API
lineage filter, 10 + 73; no published md5, local md5 kept). **54 of the 83 duplicate a GTDB rep**
(41 the same NCBI assembly, 13 more at skani ANI >= 99.9 % with AF >= 90 %); they are dropped from
every nearest call and lend their MGnify biome to the rep they duplicate, leaving 29 distinct MGnify
comparators. **SMAG skipped**: Zenodo 8223844 holds one archive split into ~370 x 100 MB parts,
no per-genome files. 538 genomes in total.

Methods. skani 0.3.1 `triangle --sparse -s 75` (same binary as `cultured_comparison_clean`).
No CompareM or EzAAI in any env, so RBH AAI is computed here (`v2_thermo_relatives_aai.py`):
Pyrodigal 3.7.1 (darkfunc env, single mode) on every genome alike (1,366,437 proteins), DIAMOND
2.1.16 blastp all-vs-all `--very-sensitive -e 1e-5 --id 30 --query-cover 70 --subject-cover 70`,
best hit per (protein, target genome), reciprocal pairs, AAI = mean identity; nearest calls need
>= 100 RBH (144,452 of 144,453 pairs pass). GTDB-Tk 2.7.2 classify_wf on the 83 MGnify genomes
against `/ibex/reference/KSL/gtdb-tk/release232` (VERSION_DATA=r232 checked in the script).

| of the 237 | all | SW-10-69-26 (185) | E2 (51) | Thermoplasmata (1) |
|---|---|---|---|---|
| any skani hit to a comparator (ANI >= ~80) | 3 | 3 | 0 | 0 |
| ANI >= 95 % (species) | **0** | 0 | 0 | 0 |
| AAI >= 65 % (genus-level heuristic) | **111** | 87 | 24 | 0 |
| AAI >= 68.33 % (same-family q95, stricter genus cut) | **35** | 23 | 12 | 0 |
| best AAI median (IQR) | 64.77 (58.76 to 66.39) | 64.79 | 64.19 | 48.68 |
| best AAI max | 80.04 | 80.04 | 76.96 | 48.68 |
| nearest by AAI is a distinct MGnify genome | 8 | 0 | 8 | 0 |

The Thermoplasmata rep (o__RBG-16-68-12) has no same-class comparator in this set, so its 48.68 is a
between-class value, not evidence of no relative.

Best-AAI bins (all 237): <50 1 · 50-55 6 · 55-60 62 · 60-65 57 · 65-70 87 · 70-80 23 · 80-95 1.
The three skani hits are ANI 84.2 to 85.3 at AF 14 to 18 %. GTDB-Tk's own r232 ANI screen gave no
closest reference for any of the 237, which agrees. Median AAI to the nearest *other desert* rep is
73.1, above the best external AAI: the radiation is closer to itself than to anything published.
The closest pairs are three o__SW-10-69-26 reps placed in g__CAYXGD01 at AAI 79.6 to 80.0 to
CAYXGD01 sp964464425 / sp964492585 (the same reps give the skani hits); nothing reaches species.

- Type genome GCA_003009755.1 (Atacama halite) is the nearest comparator of only **4 of 185**
  SW-10-69-26 reps; AAI to it median 47.3 (45.1 to 74.9). 153 of the 185 sit in o__JACQPN01 and 2 have no order; of the
  30 in o__SW-10-69-26, median AAI to it 58.4. The SW reps' nearest comparators are mostly
  DATLGO01 sp035540375 and JAJPHT01 species.
- GTDB r232 reps: for 28 of 218 (16 of 21 SW-10-69-26, 12 of 197 E2) the nearest species rep
  (desert or GTDB) is a desert rep; 19 at AAI >= 65; 40 have some desert rep at AAI >= 65.
- GTDB reps block uses desert + GTDB reps only, so MGnify duplicates cannot double-count there.
- MGnify at r232: 73 E2 and 10 SW-10-69-26, no class change. Biomes: marine sediment 60, soil 19,
  marine 4. Via MGnify twins, the nearest comparator of **142 of 237** desert reps is a soil genome
  (3 marine sediment, 92 GTDB reps with no MGnify biome).
- Calibration of the 65 % line on this pipeline (GTDB rep pairs by lowest shared rank): same genus
  median 75.4 (q05 63.1), same family 61.8 (21 % >= 65), same order 55.7, same class 47.0. So 65 is
  a lenient genus cut here; read it as "genus or close family", never as a genus call.
- Not a contradiction with the novelty ordinal above: GTDB-Tk places 228 of the 237 in a family and
  102 in a genus, but those are placeholder names, which `first_unnamed_rank` counts as unnamed.
  227 of 228 have their nearest-by-AAI comparator in the same r232 family.

Outputs: `results/v2_thermo_relatives_pairs.tsv` (237 rows), `_summary.tsv`, `_report.txt`,
`_comparators.tsv` (manifest: accession, source, catalogue, biome, lineage, bytes, md5, md5 check),
`_skani.tsv.gz`, `_aai.tsv.gz` (all 144,453 pairs), `_gtdbtk_mgnify.tsv`. The summary script reads
only these, so it reruns without the downloads. **Deleted**: `tmp/thermo_relatives/` (5.0 G: 301
comparator FASTAs, proteins, DIAMOND db and hit arrays, GTDB-Tk work dir).
