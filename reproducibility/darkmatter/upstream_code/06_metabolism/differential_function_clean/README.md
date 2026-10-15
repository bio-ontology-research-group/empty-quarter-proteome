# differential_function_clean: which functions differ between soil compartments

Owner: velazqal. Promoted 2026-09-20 from the exploratory folder
`alelopezv/eqmags_poc/04_differential_function`. This folder is authoritative; the
exploratory copy is superseded.

Carriage = 100 x sum(rel. abundance x presence) / sum(rel. abundance), the percent of the
MAPPED community that encodes a feature. Effects are differences in percentage points of
carriage, paired within a site visit and then collapsed to one value per site.

## STATUS (2026-09-20): full chain rerun, every exploratory number reproduced

| step | script | how it ran | wall | peak RSS |
|---|---|---|---|---|
| feature universes | `00_build_features.py` | login | 6.5 s | 1.0 G |
| the scan | `10_differential.py` | login | 1:22.8 | 1.4 G |
| figures | `20_figures.py` | login | 7.9 s | 0.3 G |
| known-trait sanity | `30_sanity_known_traits.py` | login | 11.3 s | 1.0 G |
| review diagnostics | `40_review_diagnostics.py` | login | 20.9 s | 1.4 G |

No Slurm job was needed; the whole chain is under two minutes on a login node. Python
`/home/velazqal/miniconda3/envs/python-utils/bin/python`. Logs in `logs/`. Note that
`/ibex/scratch/projects/c2014` and `/ibex/project/c2014` are the same directory, so the
scripts print the second spelling.

**Reproduction.** Run first with the paths repointed and nothing else changed: all seven
outputs (`feature_definitions.json`, the three `*_tests.tsv`, `passing_modules.tsv`,
`carriage_passing.tsv`, `run_log.json`, `known_traits_sanity.tsv`) and both PNGs came out
byte-identical to the exploratory copies. The audit columns below were then added and the
chain rerun: `run_log.json` is still byte-identical and every column shared with the first
run is equal value for value. **No reported number changed on promotion.**

## What changed, and why

| change | file | effect on results |
|---|---|---|
| output root located from the script path instead of a hard-coded `eqmags_poc` constant | `common.py` | none |
| pairing helpers (`site_compartment_trip`, `compartment_diffs`, `run_compartment`) moved from `10_differential.py` into `common.py` | both | none; the diagnostics now test the same code the results come from instead of a copy |
| raw `p` written next to `q` | all three test tables | none, additive |
| `n_steps`, `steps_required`, `n_kos`, `prev_species` written per feature | all three test tables | none, additive; makes the presence rule auditable row by row |
| `median_a`, `median_b` (the two carriage levels being differenced) | compartment table | none, additive; an effect in pp is unreadable without its baseline |
| `n_sites_Actinomycetota`, `n_sites_Pseudomonadota` | compartment table | none, additive; the phylum series is blanked where the phylum holds under 1 % of mapped abundance, so it needs its own n |
| compartment-only columns (`comp_*`) added beside the any-axis "best" columns | `30_sanity_known_traits.py` | none, additive; the old `best_*` columns pick the largest effect on ANY axis, which for two KOs is a non-passing trip row |
| docstring corrected: the module rule is 2/3 of STEPS, not of the raw KO list | `10_differential.py` | none; the code always did steps, the docstring said KOs |
| new script `40_review_diagnostics.py` | new | writes only `results/review_diagnostics.txt` |

Nothing was rerun differently, because no defect found below changes a number that the
paper would report or makes a reported claim false.

## Inputs (read-only; `W = /ibex/scratch/projects/c2014/rob/metagenomics`)

| input | size |
|---|---|
| `W/10_thermal/abundance/rep_relative_abundance.tsv` | 5,438 reps x 550 samples |
| `W/04_functional/repko_clean/results/rep_ko_v2.tsv` | 5,438 genomes, 11,395 distinct KOs |
| `W/databases/dram/module_step_form.20260623.tsv` | 399 modules, 3,288 module x KO x path rows |
| `W/databases/dram/kofam_ko_list.tsv` | KO names |
| `W/metadata/sample_metadata.tsv` | 551 samples, site x compartment x trip |
| `W/05_refine/_drep/out/data_tables/genomeInformation.csv` | completeness, length |
| `W/05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv` | phylum for the control |
| `R/results/darkmatter_ecology/env_matrix.tsv` (paper repo, read-only) | XRF and climate; chemistry axis only |

## Method review

Verdicts. "Changes a number" means a number the paper would report, i.e. the compartment
axis only. Evidence for every line is in `results/review_diagnostics.txt`.

| point | defect | changes a reported number | action |
|---|---|---|---|
| a. 2/3-of-steps presence rule and completeness bias | partly | no | audit columns added; sensitivity run |
| b. site-visit pairing, replicate collapse | no | no | sensitivity run on the inhomogeneous arm |
| c. phylum control is sign only | yes, as reported | no | magnitudes were already recorded and are now reported beside every effect |
| d. chemistry confounded with trip | yes | no | axis not reported; confound quantified |
| e. trip axis, 8 pairs | yes | no | axis not reported |
| f. BH across all features per contrast | no | no | pooling the two contrasts costs 0 rows; families are not mixed |
| g. size normalisation | no | no | checked, holds |
| h. mapped-community denominator | yes, for root vs bulk | no | quantified; stated as a limitation |

**a. The presence rule.** A module is present when at least 2/3 of its DRAM top-level
steps are satisfied, a step being satisfied by any one of its alternative KOs. Two real
problems. First, the rule is not comparable across modules: of the 282 tested, 58 have a
single step, so "2/3 of steps" becomes "any one of up to 79 alternative KOs", while an
eight-step module needs six. Four of the modules passing each compartment contrast are
one-step (root: M00144, M00145, M00152, M00153; surface: M00133, M00145, M00152, M00153),
and five one-step modules are a single KO that is also tested in the KO family, so the
identical test is counted in both families under different q. M00133 (polyamine
biosynthesis) is exactly K00797: +8.730 pp on surface minus deep in both tables. Second,
presence is read off MAGs that are 50 to 100 % complete, so carriage is biased low
everywhere; that only distorts a contrast if the bias differs between arms, and it does,
slightly. The abundance-weighted completeness of the mapped community is 0.66 points
higher in root-adjacent than in bulk soil (p = 0.0016) and not different between surface
and subsurface (-0.37, p = 0.19). Rerunning the whole scan on the 2,675 reps that are at
least 90 % complete reverses no sign at all (26/26, 21/21 modules; 716/716, 734/734 KOs)
but moves magnitudes: 21 of 26 root modules and 14 of 21 surface modules still clear the
5 pp floor, the drop-outs being the low-amplitude ones (M00020, M00121, M00167, M00793,
M00868 on root; M00024, M00025, M00063, M00093, M00549, M00741, M00793 on surface). Not a
defect that changes a number, but the module label for a one-step module overstates what
was tested, and effects within about 2 pp of the floor are not robust to the genome set.

**b. Pairing and replicates.** Sound. Replicates inside a site x compartment x visit unit
(316 units with one, 61 with two, 37 with three) are collapsed by median, differences are
taken within a visit so a trip1 surface is never compared with a trip3 subsurface, and the
per-site value is the median over visits, so the test is over sites (52 and 37) and not
over samples. One inhomogeneity: at the 128 root-vs-bulk visits, bulk means the median of
surface and subsurface at 64 visits, subsurface alone at 43 and surface alone at 21.
Restricting to the 64 two-sided visits leaves every sign intact and the magnitudes
essentially unchanged (median ratio 0.99), and 19 of 26 modules still pass, the loss being
power from 52 sites down to 36. Not a defect.

**c. The phylum control.** The code does record the within-phylum magnitudes; they are the
`eff_Actinomycetota` and `eff_Pseudomonadota` columns, on the same sites as the whole
community (median 52 and 52 against 52). Reporting them changes what can be claimed, and
it should. Across the 23 phylum-passing module rows the smaller of the two within-phylum
effects is a median 0.16 of the whole-community effect, 13 of the 23 sit below 0.2, and in
only 3 rows do both phyla clear the 5 pp floor on their own. (The exploratory README's
"16 of 27" is the same statistic over the 27 deduplicated module rows of all three axes;
both reproduce.) The extreme cases are
M00025 and M00024 (tyrosine and phenylalanine biosynthesis), which "pass" on an
Actinomycetota effect of +0.0047 pp against a whole-community +6.0. So the honest claim is
that the sign of the compartment difference is not created by the Actinomycetota to
Pseudomonadota ratio, not that the difference holds at comparable strength inside a
phylum. For KOs the picture is the same (median ratio 0.25, 344 of 813 rows below 0.2).

**d. Chemistry.** Confounded, and the axis is not reported. XRF coverage is 32 % of trip1,
31 % of trip3, 53 % of trip5 and 100 % of trip4 samples, and a trip is an extraction and
sequencing batch, so an element correlation on n = 249 cannot be separated from batch.
Samples are also pseudoreplicated within site, only latitude and longitude are controlled.

**e. Trip.** Eight of the ten trip pairs met the 15-unit minimum; BH is applied inside each
pair, so the eight pairs are eight separate families with no correction between them, and
the no-reversal rule that follows is a filter, not a test. Not reported.

**f. Multiple testing.** BH runs over all features of one basis within one contrast, so
modules (282 tests) and KOs (4,127) are never mixed into one family; they are simply not
comparable to each other, which is a reporting caveat, not a bug. The uncorrected
multiplicity is across the two compartment contrasts. Pooling them into one family per
basis costs nothing: 47 module rows and 1,450 KO rows pass either way, and the distinct
feature counts (32 and 1,070) are identical. The one real duplication is the five one-step
modules that are also KOs (point a).

**g. Size normalisation.** Every effect is recomputed with each genome weighted by 1/size
in Mb (completeness-corrected) and the sign compared. It holds everywhere: 47 of 47 module
rows and 1,219 of 1,219 KO rows in `passing_modules.tsv`, and 26/26, 21/21, 716/716,
734/734 on the compartment contrasts. No defect.

**h. The denominator.** Carriage is a percent of the mapped community and the mapped
fraction is itself a per-sample quantity: 50.9 % in root-adjacent samples against 45.3 %
subsurface and 44.3 % surface, a paired +5.58 pp for root minus bulk (p = 1.3e-05) and a
null -1.19 pp for surface minus subsurface (p = 0.36). Normalising by the mapped sum is
the right treatment, and it is what the code does, but it cannot fix the interpretation:
the invisible half of a root-adjacent sample is a different half from the invisible half
of a bulk sample, so every root-vs-bulk number is a statement about the mapped community
only. Related: the mapped community of root-adjacent and of surface samples carries more
KOs per genome (+48.6 and +28.7 KOs, p = 1.6e-04 and 2.3e-03), which is why the KO effect
distribution sits at +0.52 pp rather than at zero. That offset is a constant per contrast
and does not touch the module results (median effect exactly 0.000, and centring the
distribution moves no module in or out), but it does sit under the KO list: 102 and 120 of
the passing KOs are within half a point of the floor, and applying the floor to centred
effects would swap about 9 % of the list (65 out, 54 in on root; 95 out, 36 in on
surface). The module list is the safe unit to report; the KO list should be described as a
ranking rather than a membership.

## What the paper can report: the compartment axis

Stages, per contrast (rows, not deduplicated features):

| basis | contrast | tested | q < 0.05 | + 5 pp floor | + phylum sign | partial | fail |
|---|---|---|---|---|---|---|---|
| module | root - bulk | 282 | 145 | 26 | 14 | 3 | 9 |
| module | surface - subsurface | 282 | 119 | 21 | 9 | 1 | 11 |
| KO | root - bulk | 4,127 | 2,496 | 716 | 441 | 68 | 207 |
| KO | surface - subsurface | 4,127 | 2,458 | 734 | 372 | 54 | 308 |

Deduplicated over the two contrasts (the headline counts): 32 modules and 1,070 KOs clear
floor and FDR, 20 modules and 626 KOs also clear the phylum control on their stronger
contrast (22 and 688 on at least one contrast).

Modules passing floor, FDR and the phylum sign control, root-adjacent minus bulk, 52
paired sites, all with 52 sites in both control phyla:

| module | name | effect pp | q | Actino | Pseudo | root % | bulk % | steps |
|---|---|---|---|---|---|---|---|---|
| M00153 | Cytochrome bd ubiquinol oxidase | +13.04 | 2.8e-05 | +17.02 | +1.88 | 62.3 | 48.3 | 1 of 1 |
| M00021 | Cysteine biosynthesis, serine => cysteine | +10.99 | 4.4e-04 | +8.35 | +3.44 | 50.5 | 38.8 | 2 of 2 |
| M00145 | NAD(P)H:quinone oxidoreductase, chloroplasts and cyanobacteria | -10.25 | 1.4e-04 | -8.75 | -0.97 | 35.6 | 46.8 | 1 of 1 |
| M00378 | F420 biosynthesis | -10.10 | 1.5e-04 | -3.77 | -0.41 | 44.0 | 55.5 | 2 of 2 |
| M00008 | Entner-Doudoroff pathway, glucose-6P => G3P + pyruvate | +9.36 | 1.0e-04 | +14.57 | +5.75 | 44.9 | 36.9 | 3 of 4 |
| M00167 | Reductive pentose phosphate cycle, G3P => Ru5P | +8.52 | 2.9e-04 | +7.30 | +0.89 | 72.1 | 63.7 | 5 of 7 |
| M00549 | Nucleotide sugar biosynthesis, glucose => UDP-glucose | +7.94 | 6.6e-04 | +11.23 | +0.90 | 79.8 | 73.4 | 2 of 3 |
| M00144 | NADH:quinone oxidoreductase, prokaryotes | -7.83 | 1.6e-05 | -11.80 | -0.03 | 86.0 | 96.9 | 1 of 1 |
| M00053 | Pyrimidine deoxyribonucleotide biosynthesis | +7.59 | 3.1e-05 | +10.87 | +5.87 | 59.9 | 49.5 | 6 of 8 |
| M00020 | Serine biosynthesis, glycerate-3P => serine | +7.34 | 7.2e-04 | +10.47 | +1.16 | 80.2 | 73.2 | 2 of 3 |
| M00868 | Heme biosynthesis, animals and fungi, glycine => heme | +7.18 | 2.7e-04 | +8.15 | +0.92 | 55.9 | 48.5 | 6 of 8 |
| M00121 | Heme biosynthesis, plants and bacteria, glutamate => heme | +6.93 | 1.4e-04 | +11.02 | +2.48 | 74.6 | 65.9 | 7 of 10 |
| M00118 | Glutathione biosynthesis, glutamate => glutathione | +6.56 | 7.0e-04 | +0.12 | +2.48 | 24.2 | 17.0 | 2 of 2 |
| M00168 | CAM (Crassulacean acid metabolism), dark | -5.48 | 1.1e-04 | -4.82 | -3.29 | 44.1 | 49.9 | 2 of 2 |

Surface minus subsurface, 37 paired sites, 36 sites in both control phyla:

| module | name | effect pp | q | Actino | Pseudo | surf % | sub % | steps |
|---|---|---|---|---|---|---|---|---|
| M00378 | F420 biosynthesis | -10.07 | 6.7e-04 | -1.68 | -3.07 | 48.7 | 60.7 | 2 of 2 |
| M00133 | Polyamine biosynthesis, arginine => putrescine => spermidine | +8.73 | 3.7e-04 | +3.50 | +2.52 | 50.2 | 40.6 | 1 of 1 |
| M00345 | Formaldehyde assimilation, ribulose monophosphate pathway | +8.69 | 8.7e-05 | +9.12 | +4.57 | 30.2 | 18.2 | 2 of 3 |
| M00579 | Phosphate acetyltransferase-acetate kinase, acetyl-CoA => acetate | +7.01 | 2.7e-04 | +5.68 | +0.84 | 16.1 | 10.3 | 2 of 2 |
| M00793 | dTDP-L-rhamnose biosynthesis | -6.18 | 9.7e-04 | -5.53 | -6.82 | 62.0 | 68.6 | 3 of 4 |
| M00025 | Tyrosine biosynthesis, chorismate => tyrosine | +5.99 | 8.7e-05 | +0.005 | +10.21 | 14.6 | 7.2 | 2 of 2 |
| M00741 | Propanoyl-CoA metabolism, propanoyl-CoA => succinyl-CoA | -5.97 | 2.7e-04 | -2.85 | -6.51 | 75.8 | 84.9 | 2 of 3 |
| M00024 | Phenylalanine biosynthesis, chorismate => phenylalanine | +5.74 | 8.7e-05 | +0.005 | +9.66 | 13.6 | 6.9 | 2 of 2 |
| M00017 | Methionine biosynthesis, aspartate => homoserine => methionine | +5.14 | 2.5e-03 | +5.49 | +1.76 | 70.1 | 64.0 | 5 of 7 |

The KO lists are in `results/compartment_tests.tsv` (filter `passes` and
`within_phylum == pass`). Strongest five per contrast: root-adjacent K07552 MFS multidrug
transporter +21.9, K03439 tRNA guanine-N7 methyltransferase +20.8, K04047
starvation-inducible DNA-binding protein +19.4, K01849 methylmalonyl-CoA mutase -18.9,
K22205 S-adenosyl-L-methionine hydrolase -17.6; surface K03321 sulfate permease +19.2,
K01624 fructose-bisphosphate aldolase II +17.5, K03781 catalase +17.0, K06916 ZapE +15.7,
K09781 uncharacterized +15.4.

### Known traits reappear (compartment axis, `results/known_traits_sanity.tsv`)

| trait | KO | contrast | effect pp | q | Actino | Pseudo | phylum |
|---|---|---|---|---|---|---|---|
| photolyase | K01669 | surface - sub | +8.80 | 1.3e-03 | +12.82 | +12.75 | pass |
| spore photoproduct lyase | K03716 | root - bulk | -5.49 | 8.0e-06 | -5.44 | -6.42 | pass |
| betaine transport | K02168 | surface - sub | +10.44 | 1.7e-04 | +7.13 | +6.70 | pass |
| betaine transport | K05020 | surface - sub | +10.29 | 1.7e-04 | +16.07 | +6.89 | pass |
| betaine/proline permease | K02001 | root - bulk | +5.16 | 1.4e-04 | +4.76 | +7.44 | pass |
| trehalose/maltose transport | K10237 | root - bulk | -5.34 | 1.6e-05 | -3.84 | -2.46 | pass |
| hydrogenase large subunit | K06281 | root - bulk | -15.64 | 8.0e-06 | -19.78 | -3.63 | pass |
| hydrogenase small subunit | K06282 | root - bulk | -13.88 | 9.0e-06 | -16.95 | -2.93 | pass |
| CO dehydrogenase small subunit | K03518 | surface - sub | -5.43 | 1.1e-02 | -1.15 | -10.03 | pass |

The Hyp maturation set (K04653, K04654, K04655, K04656, K03605) follows the hydrogenase
subunits at -12 to -13.4 pp on root minus bulk, all phylum pass. Trehalose synthase K13057
(-8.8 pp) and choline dehydrogenase K00108 (+8.2 pp) clear floor and FDR but fail the
phylum control. BGC and AMR have no KEGG representation and are not testable on this
basis.

## Reconciliation against the exploratory README

| headline | exploratory | rerun | difference |
|---|---|---|---|
| compartment modules floor+FDR / + phylum | 32 / 20 (+3 partial) | 32 / 20 (+3) | none |
| compartment KOs floor+FDR / + phylum | 1,070 / 626 (+85) | 1,070 / 626 (+85) | none |
| trip modules / KOs | 2 / 1 and 50 / 32 | same | none |
| chemistry modules / KOs | 13 / 6 and 99 / 74 | same | none |
| any axis | 46 of 282 modules, 26 phylum pass; 1,168 of 4,127 KOs, 705 | same | none |
| trip before the no-reversal rule | 9 modules, 371 KOs | same | none |
| size-normalisation sign check | 47/47 and 1,219/1,219 | same | none |
| samples, species, mean mapped | 549, 5,438, 47.83 % | same | none |
| smaller within-phylum effect, median share, 16 of 27 module rows below 0.2 | 0.16 | 0.16, 16 of 27 | none |

`run_log.json` is byte-identical, so every count above is the same object, not a
recomputed match. Two presentation points that are not differences: the headline 32 and 20
are feature counts deduplicated over the two contrasts keeping the stronger one (15 of the
32 modules pass both contrasts, and 22 rather than 20 modules pass the phylum control on
at least one contrast); and the exploratory README's sanity paragraph gives CO
dehydrogenase K03518 on a trip pair, whereas its compartment result, which is what the
paper reports, is -5.43 pp on surface minus subsurface, q = 0.011, phylum pass.

## Denominators

| quantity | value |
|---|---|
| samples | 549 of 550 (F5Dr1 dropped, zero mapped abundance) |
| species representatives | 5,438 |
| sites | 61; 52 in root minus bulk, 37 in surface minus subsurface |
| site x compartment x visit units | 414 (316 with 1 replicate, 61 with 2, 37 with 3) |
| paired site visits | 128 root vs bulk, 69 surface vs subsurface |
| modules tested | 282 of 399 (115 absent everywhere, 2 universal) |
| KOs tested | 4,127 of 11,395 (5 to 95 % prevalence across reps) |
| mapped fraction | mean 47.83 %; PlantRoot 50.9, Deep 45.3, Surface 44.3 |

## Caveats

Mapped fraction only, and the unmapped half is not a random half, nor the same half in
each compartment (point h). Capacity, not activity. Presence is read off 50 to 100 %
complete MAGs, so an absent module is not evidence of absent biology: three of the four
methanogenesis modules (M00357, M00563, M00567) fall to zero prevalence and the fourth
(M00356) reaches 0.22 % of reps, although 269 Halobacteriota MAGs are in the catalogue. Module definitions are DRAM's step form, not KEGG direct, and KEGG names some
modules for eukaryotic contexts (M00145 chloroplasts, M00168 CAM, M00868 animals and
fungi); those rows are prokaryotic enzyme sets and must not be reported under their KEGG
title. The phylum control is sign-only by design; magnitudes are in the table and usually
much smaller. Module and KO BH families are not comparable. The chemistry and trip axes
are computed and stored but are confounded with batch and are not reportable.
