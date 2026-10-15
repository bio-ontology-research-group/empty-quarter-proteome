# Dossier 05 — Environment and stress (§2.5 dark-family ecology · §2.3 stress censuses)

Read: `paper/main.tex` (§2.3, §2.5, M.1, M.5, M.7, M.10, M.13, fig:environment, fig:stress), `paper/supplementary.tex` (S3/S4/S6), `notes/CATALOGUE_NUMBERS.md`, `notes/CORRECTIONS.md`, and the scripts+results of `M/06_darkmatter/darkfam_ecology_clean/`, `M/06_metabolism/uv_survival_clean/`, `M/06_metabolism/osmolyte_clean/` (M = `/ibex/scratch/projects/c2014/rob/metagenomics`, read-only).

## Questions (both)

- §2.5: is the 106,042-family dark repertoire ecologically organised or spread indifferently — and if organised, by chemistry or only by distance?
- §2.3: does the catalogue carry the survival toolkit a hyperarid sand sea requires, and is it organised by light (within site, across compartments) or by salt (between sites)?
- Shared payload: known and dark biology are organised by the same master gradient (salt) — the argument that the dark structure is biology, not artefact.

## Basis

| | §2.5 | §2.3 |
|---|---|---|
| object | 106,042 dark families × 551 samples (member-aware CoverM TPM) | 5,438 species reps × 548 profiled samples (rep-level CoverM) |
| env-annotated | 547 (`results/darkmatter_ecology/env_matrix.tsv`, verified 547 rows) | same table ∩ the 548-sample profile cohort |
| ordination / test set | 105,977 families at prevalence ≥5; 94,568 tested on XRF axes, 100,293 on aridity | 28-KO damage panel + 13-KO solute panel |
| XRF / climate cohort | 251 / 456 | 248 / 453 |
| spatial design | lat/lon partialling only | 44 sites with both Surface and Deep profiled |

**Why the denominators differ — two reasons, both real.** (i) 251/456 vs 248/453: §2.5 runs on the 551-sample member-aware cohort, §2.3 on the 548-sample rep-level cohort. M.7 states this ("taken here over the 551-sample member-aware cohort rather than M.10's 548-sample cohort"). Difference is exactly 3 on both axes — consistent; the paper never says which 3 (INFERENCE: arithmetic, no verified sample list). (ii) 94,568 vs 100,293: `MIN_PRESENT = 25` in `eco_1_indicators.py` is applied *inside each axis subset* — the same absolute ≥25-sample gate against 251 XRF samples vs 456 climate samples, so more families clear it on climate. Neither equals the 105,977 used for the ordination (prevalence ≥5 over 547). JUDGMENT: defensible, but three denominators in one section with no explanation is a referee magnet; one clause fixes it.

## Environmental variables

`code/build_env_matrix.py` is the only definition that exists anywhere. 547 rows × 45 columns.

- **13 XRF elements tested per family**: V, Mo, Ni, Fe, Cu, Co, Mn, Zn, S, Cl, Na, Ca, K. The matrix holds 33; "oxide duplicates" dropped "to limit collinearity".
- **The PERMANOVA geochemistry block is a different, 8-axis set**: V, Fe, Ca, K, Si, Ti, Mn, salinity_proxy (`eco_3_community.py:120`). **S, Cl and Na — the axes that carry the per-family signal — are not in the multivariate block**, while Si, Ti and V, which carry none, are. M.7's "eight geochemical axes" is literally true and substantively misleading.
- **salinity_proxy = Na + Cl + S**, an unweighted sum of three element columns (`build_env_matrix.py:64-66`). Defined nowhere in the paper. This is the composite §2.5 itself calls spatially confounded.
- **aridity_index = precip/(temp+10)** — the de Martonne index, computed in-repo from `AnnualMeanTemp`/`AnnualTotalPrecip` in `data/geodata/trip<N>_geodata.tsv`. **Not** a UNEP P/PET index. M.1's todo says its source "still needs naming" — the formula is recoverable from the builder, the climate series is not.
- **Units are self-contradicted in the repo.** The builder comments call the values "elemental %"; M.1 quotes the companion descriptor saying they "carry no documented unit assertion and must not be read as percentages".
- **XRF workflow**: the builder reads `xrf_lab_table_all_trips.tsv` (725 rows) = laboratory XRF across all five trips, not the trip-5 handheld (INFERENCE from the builder comment; the file itself is unreachable). This answers most of M.1's todo.
- **pH does not exist** — not a column, not an element, nowhere in the 45 columns. A salt-and-aridity soil paper will be asked.
- **⚠ The env matrix is not rebuildable.** Its sources (`/home/leechuck/Documents/...`, `/home/leechuck/Public/...`) do not exist on this filesystem — verified. README: "Built 2026-06-28 from rob's local sources (not rebuildable on IBEX); md5 4d6eafad897818075a0fc6ae32bdf754". The shipped TSV is the only surviving artefact behind every environmental number in §2.3, §2.4 and §2.5.

## How it ran

**§2.5** (`eco_1/2/3`, job 50575481, 3m53s). All statistics hand-written in `eco_stats.py` — no vegan, no scikit-bio. McArdle–Anderson db-PERMANOVA on the Gower-centred matrix, hat-matrix sums of squares, 999 permutations of the Gower matrix, `seed=1`; Hellinger then Bray–Curtis. Blocks: compartment (drop-first dummy), space (z-scored lat/lon), climate (temp, precip, aridity; n=456) also |space, geochemistry (8 axes; n=251) also |space, single-axis vanadium. Reported R² is `R2_marginal`; for geochem|space the partial R² is 0.0881 against the 0.0802 printed as 8.0%. Mantel is raw distance-decay (Bray–Curtis vs haversine, 999 perms), not partial. Per-family: Spearman by rank-transform with a t approximation, BH within axis; partial Spearman by rank + hat-matrix residualisation on ranked lat/lon (df = n−4), five master axes only.

**⚠ Borcard partitioning, as run, partitions only compartment and space.** The climate block enters the varpart call only `if len(clim_idx) == len(samples)` — 456 ≠ 547, so it was silently dropped; geochemistry was never a candidate block (`eco_3_community.py:150-156`). `varpart.tsv` has two rows and reproduces the PERMANOVA. M.7 and M.13 both cite it as if it partitioned the environmental blocks. This is the largest methods-vs-code gap in either section.

**§2.3** (`uv_1..6`, `osmo_1..7`, all local, no Slurm). 28-KO damage panel across 8 groups (K01859 excluded — mislabelled chalcone isomerase, zero effect) and 13-KO solute panel against `repko_clean/results/rep_ko_v2.tsv`; denominators are the 5,438 representatives. Three community metrics (M.5): occurrence, whole-community (unmapped = non-carriers, a lower bound), mapped-community. Replicates collapse to the median within each site×compartment cell; Surface vs Deep by Wilcoxon signed-rank across the 44 paired sites, Holm within metric; Kruskal–Wallis across three compartments as secondary. Gradients: partial Spearman on ranked lat/lon, BH within the salt block (n=248) and climate block (n=453). Both census scripts now read `_gtdbtk_r232` (`uv_1_census.py:28`, `osmo_1_census.py:27`) — **`PIPELINE_GUIDE.md`'s r226 warning is stale**, superseded by the 2026-08-17 fix.

## Controls: CRITICAL

**§2.5 vanadium null.** Not a separate analysis: same 251-sample Bray–Curtis matrix, same adonis, same 999 permutations, one z-scored V column as the design. R²=0.0031, p=0.780. Riding the identical pipeline is what makes it exclude "any axis returns an R²". **SURVIVES**, and it is the strongest defensive move in the section. Reinforced at family resolution: V, Mo, Ni, Co, Mn, Zn = 0 families at q<0.05 (Cu = 1).

**§2.5 space partialling.** Retention of q<0.05 hits: aridity 85.6%, Cl 82.6%, S 82.0%, Ca 71.1%, **salinity_proxy 41.4%**. **SURVIVES for S/Cl/Ca/aridity; the aggregate proxy does not** — and the paper says so, tying the claim to the mineral components, as §2.4 did for arsenic. Correctly handled.

**§2.5 — two controls NOT run.** (i) *Compartment*: the geochemistry block is conditioned on lat/lon only, while compartment alone explains 4.4% of composition and the plant-root binary axis alone is associated with 60,313 families; the XRF cohort is imbalanced (PlantRoot 127, Deep 74, Surface 50). Mitigating: median salinity_proxy is similar across compartments (0.282 / 0.351 / 0.268) — INFERENCE, my own computation on the env matrix, in no result file. JUDGMENT: residual confound probably modest but untested; conditioning on compartment is a one-line rerun. (ii) *Site*: 251 XRF samples from 59 sites, 456 climate samples from the same 59, permutations free over samples and never restricted within site. JUDGMENT: p-values anticonservative by an unknown factor, R² unaffected; every block reports p=0.001, which is just the permutation floor.

**§2.3 within-phylum replication** (`uv_4_phylum.py`). Eight most abundant phyla; a sample scores for a phylum only if ≥15 of its species are detected; a phylum is tested only if ≥15 paired sites survive (7 qualify); occurrence recomputed using only that phylum's species; Wilcoxon site-paired.
- photolyase **SURVIVES in 3 phyla**, same direction: Actinomycetota +3.47 pp (p=0.0014), Pseudomonadota +2.57 (0.020), Bacillota +2.17 (0.021). Bacteroidota p=0.30 — a fourth phylum under r226, lost when the R232 `Bacteroidota_A` merge added 107 reps (CORRECTIONS 2026-08-17: "Only ONE reported claim moves").
- NHEJ **SURVIVES in 4 phyla**, all surface-depleted: Actinomycetota 0.0025, Bacteroidota 0.0012, Thermoplasmatota 0.0034, Bacillota 0.0127.
- carotenoid replicates in only 2 (Actinomycetota 0.0015, Bacillota 0.00015; Pseudomonadota 0.094) — the paper claims within-phylum replication for photolyase and NHEJ only, and is right not to claim it here.
- **mismatch repair FAILS.** It passed the primary test cleanly (occurrence 50.6 vs 47.0, Holm 0.011; mapped 43.6 vs 39.5, Holm 0.028) and is printed in Supp S4 — but within phyla it moves in *opposite* directions: Pseudomonadota −0.65 pp (0.00058, surface lower) against Bacteroidota +0.89 (0.017) and Chloroflexota +1.28 (0.035). Not claimed. **This is the control doing its job and the best evidence it is not decorative.**
- ⚠ No multiple-testing correction inside that table: 7 phyla × 4 traits = 28 raw p-values. The Pseudomonadota (0.020) and Bacillota (0.021) photolyase results would not survive Holm across it; Actinomycetota (0.0014) would. JUDGMENT: the real argument is direction replicating in 3 of 3 testable phyla — say it that way, not with the p-values.

**§2.3 within-species ecotype test** (`uv_6_ecotype.py`). 327 species have members from both compartments, 324 with resolvable annotations; per species, fraction of surface members carrying the trait vs fraction of deep members, Wilcoxon across species, plus an HQ-only variant (77 species, ≥90% complete / ≤5% contaminated). **All six tests null**: photolyase 0.300 (HQ 0.465), NHEJ 0.866 (0.180), carotenoid 0.631 (0.889). Member completeness balanced (median Surface 86.8 vs Deep 87.4; HQ 95.4 vs 96.3), so the null is not a completeness artefact. Conclusion — species sorting, not intraspecific differentiation — **SURVIVES**. ⚠ But only 37 of 324 species show any non-zero surface−deep photolyase difference: a low-power null, not a demonstration of absence. And "all $p\ge0.18$" is wrong at the third digit — the minimum is 0.1797.

**§2.3 flat-reference control.** The near-universal repair systems are internal negatives: occurrence moves BER −0.05, NER +0.14, SOS +0.74, oxidative +0.43 pp. NER and SOS are nonetheless Holm-significant (0.0047, 0.0010) at those trivial sizes; the paper's "shift by about one percentage point or less" is accurate and wisely silent on significance.

**§2.3 compositional control on Halobacteriota** (`osmo_7_robustness.py`). CLR (pseudocount = half the minimum non-zero fraction) and detection (presence/absence, immune to closure), same partial design. S: CLR +0.376 (1.0e-9), detection +0.441 (3.3e-13); Cl +0.308/+0.369; Ca +0.238/+0.244. **Na fails both** (CLR −0.032 p=0.62; detection +0.094 p=0.14) and is not claimed. **SURVIVES.**

**§2.3 coverage control.** Unmapped fraction Surface 54.3% vs Deep 53.7%, Mann–Whitney p=0.50; PlantRoot maps better (47.3%, Kruskal p=3.9e-10), so root comparisons use the mapped denominator.

## Numbers

| quantity | value | source |
|---|---|---|
| PERMANOVA compartment (547) | R²=0.0443, F=12.61, p=0.001 | `darkfam_ecology_clean/results/permanova_blocks.tsv` |
| space lat/lon (547) | R²=0.0788, F=23.27, p=0.001 | same |
| geochemistry (251) / \| space | R²=0.1231, F=4.25 / R²_marg=0.0802 (partial 0.0881), F=2.90; both p=0.001 | same |
| climate (456) / \| space | R²=0.0582 / 0.0296, p=0.001 | same |
| **vanadium null (251)** | **R²=0.0031, F=0.77, p=0.780** | same |
| core-stratified (≥25% prev, 76,496 fams) | compartment 0.0462 · geochem\|space 0.0828 · climate\|space 0.0305 | same |
| Mantel distance-decay | r=+0.341, p=0.001 | `community_summary.txt` |
| PCoA variance | PC1 15.9% · PC2 10.7% · PC3 7.4% | `pcoa_explained.tsv` |
| Borcard varpart (2 blocks only) | compartment alone 0.0443 / unique 0.0447; space 0.0788 / 0.0792 | `varpart.tsv` |
| per-axis families q<0.05 (% negative) | S 54,227 (91.7) · Ca 47,174 (92.5) · Cl 46,972 (90.4) · salinity proxy 37,280 (86.9) · **Fe 24,110 (95.0)** · aridity 23,031 (95.9) · K 11,935 (58.2) · Na 4,603 (77.4) · **plant-root 60,313 (67.1)** | `indicators_summary.tsv` / Supp S6 |
| trace-metal nulls | V, Mo, Ni, Co, Mn, Zn = 0; Cu = 1 | same |
| space-partial retention | aridity 85.6 · Cl 82.6 · S 82.0 · Ca 71.1 · salinity proxy 41.4 (%) | `partial_summary.tsv` |
| repair floor (% of 5,438) | BER 99.6 · NER 99.5 · oxidative 93.5 · SOS 90.1 · NHEJ 74.5 · **MMR 57.7** · carotenoid 54.4 · photolyase 17.9 | `uv_survival_clean/results/group_prevalence_5438.tsv` |
| compartment, mapped community (44 sites) | photolyase 24.3 vs 12.1 (Holm 0.0019) · carotenoid 70.3 vs 61.0 (4e-5) · NHEJ 80.0 vs 85.7 (0.028) | `compartment_tests.tsv` / Supp S4 |
| compartment, whole-community bound | photolyase 11.1 vs 4.9 (1e-4) · carotenoid 31.1 vs 28.8 (0.0066) | same |
| osmolyte capacity (% of 5,438) | any 79.5 · trehalose(any) 75.4 · hydroxyectoine 18.9 · ectoine 16.5 · glucosylglycerol 3.6 · mannosylglycerate 2.8 | `osmolyte_clean/results/pathway_prevalence_5438.tsv` |
| **Halobacteriota vs sulfur** | **ρ_partial = +0.4515, p = 7.379e-14, n=248**; Cl +0.3850 (3.5e-10); salinity proxy +0.3114 (5.6e-7); Ca +0.2605 (3.3e-5); Na +0.1037 (n.s.) | `osmolyte_clean/results/salt_turnover.tsv` |
| next-strongest phyla on salt | Chloroflexota vs S −0.2416 (1.2e-4); Actinomycetota vs Cl −0.2274 (3.1e-4) | same |
| "osmolyte capacity declines weakly with salinity" | trehalose_any vs salinity_proxy marginal −0.340 (4.1e-8) but **partial −0.104, p=0.10, q=0.21**; breadth vs Na partial −0.190, q=0.016 | `gradient_tests.tsv` |
| **aridity null (group-1l vs precipitation)** | **ρ = −0.172, permutation p = 0.196, n = 58 sites** | `main.tex:283` — inside the `\iffalse` Part A fence, **not printed** |
| osmolyte breadth vs aridity | ρ_partial = +0.023, p = 0.62, q = 0.81 (n=453) | `gradient_tests.tsv` |

Computed and never cited: **Fe** (24,110 families, 95.0% negative — larger than the aridity signal) and the **plant-root contrast** (60,313). Both sit in Supp S6 without comment.

## The aridity null

There is no single aridity null — the phrase covers two results that point opposite ways.
- **Composition is not null on aridity**: climate|space R²=0.030 p=0.001, and 23,031 of 100,293 families track the index at q<0.05, 95.9% negatively, with **85.6% retention under space partialling — the highest of any axis**.
- **Trait-level aridity tests are null**: group-1l hydrogenase share vs precipitation ρ=−0.172, permutation p=0.196, n=58 **sites** (site-level, the honest unit); osmolyte breadth vs aridity ρ_partial=+0.023, p=0.62. §2.2's printed text says only "the within- and cross-biome aridity tests were negative"; the numbers are fenced and printed nowhere. `PAPER_STRUCTURE_PROPOSED.md` R8 already says: "Aridity null is range restriction inside the hyperarid core: say so, print rho -0.172 p 0.196."
- **Is it range restriction?** Partly. Over the 456 climate-annotated samples (59 sites) the de Martonne index runs **0.167 to 1.842** (median 0.702, IQR 0.440–0.828); annual precipitation **6.6 to 74.7 mm**; annual mean temperature **26.3 to 31.35 °C** (computed from `env_matrix.tsv`). Every value sits in the most extreme de Martonne class; the index never crosses a classification boundary anywhere on the 1,043-km transect. JUDGMENT: "range restriction inside a hyperarid core" is defensible **as a statement about aridity class**, not about variation — an 11-fold precipitation span is a real gradient and the dark families respond to it. The honest formulation: the gradient is wide enough to filter community composition and too narrow to move a single guild's genomic capacity; and with 58–59 independent sites a site-level ρ of −0.17 is underpowered (|ρ|≈0.26 needed for p<0.05 at n=58). Do not write "there is no aridity gradient".

## Caveats (quoted)

- M.10: "All censuses are genomic capacities from KO annotation, and medium-quality genomes under-call gene presence, **which biases compartment contrasts toward the null**."
- §2.3: "These are genomic capacities and DNA-based abundances, not measured activity"; "the KO-based panel does not resolve spore photoproduct lyase, 6-4 photolyases or mycosporine-like sunscreens, which remain for profile-based follow-up."
- M.5 (rep-level detection): "That default defines detection: a representative is detected in a sample when its relative abundance exceeds zero, which means at least a tenth of its genome is covered by that sample's reads. **It is an inherited default rather than a chosen threshold**, and the absence language of Sections~2.1 to~2.4 should be read against it."
- M.5 (member-aware detection — the rule §2.5 rests on): "A family is present in a sample when any member gene reaches a covered fraction of 0.5 ... **That fixed threshold is easier to satisfy for a short gene: replacing it with any-read presence reverses the relation between representative length and site occupancy.** Occupancy is therefore compared only within length strata."
- M.1 todo: "reports no make or model for the laboratory instrument, and **states that the concentrations carry no documented unit assertion and must not be read as percentages**. The aridity index used in M.10 is not described in either companion paper and its source still needs naming."
- M.13: "A corrected $p$-value or $q$-value from one block is therefore not comparable with one from another, and no family-wise statement is made over the union of blocks."
- §2.5 close: "this paper makes no functional claim about any dark family." Supp S4: "The mismatch-repair signal fails the within-phylum control (Methods) and is not interpreted as a compartment effect."

## Status

- **§2.3 DONE + SHIPPED 2026-08-13** (Overleaf b4bdac0, GitHub ff4b32c): Results + Methods + Supp S3/S4, registered in `CATALOGUE_NUMBERS.md` §2.3 and `CORRECTIONS.md`.
- **§2.5 DONE + SHIPPED 2026-08-16** (GitHub 32e24ae, Overleaf 4cdfcab): Results + the M.7 environment extension + Supp S6 + `fig:environment` rebuilt on the 106,042 basis; registered in `CATALOGUE_NUMBERS.md` §2.5. Compute closed on both (job 50575481; §2.3 local).
- Open in the record: M.1's XRF/aridity todo (BLOCKED on Alejandra); `fig:environment` and `fig:stress` are both `\status{candidate}`, main-vs-supplementary undecided; `PAPER_STRUCTURE_PROPOSED.md` R8 proposes merging §2.3 and §2.5 into one section — **not executed**.

## What can and cannot be claimed

**Can:** the dark repertoire is environmentally structured (geochemistry 12.3% of compositional variance, 8.0% beyond space, n=251); chemistry structures it beyond distance (distance-decay real at Mantel r=0.34, geochemistry survives conditioning on it); the structure acts by exclusion (87–96% of salt/aridity associations negative); six trace metals associate with no family and the single-axis vanadium PERMANOVA is null (p=0.78) — a same-pipeline negative holding at both resolutions; light partitions communities vertically (photolyase 24.3% vs 12.1% of the mapped community, 44 paired sites) and does so by species sorting; salt is answered by membership change (Halobacteriota vs sulfur ρ=+0.45 after geographic control, robust to CLR and to detection-only coding); the damage-repair floor is not selective (≥99.5% BER/NER, 74.5% NHEJ).

**Cannot:** any functional statement about what a structured dark family does (explicitly deferred to the annotator companion); mismatch repair as a compartment effect (fails the within-phylum control despite being Holm-significant in Supp S4); carotenoid as within-phylum-replicated (2 phyla, Pseudomonadota p=0.094); any aridity claim at guild or pathway level (ρ=−0.172 p=0.196; breadth p=0.62); the aggregate salinity proxy as a chemical axis (41.4% retention — use S, Cl, Ca); "Borcard partitioning separates the environmental blocks" (as run it separates compartment from space and nothing else); sodium as a Halobacteriota axis; any measured-activity statement in either section.

## Reviewer attack surface

1. **Pseudoreplication.** 251 samples from 59 sites, 456 from 59, three compartments and up to three replicates per cell, permutations free over samples. Every PERMANOVA p is 0.001 — the floor. §2.3 does the right thing (collapse to site, pair within site); §2.5 does not, and a referee will notice the asymmetry inside one paper. Fix: restricted permutation within site.
2. **The Borcard claim does not match the code.** Methods promise a partitioning of the environmental blocks; `varpart.tsv` has two rows and neither is chemistry or climate. Rerun properly or delete the sentence.
3. **The geochemistry block omits S, Cl and Na** and includes Si, Ti and V. The 12.3% is not the salt result it reads as.
4. **Compartment is not conditioned out of the geochemistry block**, and the XRF cohort is 51% root-adjacent against 20% surface.
5. **Detection-threshold direction artefact — the strongest attack, currently unanswered.** 87–96% of associations being negative is exactly what a covered-fraction-0.5 presence rule produces if biomass or library complexity falls toward the salty/dry end: families drop below detection, TPM goes to zero, Spearman reads negative. No depth or biomass covariate is in the model.
6. **Units.** The paper cannot state what the XRF numbers are, and the salinity proxy sums three of them.
7. **No pH**, in a paper whose organising variable is salt.
8. **Not rebuildable.** `env_matrix.tsv` cannot be regenerated on IBEX; any availability statement must ship the matrix itself.
9. **Uncorrected within-phylum p-values** (28 tests; two of the three photolyase replications at p≈0.02). 10. **Underpowered ecotype null** (37 of 324 species informative) presented as evidence of no intraspecific structure. 11. Minor: "all $p\ge0.18$" (actual 0.1797); Supp S6 lists the plant-root contrast at n=551 while 4 of those samples have no environmental annotation and are silently coded as non-root in `eco_1_indicators.py`.

## Open questions for Alejandra

1. **XRF units and instrument** — M.1's todo is still BLOCKED on you: what units, which laboratory instrument, and is the trip-5 handheld excluded from `xrf_lab_table_all_trips.tsv`? Nothing else here can be fixed without it.
2. **Aridity index** — name it as the de Martonne index P/(T+10) computed in `build_env_matrix.py`, and what is the source of `AnnualMeanTemp`/`AnnualTotalPrecip` (which reanalysis product, which years)?
3. **Does pH exist anywhere** in the companion ecology data? If it does it belongs here; if not, say so once.
4. **Borcard** — rerun it properly (compartment | space | chemistry on the XRF subset) or drop the claim from M.7 and M.13? Recommendation: drop it; it adds nothing the PERMANOVA does not give.
5. **Rerun the geochemistry PERMANOVA with S, Cl and Na included and compartment conditioned out?** ~4 minutes of compute; it converts the weakest number in §2.5 into the strongest.
6. **The detection-threshold objection** — is a per-sample depth or biomass covariate available to add, or do we concede it in the text? Currently unaddressed.
7. **§2.3 + §2.5 merge (R8)** — execute, and print ρ=−0.172 p=0.196 as the aridity null? Note that number belongs to §2.2's hydrogenase test, and the dark families are *not* null on aridity; the merged section must keep those straight.
8. **Fe (24,110 families, 95.0% negative) and the plant-root contrast (60,313)** are computed, tabulated and uncited — mention Fe in one clause, or state in the caption why the table carries axes the text does not discuss?
