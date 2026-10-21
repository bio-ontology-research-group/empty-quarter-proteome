# Corrections tracker — what changed and why

Division of records: `CATALOGUE_NUMBERS.md` = what is current · **this file** = what changed (numbers,
methods, basis decisions), old → new, with evidence · `reproducibility/claims.tsv` = the PI's
machine-checked mirror — never edited by us, rebindings batched on the PI list in `PLAN.md`, may lag
this file.

**Status:** FROZEN — current result, locked; cite as is · HISTORICAL — superseded, audit-only; do not
recompute, do not headline · RECOMPUTE — needs the current basis · VERIFY — check the basis first ·
RUNNING — job in flight · DONE — corrected in paper + logged.

## The correction loop (per number)

1. **Verify or recompute** on the current basis (5,438 catalogue; 106,042 family universe).
2. **Sweep** for the old value (`grep -rIn` over `paper/ manuscript/ notes/`); fix in `paper/`, then
   regenerate the mirror: `bash reproducibility/sync_manuscript_snapshot.sh`.
3. **Record** — three writes, same sitting: a dated entry appended to the log below · the section's
   table in `CATALOGUE_NUMBERS.md` (old value → its superseded map) · the block row here.
4. `reproducibility/` (claims.tsv, verify.py, ledgers) is never edited in this loop — stale bindings
   go on the PI coordination list.

## Open blocks

Moved to `OPEN_ITEMS.md` on 2026-09-20 (the one open-items register). This file records what changed
and why; it no longer lists what is open.

## Closed blocks — resolved (evidence in the log)

| Block | Outcome | Log entry |
|---|---|---|
| A0 scale / taxonomy | FROZEN — 5,438 / GTDB R232; cite as is | — |
| 573 remnants | DONE — tab:scale row → 551; historical floats fenced | 2026-08-23 |
| B1a structural dark | DONE — recomputed end-to-end on the current dark set; historical 608/1,875 superseded and fenced; §1.3 written | 2026-08-18 |
| C3 global specificity | DONE — recomputed + written on the 106,042 universe with a matched annotated control; SPIRE added as a fourth catalogue; June 21.4% superseded by 28.74% and the abstract's pending-SPIRE hedge closed | 2026-08-17 |
| A4 BGC census | DONE — recomputed + written on 5,438; novelty never computed (deferred) | 2026-08-15 |
| A5 AMR / metal | DONE — recomputed + written on 5,438; STRESS layer recovered; pathogen census restated | 2026-08-15 |
| B1 predicted-3Di screen | HISTORICAL — withdrawn; never revive "89% novel folds" | — |
| Sequence dark, assembly-wide | 68.6M/15.3% → **69.9M/15.6%**; confident (≥100 aa) 5.4% | 2026-08-06 |
| Sequence dark, genome-resolved | 51.9M/2.84M (12.1%) → 88.2M/43.4M/**4.29M (9.89%)** | 2026-08-06 |
| A1 life on air | recomputed on 5,438 then rep_ko v2; §2.2 shipped | 2026-08-11 · 2026-08-12 |
| rep_ko annotation basis | → `rep_ko_v2.tsv`: 0 unannotated, `_sub` contig-subset, 8,134 spurious KOs removed | 2026-08-12 |
| A2 UV survival | recomputed on 5,438/v2, site-paired design; §2.3 shipped (`verify.py` rebind = PI list) | 2026-08-13 |
| A3 osmolyte | first proper census; phantom OR removed; §2.3 shipped | 2026-08-13 |
| B2 dark dN/dS | June 23,382 → rebuilt **106,042**; 0.425 vs matched control 0.388; §1.4 shipped | 2026-08-14 |
| B3 occupancy / abundance | rebuilt member-aware: 98/551, 79.4% all-5-trips; §1.4 shipped (§2.5 use → Open) | 2026-08-14 · 2026-08-15 |
| B4 positive selection | SUPERSEDED — tail = rebuilt universe's own 3.3% vs 0.4%; June FUBAR count retired | 2026-08-14 |
| B5 dark enzymes | CUT → dark-matter-annotator sibling paper; bridge sentence only | — |
| §2.5 gene–environment | DONE — rerun on 106,042 (`darkfam_ecology_clean/`); §2.5 written | 2026-08-16 |
| C1 reference bias | DONE — recomputed on 5,438 / 448 samples (`refbias_clean/`); ceiling 56.7% → **55.6%**; §2.1 written | 2026-08-16 |
| Citation sweep | DONE — 80/80 bibkeys real, 21 added, 3 defects fixed; full record + decision list = `CITATION_AUDIT.md` | 2026-08-18 |
| §1.4 write + Methods M.9 + Supp S5 | DONE — all numbers verified on disk at write time; tool versions in M.9 | 2026-08-14 |

---

## Log — one dated entry per correction, oldest first; append new entries at the end
_(CHANGELOG.md merged in here 2026-08-10.)_

### 2026-08-06 — Whole-catalogue sequence-dark completed
An earlier assembly-wide UniRef50 search left 37,828,539 representatives (8.4%) unsearched. Completing the
search: dark 68,646,640 (15.3%) → **69,938,628 (15.6%)** of 449,265,306 representatives; uncovered → 0.
(source: `darkmatter_clean/04_catalog_dark/dark_summary.tsv`)

### 2026-08-06 — Dark set profiled by length and completeness
`darkmatter_clean/06_dark_lengths/SUMMARY.txt`. The complete 69,938,628-rep dark set contains
**93,335,016 proteins** (corrects the incomplete-search figure 89,414,829). Dominated by short/incomplete
sequence: 52,967,341 reps (75.7%) <100 aa and 60,612,661 (86.7%) incomplete ORFs; only 16,971,287 (24.3%)
≥100 aa and 9,325,967 (13.3%) complete.

### 2026-08-06 — Dark length denominator + §2 reframing
`06_dark_lengths/all_reps_length_summary.txt`. Background over all 449,265,306 reps:
315,606,983 (70%) ≥100 aa; 393,114,691 (87.5%) incomplete. So **incompleteness does NOT distinguish dark reps**
(dark 86.7% vs all 87.5%) — only **length** does (dark 76% <100 aa vs all 30%). Length-filtered dark fraction =
16,971,287 / 315,606,983 = **5.4%** (vs raw 15.6%). §2 corrected: removed the misleading "mostly incomplete
fragments" claim.

### 2026-08-06 — MAG-mapped dark composition
`06_dark_lengths/mag_dark_length_summary.txt`. Of 4,292,446 genome-resolved dark reps:
1,466,943 (34%) ≥100 aa and 1,871,630 (44%) complete ORFs — cleaner than assembly-wide (24% / 13%).
NB protein-count caveat: darkclean_7 `mag_dark_proteins` = 12,674,790 (all cluster members incl. non-MAG);
darkclean_5 `mag_mapped_dark_proteins` = 6,151,650 (MAG-resident only). Use the rep count (4,292,446) as the
headline.

### 2026-08-11 — A1 life on air (§2.2) — **DONE (provisional)**
Recomputed on the 5,438 catalogue (workbench `M/06_atmospheric/lifeonair_clean/`, scripts
`lifeonair_1..4`). Old values were the 15,469 single-sample basis (plain `__SemiBin` naming).

| number | old (15,469, HISTORICAL) | new (5,438, current) |
|---|---|---|
| high-affinity [NiFe]-hydrogenase (HydDB 1h/1l) | 3,386 MAGs (22%) | **946 reps (17%)** |
| form-I CO dehydrogenase | 3,828 MAGs (25%) | **1,147 reps (21%)** |
| guild, either marker | 6,010 MAGs | **1,822 reps (33%)** |
| guild core, both markers | 1,204 MAGs (~8%) | **271 reps (5%)** |
| HydDB group 1l vs 1h (proteins) | 5,292 vs 986 | **1,097 vs 431** |
| community abundance | 13.5% median (16S-weighted, 1,227 amplicon samples) | **18.0% median whole-community; 42.9% of mapped** (MAG-direct CoverM, 548 canonical soil samples) |

Method change registered: abundance is now **MAG-direct** (CoverM over the catalogue) instead of 16S
prevalence-weighting; the 16S cross-check is BLOCKED (amplicon parquet cache not on IBEX). Both denominators
are reported (whole-community incl. ~51% unmapped, and mapped-only). Capacity/prevalence, not flux.
Table 2 (`tab:atmo`) and `fig_lifeonair_by_phylum` remain on the dated 15,469 snapshot by design.
**Provisional**: to be recomputed once `rep_ko_v2` lands.

### 2026-08-12 — rep_ko v2 rebuild + §2.2 recompute — **DONE**
`04_functional/repko_clean/results/rep_ko_v2.tsv` built from 3 fixes: (1) parser handles `48PRr1__SemiBin_N` naming,
(2) 124 previously-unannotated raw bins annotated with the identical stack (Bakta 6.0, eggNOG-mapper 2.1.12), (3)
all 1,010 `_sub` reps switched from parent-inheritance to contig-subset KOs. Unannotated reps 168 → 0; identical
rows 4,580 / 5,438; 690 reps shed at least one spurious KO (8,134 total). §2.2 recomputed on v2:
| number | v1 (was in paper) | v2 (new) |
|---|---|---|
| high-affinity hydB (Group 1h/1l) | 946 (17%) | **968 (18%)** |
| form-I coxL | 1,147 (21%) | **1,170 (22%)** |
| guild (either) | 1,822 (33%) | **1,860 (34%)** |
| both H₂+CO | 271 (5%) | **278 (5%)** |
| HydDB proteins group 1l vs 1h | 1,097 vs 431 | **1,124 vs 445** |
| abundance median whole-community | 18.0% (IQR 12.6–24.7) | **18.4% (IQR 12.8–25.6)** |
| abundance median mapped | 42.9% | **43.9%** |
Post-typing filter dropped 5 `_sub` reps whose marker sat on contigs refinement removed
(hydB: `S54PRr1_..._2_sub`, `V5Dr1_..._31_sub`; coxL: `30PRr1_..._3_sub`, `38Sr1_..._32_sub`, `V10PRr3_..._24_sub`).
Supplementary Note 1 unchanged (guild fraction still rounds to 0.18 at 2 sig figs); Methods M.7 gained one sentence
documenting the correction. PI's `04_functional/refs/rep_ko.tsv` and `reproducibility/claims.tsv` untouched.

### 2026-08-13 — A2 + A3 (§2.3) recompute — **DONE**
Recomputed on the 5,438 catalogue / `rep_ko_v2` (workbench `M/06_metabolism/uv_survival_clean/` and
`M/06_metabolism/osmolyte_clean/`; local runs, no SLURM).

**A2 census, old (15,469 MAGs) → new (5,438 species):**
| group | old % | new % |
|---|---|---|
| BER | 99.5 | **99.6** |
| NER | 99.4 | **99.5** |
| oxidative protection | 91.9 | **93.5** |
| SOS/translesion | 94.2 | **90.1** |
| NHEJ | 78.4 | **74.5** |
| mismatch repair | 48.8 | **57.7** |
| carotenoid pigment | 58.5 | **54.4** |
| photolyase | 17.5 | **17.9** |

**METHOD CHANGE — compartment test.** Old: per-MAG counts, 15,469 genomes treated as independent (MWU
p≈1e-14; photolyase 22.4% vs 13.1% "of MAGs") — pseudoreplicated and recovery-skewed. New: sample-level
community metrics (occurrence / whole / mapped), replicates collapsed per site×compartment, **Wilcoxon
signed-rank paired across 44 sites**, Holm; within-phylum replication control; full grid = Supp Table S4.
New values (mapped community): photolyase **24.3 vs 12.1**, carotenoid **70.3 vs 61.0**, NHEJ **80.0 vs
85.7** (deep-enriched); whole-community conservative bound: photolyase 11.1 vs 4.9, carotenoid 31.1 vs
28.8. **Mismatch repair fails the within-phylum control (opposite signs in different phyla) → NOT
claimed.** Within-species ecotype test (327 spanning species) null → partitioning is species sorting.
Read ladder (546 samples): never-assembled ~44% / assembled-unbinned ~40% / binned ~12%; Surface≈Deep at
every step.

**A3 — first proper census** (old = 4 KOs in `82_metabolism_marker_census.sh` side table, never tested):
13-KO synthesis panel (Suppl. Table S3). any-pathway **79.5%**, trehalose-any **75.4%** (otsA old
67.8 → 66.3 otsAB), ectoine **16.5%** (old 17.0), glucosylglycerol **3.6%**, mannosylglycerate **2.8%**
(both largely Pseudomonadota-only); flat by compartment (trehalose occurrence signal fails the phylum
control). **main.tex:242's old sentence was WIDER than any census — GG/MG had never been censused;
rewritten to match the data.** Salt is turnover, not gene content: Halobacteriota vs S rho +0.451
p=7.4e-14 / Cl +0.385 / salinity proxy +0.311 (n=248 XRF, lat/lon-partialled); robust to CLR and
detection-only variants; detected in 36/248 samples; haloarchaea lack ectoine (5/269) and GG/MG (0/269),
retain trehalose (172/269). Aggregate osmolyte capacity declines weakly with salinity.

**PHANTOM NUMBER REMOVED — "osmolytes not enriched as a class, OR=1.10, p=0.45".** Its cited source
`results/darkmatter_ecology/OSMOLYTE.md` contains no odds ratio and reports MWU **p=0.11** (point-biserial
r median +0.001 vs background −0.015, n=212); no code computes an OR; it is a §2.5 decode result, not an
A3 census result. Removed from §2.3 and `notes/PAPER_STRUCTURE.md`. Pending PI coordination: `verify.py`
`a3.osmolyte.or` carries the phantom value, and `a2.ner/ber/nhej/photo.*` still bind to the superseded
15,469 tables in `results/uv_survival/` — not edited here per the claims.tsv coordination rule.

Written: §2.3 Results, Methods M.8 (stubs M.8–M.11 renumbered to M.9–M.12), Supp Tables S3 + S4.
**Panel correction (zero-effect):** K01859, labelled "phrA photolyase/cryptochrome" in the inherited
`uv_markers.tsv` panel, is **chalcone isomerase** (EC 5.5.1.6, KOfam) — removed from the census panel
(now 28 KOs); its 1 carrier in 5,438 also encodes K01669, so photolyase stays 975 (17.9%), all counts
unchanged. Basis captions added to tab:atmo, fig:bgc, fig:amr ("historical 15,469-MAG single-sample
snapshot").

### 2026-08-14 — Dark-family universe rebuilt (§1.4/§2.5/§3.1) — **DONE**
Basis verification of the June 23,382-family tables found the genome basis correct (current 5,438;
three-way ID match) but the darkness call not the paper's: the decode-era rule (−s 5.7, "no hit OR
<30% identity"), under which 4,701 families (20.1%) are dark only via the identity clause; the
screen had covered only decode candidates — 68% of family-eligible clusters never tested.
Decision: rebuild, not patch. Families derive from §1.2's own dark call (darkclean, 4,292,446
MAG-mapped dark reps), clustered 50%/80%; a family = ≥3 diverged members spanning ≥3 canonical
samples and ≥2 sites.
**New universe: 106,042 families** (536,437 member CDS; paralog-only 2,784 dropped; young
1–2-variant class 782,616). All four legs measured:
- recurrence: median 12 samples / 9 sites; 97.9% single-phylum (Actinomycetota 32%)
- selection: dN/dS 0.425, 96.7%<1 (n=97,200) vs length+size-matched control 0.388, 99.6%<1;
  MWU p≈0. June's larger gap (0.391 vs 0.286) was partly the length confound — the matched
  control is the defensible comparison
- abundance: member-aware recruitment; median family in 98/551 samples; 4 single-sample;
  median peak TPM 99
- persistence: 79.4% present in all 5 expeditions (94.8% ≥4); old 70.3%
- robustness: 130,843 families at 30% vs 108,826 at 50% — threshold-stable
June object superseded for this paper; PI to be briefed, his folders untouched. Review: 8
must-fixes applied pre-run (incl. canonical-sample fix, 106,064→106,042). Workbench:
M/06_darkmatter/darkfam_basis_clean/ (README = audit + job log).
Footnotes: catalogue holds 549/551 canonical samples (51Sr1, F5Dr1 absent); "≥3 dark members"
now true by construction. verify.py rebinding (b2/b3/c3) = PI coordination.

### 2026-08-15 — §1.4 abundance denominator completed, 550 → 551 — **DONE**
F13PRr1's member-aware CoverM profile landed and the family × sample matrix was refreshed to 551
columns. **Every §1.4 digit unchanged** (98 / 4 / 99 TPM / 79.4% / 94.8%); only the denominator
moves to "98 of 551" and the incomplete-profile caveat is retired. B4 retired with the §1.4 write
(fast tail = the universe's own 3.3%-vs-0.4%; June 148-family FUBAR count not used).

### 2026-08-15 — §2.4 recomputed + written on 5,438 — **DONE**
Both censuses aggregated on the completed arrays (antiSMASH and AMRFinder; 5,438/5,438
each, logs clean) and §2.4 written (three paragraphs: biosynthetic census / intrinsic resistome +
recovered metal layer / biosafety control), Methods M.10, both figures rebuilt on the current basis
(`code/make_bgc_amr_figures.py`). Old → new (old = n_MAGs/15,469 redundant bins, 58% of current reps
never scanned): terpene 64.5% → **62.14%**; NRPS 25.3% → **25.06%**; ectoine 17.0% → **16.53%**
(kept: terpene-precursor separate, 83.78%); β-lactam 24.67% → **24.55%** — remarkably stable across
the basis change. New measured numbers: 97.48% of species ≥1 BGC; 44,002 regions; **contig-edge
85.29% measured** (replaces the inherited "90%"); AMR any-class 37.40%; STRESS layer recovered
(As 17.34%, Hg 7.52% — the historical promotion dropped Type==STRESS, METAL.md). Honest correction:
"arsenic tracks the salinity proxy (r=+0.26)" attenuates under geographic control (partial ρ
0.24→0.12, q=0.11); the FDR-robust associations are Ca (ρ=+0.224, q=0.011) and aridity (ρ=−0.243 vs
precipitation, q=2.4e-6) — §2.4 states the Ca/Na + aridity version. Pathogen census RESTATES on the
species basis: old "7 MAGs / 0.05%" → **2 named pathogen species of 5,438 (0.04%)** (dereplication
collapses 5 E. coli + 2 A. baumannii bins to 1+1). VFDB reused by reference from the PI's
current-basis census (06_resistome_safety/vfdb/, verified by MAG naming): 5,435 → 202 MAGs at
id≥40 → ≥90. Verification armor: ectoine antiSMASH×KO Jaccard 0.921 (860/895/899); trehalose
KO×gapsmith 0.795; gapsmith binary ectoine call insensitive (key-reaction gated), continuous
completeness concordant (median 60 KO+ vs 40 KO−); gapsmith all-Pathways.tbl verified pre-gap-fill
(log step [1/5] vs [4/5]) → safe for capacity claims. Workbench: bgc_clean/, resistome_clean/,
pathway_check_clean/ (new). Deferred, explicitly: BGC novelty (deferred list), VFDB rerun not needed,
nitrogenase/rhodopsin recompute → §2.2's story, PI list (vnfD absent on old basis).

### 2026-08-16 — §2.5 gene–environment rerun on the 106,042 universe — **DONE**
The June B3 environment analysis (23,382-orthogroup object, repo `results/darkmatter_ecology/`,
2026-06-28) is replaced by a rerun on the rebuilt dark-family universe: workbench
`M/06_darkmatter/darkfam_ecology_clean/`, scripts = verbatim copies of the June
`phase_e_indicators.py` / `phase_e_partial_osmolyte.py`(a) / `phase_b_ecology.py` + `eco_stats.py`
repointed to `darkfam_basis_clean/results/fam_abundance.tsv` (106,042 × 551) and the shipped §2.3
`env_matrix.tsv`. Old → new (all June values on the superseded object): PERMANOVA geochem R²
0.066 → **0.123** (beyond-space 0.051 → **0.080**); climate|space 0.033 → 0.030; compartment
0.020 → 0.044; distance-decay Mantel 0.161 → **0.341** (all p=0.001); vanadium block null in both
(p=0.78). Family-level q<0.05 counts: Ca 4,066 → **47,174** (92.5% negative) · S 3,964 → **54,227**
(91.7% neg) · Cl 3,594 → **46,972** · salinity proxy 3,179 → **37,280** · aridity 491 → **23,031**
(95.9% neg), of 94,568 testable (XRF n=251) / 100,293 (aridity n=456); ~10× June because the
universe is 4.5× larger and median prevalence 98/551 clears the ≥25-present gate. Space-partial
retention: aridity 85.6 / Cl 82.6 / S 82.0 / Ca 71.1 / salinity proxy 41.4% (aggregate proxy most
spatially confounded — matches §2.4's Ca/Na+aridity restatement). Direction-stability: axis
ranking, ~90%-negative filter shape and V/Ni/Co/Mn/Zn nulls all reproduce. CHANGE: June's Mo
signal (177 families, 100% positive) does NOT reproduce (0) — never a paper claim. Verified:
top-Ca-family Spearman independently reproduced with scipy.stats.spearmanr (−0.525 exact);
eco_1/eco_2 marginal counts cross-agree. Scope (forced): EC-decode layers + osmolyte vignette not
ported (no decode on this basis; sibling annotator paper — `a3.osmolyte.or` stays a PI verify.py
item); plant-host genus axes not rerun (amplicon cache off-IBEX); plant-root contrast retained
(60,313/102,005 q<0.05). §2.5 written this session; fig:environment rebuilt on the new basis
(`code/make_darkfam_env_figures.py`); Supp Table S6 = per-axis association grid. June dirs
untouched (banner-marking = PI item 2b).

### 2026-08-16 — C1 reference bias (§2.1) recomputed — **DONE**
The §2.1 "BLOCKED on the PI's PICRUSt rebuild" label was **inherited and stale** — every input already
existed: CoverM profiles (2026-08-06, `picrust2_current_20260803/coverm_current/profiles_full_recompute/`),
the corrected MAG–KO table (2026-08-11, `rep_ko_v2.tsv`), and the staged, checksummed 16S side
(`picrust2_current_20260803/inputs/`; PICRUSt2 v2.4.1, 10,543 KOs × 1,270 profiles + per-profile NSTI,
`INPUTS_COMPLETE`; `ko_metagenome_unstrat.tsv` sha256 `c71c2256…`, byte-identical to the companion amplicon
repo's `input_manifest.tsv` pin). What is stalled is the PI's own **fail-closed 573/573 gate**: 4 of 573
CoverM tasks failed (`51Sr1`, `F5Dr1`, the `SPositive` control, `F13PRr1` on a time limit), his downstream
jobs sit in `DependencyNeverSatisfied`, and his `mag_ko/` + `comparison/` were never produced — not a
dependency for us; his chain was not touched. **The amplicon parquet cache is NOT a §2.1 blocker** (it
blocks only §2.2's 16S abundance cross-check); the chain reads real `.biom`/`.fasta` inputs under
`/ibex/project/c2014/EmptyQuarter_Data/`.

Recomputed in `M/06_reference_bias/refbias_clean/`; all 5,438 reps resolved, 0 failures.

| number | old (15,469 bins / 436 samples, HISTORICAL) | new (5,438 reps / 448 samples) |
|---|---|---|
| no-KO ceiling | 56.7% (CheckM2 denominator) | **55.6%** of 19,186,853 protein-coding genes |
| ceiling, same CheckM2 denominator | 56.7% | **58.4%** of 20,445,320 |
| ceiling, HQ only | — | **53.9%** (2,259 reps) |
| comparison cohort | 436 | **448** |
| weighted NSTI | 0.16 | **0.159** |
| presence-recall of MAG KOs | 93% | **86.0%** |
| MAG KO richness / sample | 4,482 | **8,666** (PICRUSt2 constant 10,543) |
| Spearman ρ, shared KOs | — | **0.839** |
| encoded KOs outside the reference | — | **2,517 of 11,374 (22%)**; 14.0% per sample |
| over-prediction ratios | 2.40×/1.11×; 31.41×/14.63× | **withdrawn** (asymmetric denominator, unchanged) |

Denominators are not interchangeable: the numerator comes from Bakta gene calls scored by eggNOG, so 55.6% is the
proteome-consistent value we quote, while CheckM2's own calls run ~6% higher — the ceiling survives on either
convention, and is not a genome-quality artefact. `_sub` reps use `repko_2_submap.py`'s contig-subset
construction, not the parent superset. Cross-check: our numerator for `10Dr2__single_SemiBin_0` (1,754)
reproduces the historical `dark_fraction.tsv` row exactly — only the denominator differs, by the caller gap.
The 93% → 86.0% recall move is a different cohort *and* basis; no continuity claimed.

**Sample-cohort note — the counts in this log are not interchangeable.** §2.1–§2.3 use the **rep-level** CoverM
run (`profiles_full_recompute`, mapping to the 5,438 representatives), where `51Sr1`, `F5Dr1` and `F13PRr1` all
failed → **548** usable of 551 canonical. §1.4's abundance uses the separate **member-aware** run against family
members, which reached **551** once F13PRr1 landed (entry 2026-08-15). Both are current; they are different runs.

Written: §2.1 + new **Methods M.11** (old M.11/M.12 renumbered M.12/M.13). C1 scaffold retained, re-tagged
HISTORICAL; SUPERSEDED (its table/figure keep the old numbers, now explicitly labelled). Companion-paper
PICRUSt2-vs-shotgun reconciliation → `DISCUSSION_NOTES.md`, not Results.

### 2026-08-17 — Read-mapping presence rule is length-biased
**No number changes.** `fam_presence.tsv` calls a dark family present when any member gene reaches CoverM
covered fraction ≥0.5 (M.9). A fixed breadth threshold is easier to satisfy on a short gene, so occupancy
from this matrix partly measures length: median representative length falls from 145 aa at 1–3 sites to
**37 aa** at 56–61 sites. Demonstrated, not inferred — swapping the rule for `TPM>0` reverses the
relation (76 → 175 aa). **Any distribution claim from this matrix must be length-stratified.**
- **§1.4 is the only affected claim and it survives.** Stratified <60 aa → ≥250 aa: median **201 → 55**
  of 551 samples; all-five-expedition persistence **89.7 → 67.1%**; ≥4 expeditions ≥90.8% in every
  stratum; all four single-sample families are <60 aa. `98/551`, `4` and `79.4%/94.8%` **stand as
  computed** — 98 is composition-weighted and spans 201–55, so quote it with that floor. The
  anti-relic-DNA inference is strengthened, not weakened.
- **§1.4's other legs, §2.5 and §2.1–2.4 are unaffected**: universe, recurrence, dN/dS and phylum are
  assembly- or annotation-based; `fam_6c_dnds.sh` length gates cut dark and its matched control
  symmetrically (97,200 vs 95,647); §2.5 uses `fam_abundance.tsv` gated on `TPM>0`; §2.1–2.4 use
  genome-level mapping — tested: median samples detected by genome size 7 / 35 / 43 / 44 / 42 from
  <1.5 Mb to ≥5 Mb, flat above 2.5 Mb and *lower* for small genomes, so no size analogue.
- **Corrected here:** §2.5's gate was described three ways; Results "detectable" and Supplementary
  "present" → **"non-zero abundance"** (M.9's "non-zero" was already right).
- **Owed:** M.5 is a `\todo` stub governing §2.1–2.4; its 10% covered-fraction floor is CoverM's
  never-passed default. Parameters and evidence: `M/06_darkmatter/specificity_clean/README.md`.
- **Dropped (never in the manuscript):** the planned §3.1 occupancy × dN/dS analysis — dN/dS coverage
  falls 96.0% → 59.7% across occupancy bins (short families dropped, bins differentially filtered);
  reinforces the `RQ/B6/result.md` walk-back. Replacement result: reference-absent families are **not
  local rarities** (median 36 of 61 sites vs 38, holds in all five length strata).

### 2026-08-17 — Pipeline-guide provenance corrections
**No number changes.** `notes/PIPELINE_GUIDE.md` rewritten from the thirteen `*_clean` dirs, verified
on disk — it owns the tool parameters. Paper-facing facts:
- **GUNC never ran** (no artifacts); **DAS_Tool DID** (v1.1.7; 542 → 539 → 14,426 refined bins).
- **R232 does NOT come from `_gtdbtk_refined`** (that is GTDB-Tk 2.4.1 / r226); the 2.7.2/R232 census
  comes from the PI rebuild's `legacy_taxonomic_screen_r232/`. `CATALOGUE_NUMBERS.md` §1.1 had cited
  the r226 path — corrected; numbers unaffected.
- **All four CoverM thresholds are inherited defaults**, never passed (`--min-covered-fraction` 10%,
  `--contig-end-exclusion` 75, identity/aligned-percent 0) — M.5 must not present them as choices.
- Cohort wording: **"548 usable of 569 profiled"**, never "of 573" (4 abandoned `.tmp` dirs); all 569
  are fresh recomputes, zero legacy reuse.
- `rebuild_20260729` is not blanket-dead — the "ignore" rule covers its *coassembly numbers* only.

### 2026-08-17 — §3.1 specificity recomputed — **DONE**
Basis `M/06_darkmatter/specificity_clean/`; registry `CATALOGUE_NUMBERS.md` §3.1. The June object (23,382
orthogroups) is superseded in full. New in this recompute: the 106,042 length- and size-matched annotated
control families are searched in the *same* DIAMOND pass, so every catalogue reports a known-protein
absence floor; SPIRE is added as a fourth catalogue (its June failure was an `ARG_MAX` concatenation bug,
not missing data).

| number | old (23,382 basis) | new (106,042, dark / matched control) |
|---|---|---|
| **no homolog in any named catalogue** | 21.4% (5,011) | **28.74% (30,479) / 0.41% (435)** |
| MGnify90, 717,738,164 | 36.6 | 52.68 / 1.36 |
| GEM, 122,978,191 | 47.4 | 64.38 / 2.97 |
| five deserts combined, 147,531,178 | 31.9 | 42.60 / 2.07 |
| SPIRE representative genomes, 269,496,580 | never ran | 39.90 / 0.98 |

**The abstract's "pending recompute / pending SPIRE" hedge is CLOSED.** Two controls now carry the claim.
(i) An E-value-corrected rarefaction gives the absence a size-matched random MGnify90 slice would produce;
at full size it reproduces the observation to 0.15 points, and against it the desert union leaves 24.33
points less of the dark proteome unexplained than its size predicts, SPIRE 21.30 points less, with control
gaps near 1 point. No subsampled databases were built. (ii) Against the complete 1,661,142,385-sequence
MGnify90 2026_07 set (MMseqs2, a different aligner, so reported separately and never unioned with the
DIAMOND sets), 27.08% remain unmatched; that reference resolves 34.3% of the control's residual absences
but 5.8% of the dark set's.

Measured limits, both stated in M.12: searching the five desert catalogues separately rather than
concatenated returns hits for 2,544 further families (2.40%) from the same sequences, so raw
cross-catalogue ordering is descriptive only; and `atacama`'s gap is not interpretable as content, its
*control* gap being −29.97 because a single 35-sample study loses to a random draw even for conserved
proteins.

Written: §3.1 (reframed from exceptionalism to reference coverage per `PAPER_STRUCTURE` 2026-08-17) + new
**Methods M.12** (old M.12/M.13 renumbered **M.13/M.14**); `tab:specificity` rebuilt with a control column
and SPIRE; new **Supp Table S7** = the size-matched expectations; M.9 now discloses its length-dependent
presence rule; the abstract's stale "recomputation pending" sentence deleted (§2.1 shipped 2026-08-16).
Display (Alejandra delegated the call): one table + one figure + one supp table, no duplicated data;
**`fig:endemism` CUT** (its dangling `\ref` in the frozen B3 scaffold repointed to `fig:specificity`; no
scaffold claim altered); generator `code/make_specificity_figures.py`; the superseded generators are left
in place for audit and **must not be rerun**; `reproducibility/regenerate_figures.sh` deliberately NOT
modified (the PI's — coordination). Owed bibkeys `bay2021aridity` + Aplakidou 2025 → resolved 2026-08-18.

### 2026-08-17 — Methods written; M.7–M.11 renumbered
**No result number moves.**

**Renumber — frozen, last one. Entries above use the OLD numbering:**
M.7→**M.9** (§2.2) · M.8→**M.10** (§2.3) · M.9→**M.7** (§1.4) · M.10→**M.11** (§2.4) ·
M.11→**M.8** (§2.1) · M.12–M.14 unchanged. Block is M.1–M.5 general, M.6–M.12 in section order
(M.6 = §1.3), M.13–M.14 back matter.

- **GUNC deleted from M.2** — it never ran. Contamination control is CheckM2 + DAS_Tool scoring.
- **§1.1 site count 60 STANDS** — a change to 61 was made and reverted. 61 counts field-notebook
  labels (59 transect sites + 2 trip1-only named locations); both companion 16S papers say 60.
  Do not re-derive 61 from the metadata.
- **542 → 539 is a funnel**, not a discrepancy: 542 submitted to DAS_Tool, 539 produced bins.
- **M.5 defines "detected" once**: relative abundance > 0 = **≥10% genome breadth, an inherited
  CoverM default, not a chosen threshold**. Cohort is **548 usable of 569 profiles**; 573 is a
  directory count (4 abandoned `.tmp` dirs), not a cohort.
- **168 unannotated reps decompose 124 + 44** (the earlier draft had the causes wrong): 124 never
  annotated → re-annotated 124/124; 44 `48PRr1__SemiBin_N` only appeared unannotated because the
  naming convention broke the join. Do not merge the groups.
- eggNOG-mapper executable **2.1.13**, reference data **2.1.12** — both quoted, not harmonised.
- Five tool-parameter facts the guide had wrong (foldseek cluster not MMseqs2; two deliberate Foldseek
  builds; pyrodigal, ProstT5 and ESMFold parameters) → corrected in `PIPELINE_GUIDE.md` §6.1, §8, §13.

### 2026-08-17 — GTDB r226 → R232 everywhere (Alejandra: r226 was an error)
Fixed at the root: R232 output subset to our 5,438 reps at `M/05_refine/_gtdbtk_r232/`, 17 taxonomy
paths in 8 consuming scripts repointed from `_gtdbtk_refined` (= r226); all reran clean.
**Only ONE reported claim moves:** §2.3's within-phylum photolyase control replicates in **three**
phyla, not four — **Bacteroidota falls to p=0.30** (the `Bacteroidota_A` merge adds 107 reps). All
other shipped numbers verified unchanged; the pathogen census is R232 and true. §3.1 wording:
"61 sites" → "61 sampled locations"; the site count stays **60**.

### 2026-08-18 — B1a structural novelty (§1.3) recomputed — **DONE**
The structural chain finished on the current basis (workbench `M/04_structure/00_structure_clean/`).
It replaces the PI's Aug-2 result as the paper's structural
number: his was the June 68,527-model cohort on the old May tier1 basis re-scored against current references,
and only 0.87% of our queries appear in it.

| number | old (B1a, HISTORICAL) | new (current dark, 2026-08-18) |
|---|---|---|
| model cohort | 68,527 historical candidate-family reps | **17,221** admitted fold families (13,770 folded + 3,451 reused) |
| gate survivors (Q) | 1,875 | **2,838** |
| qTM dual no-hit (K) | 608/1,875 (32.43%) | **1,031/2,838 (36.3%)**, 95% CI 34.6–38.1 |
| alignment-TM dual no-hit | 97/1,875 (5.17%) | **95/2,838 (3.35%)** |
| control floor, qTM / alnTM | 101/3,013 (3.35%) / 68/3,013 (2.26%) | **unchanged — reproduced exactly** |

The control floor reproducing to the digit is a genuine cross-check, not a copy: the same 3,684 blinded models
went through our own invocation with the same references and parameters, so the instrument is demonstrably the
one his calibration describes.

**What the recompute changes scientifically.** The claim is now explicitly normalisation-dependent and that is
the finding, not a caveat: whole-protein novelty is common and far above its floor (36.3% vs 3.35%, 10.8×,
p < 10⁻²⁰⁰; ≈936 models after subtracting the floor), while *local* novelty is rare and only marginally above
its own floor (3.35% vs 2.26%, 1.5×, p = 0.011). Consistently, 96.7% of confident models contain an aligned
region matching a known structure but only 63.7% are covered end-to-end, and candidates are longer than other
qualified models (301 vs 215 aa) at the same confidence (pLDDT 0.763 vs 0.782). §1.3 therefore reports
structural novelty as **architectural rather than local**.

Also corrected this session: the **abstract** carried the historical 608/1,875 as the paper's structural result
(now the current numbers); the **B1a scaffold paragraph was live, not fenced** — it is now fenced with the
other ported blocks; the **spine figure** panel 3 carried the historical cohort (now current). The historical
Table (B1/B1a) is explicitly captioned historical and is left as the audit record.

### 2026-08-19 — Methods trim (M.1–M.6)
**No result number moves.** Three factual errors fixed:
1. **M.4 denominator**: 9.89% is of the **43,423,322** reps in an HQ/MQ bin, not "of the catalogue"
   (against 449,265,306 it would be 0.96%).
2. **M.3 dRep**: contamination ×5 **is** the default, not a deviation — the run script passes no
   scoring weights.
3. **M.4 pyrodigal**: "translation table 11" was never passed; it is the meta-mode default read off
   the GFF header.
**Sample-count policy**: the paper states only **551** (cohort) and **549** (assemblies); 573 / 552 /
569 dropped from prose (552's 3 extra libraries contribute 4,425,298 of 753,305,533 proteins —
workbench record).

### 2026-08-19 — Abstract rewritten on current numbers
Dead R226 block (573 / 51.9M / 23,416,070 / 12.1%) replaced by the current numbers throughout
(registry: `CATALOGUE_NUMBERS.md`); dN/dS 0.425-vs-0.388 flagged as not stronger constraint.

### 2026-08-20 — Introduction rewritten
Fixed a live self-contradiction (pre-reframe "specific to the Rub' al-Khali" vs the reframed §3.1);
now on framing (b), PI skeleton fenced not deleted. Two factual corrections (Alejandra):
**"contiguous" RESTORED** in body ("largest contiguous sand sea"; title left for the end-game
true-up) and **551 = shotgun SEQUENCED, not collected** (>1,000 collected; exact count owed to M.1).

### 2026-08-23 — Leftover superseded numbers stopped printing
**No result number moves.** tab:scale's "573 / 551 / 550" row → "551 (60 sites, five expeditions)"
and its unregistered "coassembly 527/126" row removed; fencing the historical floats stopped the
superseded 608/1,875 from printing.

### 2026-09-17 — §3.1 occupancy recomputed on the corrected union (36 → 37)
The `b6_*` outputs had been written before the SPIRE leg landed, so they split the universe
35,158 / 70,884 (the pre-SPIRE union) while `dark_nohit_union_all.ids` already held the corrected
30,479. Rerun of `b6_occupancy_join.py`, same inputs and definitions (551 columns, 61 sites,
106,042 families). **Occupancy of the reference-absent set: median 36 → 37 of 61 sites**, against 38
for families with a hit; length strata 44/50 · 38/40 · 33/**36** · 31/32 · **26/27** (was 33/37 and
25/28). Class sizes 30,479 / 75,563; median dN/dS 0.4718 / 0.4232. The gap narrows and the
conclusion is unchanged: absent families are marginally less widespread, confined in none, so the
caption wording at main.tex L1331 still holds. Registry row updated.
**Still stale, not yet fixed:** main.tex L234 prose ("a median of 36 of the 61") and
`figures/fig_specificity_controls.pdf` panel b, whose generator
`code/figures/make_specificity_figures.py` reads the corrected TSV and needs a rebuild.

### 2026-09-17 — cultured-vs-MAG dark comparison paired by species (PoC, not a paper number yet)
Both sides were pseudo-replicated: 132 cultured genomes resolve to 23 distinct MAG representatives,
one of them 31 times, so each side was weighted by cultivation effort. `step3_04_summarise.py` now
pairs by species. **Cultured lower in 16 of 23 pairs, median -0.33 pp (Wilcoxon p=0.042, sign test
p=0.093);** distinct-rep MAG side 1.73 % pooled against cultured 0.99 %. Superseded: 1.17 % pooled,
1.03 % median, "cultured lower in 107 of 132". **The residual gap is cultivation depth:** pairs with
≥10 cultured strains (n=4) give median +0.07 pp, single-strain pairs (n=10) give -0.49 pp, so the
difference is strain sampling, not assembly. Pooling vs averaging strains changes nothing (-0.33 vs
-0.39 pp, p=0.042 vs 0.035). Darkness is not an artefact; qualify "the same rate" with the
single-strain caveat.

### 2026-09-17 — §1.3 known-fold floor standardised on confidence (10.8x → 6.6x; alnTM 1.5x withdrawn)
Controls are higher-quality models than candidates (pLDDT 0.867 vs 0.775) and their false-call rate
falls from 5.49 % to 1.53 % across the pLDDT range, so the pooled floor was measured where errors are
rarest. `05e_stratified_floor.py` reweights the control rate onto the candidate distribution.
**qTM: 10.84x → 6.56x (4.97–9.08) on pLDDT, 8.09x (5.96–11.49) on pLDDT × length** — report
standardised. **alnTM: 1.48x → 1.38x, CI 0.87–2.41, crosses 1** — withdraw as a positive result.
Crude values reproduce exactly. Architectural-not-local claim unaffected. Not yet registered.

### 2026-09-17 — §2.5 geochemistry PERMANOVA block corrected (S, Cl, Na added)
The block omitted sulfur, chlorine and sodium, which carry the family-level signal, while including
Si, Ti and V, which are associated with no family and where V is the declared null block; and the
Borcard partition admitted climate only when complete for every sample (456 of 547) and never
offered geochemistry, so it partitioned compartment and space alone although M.7 and M.13 describe
a partition of the environmental blocks. Both fixed in `eco_3_community.py` and rerun.
**Geochemistry R² 0.123 → 0.132 marginal, 0.080 → 0.089 beyond space** (n=251, p=0.001); the old
eight-axis block is retained as `GEO_ELEMS_V1` and reproduces the registered values to four
decimals. Compartment is not a confound (0.0893 → 0.0916 conditioning on space plus compartment).
Partition on the XRF cohort: geochemistry unique R² 0.106, space 0.055, compartment 0.057.
Coverage caveat: Mn measured in 54 of 251 samples, V in 163. **Not yet registered** pending
Alejandra's choice of which block the paper reports.

### 2026-09-20 — lineage darkness, archaeal matched null and cultivation control REGISTERED
The exploratory analyses `eqmags_poc/01_dark_by_lineage` (round two), `03_genome_size_function` and
`05_rakpr_overlap` were promoted to `M/06_darkmatter/lineage_darkness_clean/` and
`M/06_darkmatter/cultured_comparison_clean/`, path constants repointed at the root, and rerun; every
result file is byte-identical to the exploratory copy. Registered as new sections R2 (cultivation
control), R3 and R4 in `CATALOGUE_NUMBERS.md`. Superseded on the way in: the 5,359-rep rename basis
(round one), the pseudo-replicated cultured comparison (1.17% / 1.03%, 107 of 132 pairs → paired by
species, n = 23, 16 of 23, −0.33 pp), and the novelty-matched Thermoplasmatota contrast (0.310 vs
0.125) as the headline contrast → 1:1 matched null (6.5x vs archaea, 4.5x vs bacteria). The
UniRef50 leg of the cultured comparison reuses the hit tables of job 51696940 until rerun 52131486
lands. Same session: §2.5 corrected block (0.132 / 0.089), the standardised structural floor
(6.6x / 8.1x), 58.3%, 419 of 448 and the 37-of-61 occupancy are now in the registry and in the v2 draft.

### 2026-09-20 — dN/dS gains a calibrated null (frame-shifted arm) REGISTERED
`darkfam_basis_clean/scripts/fam_6c_dnds.sh` gained ARM=null_frame1 (same dark member CDS read in the
+1 frame, stop codons masked to NNN, all other filters identical) and ARM=dark_subset (the same seeded
10,000 families in the true frame); dark and ctrl code paths verified byte-identical by regression.
Null median **1.087** vs true frame **0.427** (full arm 0.425); paired shift +0.655, null higher in 92.4%.
The R2 sentence "the distance from 1.0 is uncalibrated" is withdrawn from the draft. Caveat carried: a
+1 frame of a coding sequence inherits its composition and is a null for "is the annotated frame the
selected one", not for neutral evolution in general.

### 2026-09-20 — PERMANOVA: within-site permutations added, no R² changes; every block keeps p ≤ 0.002; the Methods sentence conceding free permutations is replaced.

### 2026-09-20 — ⚠ `rep_relative_abundance.tsv` and `v2_reps.tsv` differ in row order: join by genome name, never by position.

### 2026-09-21 — cultivation control rerun on all 722 cultured genomes (was 651: 71 soil-batch proteomes live in
`EmptyQuarter_Data/cultures/annotation/thermophilic/bakta/`, not in `eqmgc`). 23 pairs unchanged; cultured lower
16 → **18 of 23**, Wilcoxon 0.042 → **0.021**; shared-species rate 0.99 → **0.86%**; re-found 37.8 → **40.9%**;
control 46.5 vs 0.16 → **46.4 vs 0.18%**; cultured dark median 60 → **59 aa**. R2 and Methods updated.

### 2026-09-22 — R4's absences-only physiology paragraph replaced by the encoded-metabolism paragraph
(`v2_thermo_metabolism.py`, completeness-matched comparators); the mcr/amo/rhodopsin absences stay inside it.

### 2026-09-22 — culture–MAG pair anatomy REGISTERED (step3_06): within a species, annotated proteins share a
reciprocal partner 80% of the time and dark proteins 16%, in 23 of 23 pairs both ways; MAG-unique remainder =
short-read fragments, culture-unique = accessory content. Two sentences added to R2.

### 2026-09-20 — shotgun blanks profiled; two reagent contaminants CONFIRMED (contamination_clean)
The negative control SNTC puts 8.56% of its reads on exactly two representatives, the *E. coli*
genome (4.30%) and a *Cutibacterium acnes* genome (4.26%, 290/550 soil samples, DNA rho −0.47, never on
the pathogen panel). The positive control (*Hafnia paralvei*, 120 M reads, not in the catalogue) puts
zero reads on any representative, ruling out index hopping. The exploratory contamination check
reproduces byte-identical. Consequences: both representatives carry a contamination annotation in the
draft (removal from the 5,438 is the PI's call); the biosafety paragraph is not written; the SNTC
assembly is a laboratory plasmid with three resistance markers and matches contigs in 10
representatives, so the AMR census must be checked for vector-derived hits before it is defended
(open item).

### 2026-09-20 — Pfam on the structural candidates: "known parts in new arrangements" WITHDRAWN
Pfam-A 38.2 (gathering thresholds) on the 4,880 members of the 1,031 candidate fold families
(`M/04_structure/pfam_candidates_clean/`, job 52136381): 91.25% carry no Pfam domain; 86.0% of families
have no Pfam-bearing member. The locally familiar region that Foldseek finds in 96.7% of confident
models is structural resemblance, not a recognisable sequence domain. Where a domain exists (8.75%) it
covers a quarter of the chain and leaves a median 159 aa unmatched, and K exceeds the length-matched
other confident models in every 150–599 aa bin, so the multi-domain reading holds only for that
minority. R5's closing sentence and the Discussion's "recombination of known domains" reading are
rewritten to: whole-protein architectures with no covering match, whose parts resemble known folds
structurally but are not recognisable sequence domains. Caveat: K is sequence-dark by construction, so
the K-vs-control Pfam gap is partly built in; the K-vs-Q-not-K contrast carries the signal.

### 2026-09-24 — genome-mapping cohort 548 → 550 (census and stress analyses rerun)

The three canonical samples without a CoverM representative-level profile (`51Sr1`, `F13PRr1`, `F5Dr1`) now
have one, but `F5Dr1` (1,232 read pairs after filtering) maps no read and its profile is all-NaN, so the
usable cohort is 550 of 551, not 551. `lifeonair_clean`, `uv_survival_clean` and `osmolyte_clean` were rerun
on 550 into `results/v551/` (scripts gained a `RES_OUT` override and a guard that skips zero-mapped samples;
registered outputs untouched). Promoted: paired sites 44 → 45 (site 51 gains its surface profile); photolyase
24.3/12.1 → 25.0/12.3, carotenoid 70.3/61.0 → 70.7/60.9, NHEJ 80.0/85.7 → 79.7/85.6 (mapped share, %);
unmapped 54.3/53.7 (p 0.50) → 54.4/53.7 (p 0.47); Halobacteriota vs S ρ 0.451 → 0.456 (p 3.0e-14), Cl 0.385 →
0.392, salt cohort n 248 → 250. Unchanged: guild share 18.4 %, species counts 968/1,170/1,860. Not rerun: the
climate-axis gradients (n = 453) and the arsenic gradients in `resistome_clean` (n = 248); Figure 6c and 6e still
show the 548 basis until rebuilt. Superseded: 548, 44 sites, 24.3/12.1, 70.3/61.0, 80.0/85.7, 0.451, 0.385, n 248.

### 2026-09-24 — known relatives of the 237 Thermoplasmatota REGISTERED (lineage_darkness_clean, v2_thermo_relatives)

Comparators: the 218 GTDB r232 species representatives of c__SW-10-69-26 (21) and c__E2 (197) plus the 83
MGnify genomes in those classes (54 duplicate a GTDB rep and are excluded from nearest calls). SMAG's 49 genomes
could not be fetched per genome and are not included. Result: 0 of 237 at ANI ≥ 95 %; best AAI median 64.8 %
(same-genus median among GTDB reps 75.4 %); 24 reach 70 %; median AAI to the nearest other desert rep 73.1 %.
New R4 paragraph and Methods sentences written. AAI is a new small step (reciprocal-best-hit DIAMOND on
Pyrodigal proteomes) because no CompareM/EzAAI exists on the system. Downloads (5.0 GB) deleted; manifests kept.

### 2026-09-24 — MIMAG rRNA/tRNA on the 5,438 REGISTERED (new folder `M/mimag_rrna_trna_clean/`)

barrnap 0.9 + tRNAscan-SE 2.0.12 in domain mode (Slurm array 52468548). Under the full MIMAG rule 126 of the 5,438
representatives are high quality (16S assembled in 495; all three rRNAs in 157; ≥18 tRNA isotypes in 3,629). The
paper keeps completeness/contamination tiers (2,259 / 3,179) as operational and now discloses the MIMAG count in
R1 and Methods; the former "not assessed" sentence is replaced.

### 2026-09-24 — SMAG added to the reference-absence search (specificity_clean, SMAG leg)

SMAG (Ma 2023; 40,039 genomes as listed in its Supplementary Data 2, from the Zenodo 8223844 archive, which holds
40,352 files including 312 timestamped duplicates) was gene-called here (Pyrodigal -p meta, 129,550,340 proteins)
and searched under the identical DIAMOND 2.0.15 rule. Dark 40.47 % absent vs control 1.57 %; +27.69 points over the
size-only expectation, the largest gain of any multi-study catalogue (desert set +24.33, SPIRE +21.30). Union with
SMAG 26.83 / 0.36 %. The shipped headline union (28.74 / 0.41 %, 1.258 G) is UNCHANGED; SMAG is written as a
separately reported sixth catalogue in R6 and Methods until the depth check (spec_5) is rerun on the with-SMAG set.
Genomes and intermediates deleted; only the DIAMOND db (54 GB, /ibex/user/velazqal/eqmags_31/smag/) is kept.

### 2026-09-24 — depth check rerun on the with-SMAG absent set (specificity_clean, spec_5 with TAG)

Job 52498970. Of the 28,446 dark families absent from the union with SMAG, 1,546 (5.43 %) gain a hit in MGnify90
2026_07 against 138 of 385 control absences (35.8 %); 25.37 % of the 106,042 dark families remain unmatched by any
reference searched (control 0.23 %). The hits are identical to the original run's on the same families. Registered
beside the SMAG row; the paper still reports the old set (5.8 / 34.3 / 27.08) until decision A11.

### 2026-09-27 — SMAG JOINS the headline reference-absence union (decision A11, Alejandra)

Headline absence is now over MGnify90 ∪ GEM ∪ SPIRE ∪ SMAG ∪ five desert catalogues (1.387 G sequences): dark
26.83 % (28,446) / control 0.36 % (385), replacing 28.74 % (30,479) / 0.41 % (435) on the 1.258 G union. Depth
sentence now 5.4 / 35.8 / 25.37 % (was 5.8 / 34.3 / 27.08). Discussion gain range "24 to 28 points" (desert set,
SMAG). Occupancy of the absent set recomputed on the 28,446 (see registry). Figure 5b needs the SMAG point and 5c
the new occupancy values. Superseded: 28.74, 30,479, 0.41, 435, 1.258 G, 27.08, 5.8, 34.3.

### 2026-10-12 — P3 life on air: review round 1 (Astra) corrections
Source: `lifeonair/reviews/round1/` (33 findings) and `response_round1.md`; code `lifeonair/code/fix_coxL_v2.py`,
`h123_v2.R`, `h1_mechanism_v2.R`, `apply_sub_fix.py`, `b1_classifier_benchmark.sh`, `s1_sub_reconcile.sh`;
workbench `M/06_atmospheric/lifeonair_review_r1/` (jobs 53611804 genome size, 53613713 benchmark, 53614016 `_sub`).
Superseded P3 result folders: `results/h123/`, `results/h1_mechanism/` (v1). Current: `results/h123_v2/`,
`results/h1_mechanism_v2/`, `results/a5_genome_traits/`, `results/benchmark/`.

| number | old | new | basis |
|---|---|---|---|
| CoxL form-I rule | AYRGAGR = form I (inverted) | **AYXCSFR = form I; AYRGAGR = form II** | King & Weber 2007; Cordero refs carry AY[RA]CSFR |
| form-I CoxL proteins / reps | 1,901 / 1,170 | **284 / 245** | corrected motif + retained contigs |
| form-I coxMSL co-occurrence (complete context) | 73% | **96%** | operon_summary_v2.tsv |
| bulk median form-I GE | 0.616 | **0.105** | GE v2 |
| group 2a | in other-NiFe control | **own column (bulk 0.03 GE)** | Islam 2020; R1-02 |
| sampling frame | 550 samples, 61 locations | **546 samples, sites 1–60 (59 with data)** | sites 62/64 special-interest |
| root/bulk visits | 107 at 48 sites | **106 at 47** | sites 1–60 |
| root/bulk 1l log2 | −2.23 | **−2.27 (CI −2.90 to −1.59)** | v2 |
| root/bulk 1l excess | −1.97 | **−1.83 (CI −2.62 to −1.00)** | controls without 2a, CoxL other incl. form II |
| surface/bulk form I | −0.02 ("unchanged") | **−0.51; excess −0.24** | v2 |
| H2 claim (CO capacity independent of position) | stated | **WITHDRAWN** | corrected CoxL |
| MAP product | Open-Meteo default best_match (ERA5+IFS), 14–48 mm | **pure ERA5 1995–2024, 13.3–19.2 mm** | best_match re-query not reproducible |
| bulk 1l vs MAP rho | −0.51 (51 sites) | **−0.37 (52 sites); 6-block CI −0.65 to 0.16** | pure ERA5, sites 1–60 |
| temporal window | 71 visits/27 sites, 17 months (trip 5 missing) | **103 visits/35 sites, 31 months** | ERA5 to 2025-10-31 |
| temporal claim "unrelated to antecedent rain" | stated | **WITHDRAWN** (not equivalent; surface/root 1l lower after wet 90 d) | within/between model |
| `_sub` reps | parent-bin annotations | **retained contigs only**: 1l carriers 599→597, 2a 136→134; 891/1,010 growth re-estimated | s1_sub_reconcile |

### 2026-10-12 — P3 review round 2 corrections

Sources: `lifeonair/reviews/round2/response_round2.md`, regenerated `h123_v2/`,
`h1_mechanism_v2/`, `operon/`, and new `benchmark_r2/` results. R2 supersedes the
corresponding R1 values above.

| Quantity | R1 | R2 | Source |
|---|---|---|---|
| Spatial bulk/surface/root site counts | 52/46/60 (duplicated site 19) | 51/45/59 | H3_spatial_by_marker_position.tsv |
| Bulk group 1l ERA5 correlation | -0.37 | -0.401; six-block CI -0.64 to 0.07 | same |
| Surface group 1l POWER/blended correlation | described as robust | -0.292/-0.255, block CIs span zero | same |
| Root decomposition turnover share | unqualified mechanism | 97.4%, bootstrap 96.7–98.3%, catalogue accounting | A2_root_bulk_decomposition_summary.tsv |
| Bulk aridity observed/turnover-only correlation | -0.380/-0.367, missing zero cells | -0.430/-0.433, compositional closure | A2_aridity_decomposition_summary.tsv |
| Within-site 90-day rainfall cross-product agreement | untested | rho 0.193 | H3_temporal_product_agreement.tsv |
| All-position expedition/site/position R2 | 2.0/28.8/4.6% | 2.1/30.4/4.5%, visit-aggregated | stability_permanova_restricted.tsv |
| Complete-context group 1l loci | included removed-contig proteins | 445; small-subunit 72.8%, maturation 87.4% | operon_summary_v2.tsv |
| GE calibration | assumed copies/genome | marker/SCG estimator; 1l mock log2 bias +0.67, paired-ratio mean error -0.11 | benchmark_r2/mock_ge_* |

The group 1l root/bulk median (-2.275 log2 units; 106 visits, 47 sites) and main
excess contrast (-1.827) are unchanged. Literature and release gaps are documented
in the round-2 response; frozen-evidence verification is not a public raw-data rebuild.

## 2026-10-12 — P3 round 3 (Gemini Pro)

- Missing expedition-4 calls → seven candidates with sequences recovered from the September-9 control replay (source checksum in lifeonair/results/review_r3/trip4_provenance.json). Four additional MAG exclusions, union 383 → 387 including two orphans.
- Previously untested MAG exclusions → genus/trait/CBB refits saved in review_r3; group-1l genus coefficient -0.812 → -0.770 with phylum, -1.035 → -0.964 with order. Main interpretation retained.
- Expanded read exclusion: root/bulk group-1l -2.275 → -2.285; excess -1.827 → -1.760 across the same 106 visits/47 sites. Primary estimates unchanged; these are sensitivity medians.
- Claimed sequence-disjoint benchmark/coverage-only explanation → graph-component holdout with residual homology: 34/6478 comp70 queries exceed 70% at half-length coverage, max 83.3%; permissive verification detects omitted graph edges. Caption and Methods corrected; no benchmark rerun claimed.
- Wang 2025 “global soil survey” → northern-China gradient, checked against author-uploaded full text; comparison bounded across lineage/functional definitions and growth methods.
- Explicit climate-grid/finite-bootstrap qualification. Public deposition and publication metadata remain open. Evidence and point-by-point dispositions: lifeonair/reviews/round3/response_round3.md.

## 2026-10-13 — dark-matter three-round review, rounds 1–2 dispositions

User authorised three external Opus/Flash reviews with root revisions and chose Nature Communications for rounds 2–3. Evidence and reports are under `reviews/darkmatter-2026-10-13/`; current checks under `reproducibility/darkmatter/`. Superseded submitted drafts and outputs are retained.

- **DONE: coding-evidence scope.** 106,042 admitted families contain 536,437/4,292,446 representatives (12.5%). Constraint and recurrence findings apply to this selected subset. Genome-resolved short/partial fractions are 65.8%/56.4%, distinct from assembly-wide 75.7%/86.7%.
- **DONE: genome-quality nomenclature.** Completed MIMAG table gives 126 high/5,312 medium, while the CheckM2-only near-complete tier has 2,259. The statement that RNA content was unassessed was false. Source: frozen `mimag_tiers.tsv`; original barrnap/tRNAscan code preserved.
- **DONE: selection comparisons.** Pooled family medians .4253 dark/.3884 control are divergence-confounded. Coarsened matching retains61,427 per arm (.3862/.4229); direction depends on target weights and divergence. Independent1000-family bootstrap and target-weighted contrasts are documented; no universal comparability claim. Jobs53678755 and53679366.
- **DONE: frame null.** Actual median1.08645 (rounds to1.086),37.6% below1,9,148evaluable null/9,182annotated/9,029paired. Audit50,898CDS found21internal stops, not widespread contamination.
- **DONE: culture comparison.** Excluded confirmed E.coli pair, leaving22. Protein shares65,478/81,723 annotated vs238/1,507dark. First nucleotide job53678768 lacked minimap2 `-c`; approximate-identity results are SUPERSEDED. Base-alignment job53679350 gives dark-unshared516/1,105, annotated-unshared4,877/16,056, dark-shared219/236, annotated-shared62,925/65,477. Locus audit53680202 retained. No accessory-content inference from protein no-hits alone.
- **DONE: spurious-call sensitivities.** AntiFamv8 flags141admittedmembers/56families; opposite-strand non-dark overlap flags70members/46families. Union102families; exclusion leaves97,109evaluable darkfamilies,median.4251,96.7%below1. Jobs53678944/53680096. Overlap is conditional on existing Pyrodigal calls; it is not exhaustive alternative-ORF validation.
- **DONE: structural controls.** Job53679366 independently verifies22dark/66controlmodels have no returnedrows. Conditional local-score failures73/2,816 vs2/2,947; oldstandardisedratio1.38 is not used to infer equivalence. Of1,031query-thresholdfailures797have>=100aa,>=50%coverage localmatches. Novel-fold/whole-protein-absenceclaimsdemoted; maximaindependent. Reviewer prefiltermechanism remainsuntested.
- **DONE: taxonomy/cohort.** 61sampledlocations=59core+2auxiliary. Firstemptyrank and firstempty/digit-bearinglabel are separate rules; tableadded. Thermoplasmatota do not own the largest absolute dark-protein count.
- **DONE: environment sources.** Metadata, geodata and laboratoryXRFtables retained; portablebuilder reconstructs547-rowmatrix byte-identically. Open-Meteo visit-centred12-monthsummary provenance identified. Laboratoryunits/LOD/zeros andoriginalclimatemodel/rawresponses remain unresolved. Exploratoryecology movedtoSupplementaryNoteS1; siteaggregate sensitivity53679245doesnotresolveallconfounders.
- **VERIFY before submission:** publicENA/revieweraccess, MAGdeposit/archiveDOI,78missingrawpaths/302permissionfailures, exactduplicate-runmergecommands, remainingwhole-paperregistrymigration. Current35revisionchecks areexplicitlyscoped; legacyverifygatesremainstaleandarenotclaimedtopass.

## 2026-10-13 — final (round 3) revision

- All three Opus5.5reports completed; final dispositions are in `reviews/darkmatter-2026-10-13/round3/response.md`.
- Flagged-family calibration added:91evaluablefamilies, median dN/dS0.7622,70/91below1; AntiFam alone47families,median0.8307,36/47below1. Excluding flags leaves the original pooled median almost unchanged but does not validate coding status.
- Final culture job53680886 extends the locus audit to all DNA-supported calls:344/735dark vs607/67,802annotated have no BaktaCDS over half the aligned span. Protein-or-DNA support is separately752/1,341dark vs70,354/81,533annotated. Only65eligible dark calls map into47admitted recurrentfamilies.
- Retained Thermoplasmatota enrichment4.6-fold: exact source value4.551874 rounds to4.6; the review's4.5 used rounded inputs.
- Removed persistent-detection inference from the Discussion lead, narrowed rarefaction to the named catalogue results, fixed figure labels and cohort descriptions, and replaced publication TODOs with explicit pending-access/deposition statements.
- Current compact verification covers49revisionclaims and manuscript regression guards. It is separate from the historical whole-paper registry. No submission or public release has occurred.

## 2026-10-13: Gemini 4 Argon fourth review

See `reviews/darkmatter-2026-10-13/round4/response.md` for accepted and rejected
findings. PICRUSt reference saturation is now explicit (419/448 full-set outputs).
Thermoplasmatota metabolic contrasts use the same 181 genomes in each matched arm;
full-lineage prevalences remain separate. Depth medians are paired-site statistics.
IBEX job53695688 corrected soil-position filtering:547samples,101897testedfamilies,
20136enriched/41046depleted; propagated phylum and detection summaries replace prior
551-column binary-screen values. Corrected partial-correlation df leaves57/224FDR
hits (previously58), with identical rho values. rRNA complete-hit preference changes
zero calls in5438GFFs. PERMANOVA core-model logged sample counts are251/456. Original
outputs remain frozen; corrected artifacts are in reproducibility/darkmatter/
round4_corrections/. Fulltext count remains77. The compact driver now has89checks;
239input/source hashes pass. Historical whole-paper gates remain unreconciled.

## 2026-10-13: fifth-review exploratory extensions

IBEX jobs 53697241/53697281 add three post hoc tests, with designs recorded before
execution, frozen inputs/results and seeds under `reproducibility/darkmatter/round5_extensions/`.
Caller substitution at fixed DNA-supported culture loci changes dark no-CDS calls
344/735 -> 235/735; recurrent-family subset stays 18/42 (one gain, one loss).
All-five-expedition detection 79.4% compares with a site-by-position label-null
mean 82.7%, an observed-minus-null difference of -3.36 percentage points.
PICRUSt fixed-common-universe mean rho is 0.8431 versus 0.8245 under the
expedition-by-position null; adding site conditioning gives 0.8431 versus 0.8408
on 146 exchangeable samples. These are distinct estimands from median shared-KO
rho 0.839. At KO-fraction threshold 1e-3, 446/448 MAG profiles have zero retained
KOs and undefined recall. No coding validation, biological persistence or
predictor false-positive claim follows from these exploratory diagnostics.

A dangling Supplementary Table S3 reference DOES occur in the osmolyte Methods;
the round-4 response's assertion that it was absent was incorrect. It is replaced
with the exact 13-KO panel and marker/strict-route rules from the retained source.
Seven bibliography entries gain verified volume/issue/pages from their retained
DOI-matched primary XML, with an itemised audit under round5/.

The fifth Argon review prompted a whole-visit temporal-label sensitivity (IBEX
53698522): observed79.38%, nullmean80.54%, central95%range78.95–82.11%, lower-tail
p0.11. The deficit seen when permuting individual samples is sensitive to within-visit
replication. Both results remain reported. Added independently aggregated confidence
imbalance and Pfam dark-representative comparison, corrected the second dangling
Extended Data Table1 reference, and reported the12flagged-representative instance
counts and median-darkness exclusion sensitivity. Abundance-weighted effects remain
unassessed by this compact exclusion.

## 2026-10-14: H1–H4 extensions and structural-score disposition

The round 7 extension source is `manuscript/main.tex`; evidence and independent review
are under `reviews/darkmatter-2026-10-13/round7/`. H1 reports conditional within-lineage
PGLS associations; H2 retains its inconclusive adjusted between-group recurrence
contrast; H3 reports exact-length, equal-count, deduplicated SMAG/GEM searches over
the full 212,084-query census. H3 raw audit 53796328 reproduced all 16 summary and 80
stratified rows from 16 raw hit files and verified 14 matched databases.

Historical alignment-normalised structural threshold counts (95 dark, 68 controls;
797 high-coverage local matches) are retired from manuscript interpretation. The
Foldseek 10 source audit identified missing inclusive-span normalisation, which also
changes length-dependent scoring parameters. Neither clipping nor multiplying by a
length ratio is an exact repair. Original scores remain in the diagnostic archive.

The accepted structural comparison uses qTM throughout and pools the 1,031 initial
qTM-low dark models independently of alnTM. Same-run CPU PDB searches recover crop
matches for 714/964 whole-chain-unmatched models; 775/1,031 have any matched crop.
The GPU prefilter sensitivity gives 731/961 and 794/1,031, respectively. CPU source-entry-
excluded controls have 2,889/3,013 crop matches. Independent GPU raw audit 53796103 passed.
Figure 4 is regenerated by `reproducibility/darkmatter/round7_figures/make_structure.py`.
Cropping changes content, length-dependent normalisation and search opportunity jointly.

## 2026-10-15: family follow-up after independent Argon review

F3 is now a bounded known-family case study, analysed separately from the UniRef50-dark classification. The expanded Bulzu 2021 panel has 1,566 profile/site-qualified sequences, and all 45 EQ site-conserving HeRs have qualifying homologues. Read support is 43/45 loci; the two failures remain flagged. Profile-core addition to the expanded reference persists, with mean 7.33 clusters per ten Saccharimonadia sequences, but GEM overlaps the reference directly. New reference-free clustering gives 8.04 EQ versus 4 GEM core clusters at 50% identity and 9.605 versus 10 at 80%; threshold dependence is reported. These metrics do not establish biochemical function or model-training novelty.

F1 contributes a strict CPF census and exact B+W+R decomposition, not independent replication of the existing photolyase association. The median surface contrast is 8.47 percentage points (Holm across six contrasts, p=.00131). Unclassified genus groups comprise 25.7% of mean abundance and have zero within-group contribution by construction.

F4 supersedes the earlier 623 A/B-neighbour count with 616 loci requiring complete same-strand A/B calls, in 608 MAGs. Of these, 609 have distinct A/B ORFs and 7 have non-overlapping A/B domains in a single ORF. Corrected marker-order labels are in results/review_checks/neighborhood_gene_orders.tsv; tie-sorted B-A labels in the original table are not physical gene order. The ecology contrast remains unresolved (median +0.56 points, Holm across six contrasts, p=.389). Five EctCD fusion predictions have original-read junction support; one is FEYH and outside the WEYH subset.

Authoritative follow-up root: /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean. Sources, parameters and accepted jobs are documented in reproducibility/darkmatter/family_followup_20261015/REPRODUCIBILITY.txt; independent reviews and dispositions are in reviews/darkmatter-family-followup-20261015. New values are FROZEN compact-result extractions in verify_family.py and the central claims registry. The separate Bakta and Prodigal protein denominators are explicitly distinguished in Methods.

## 2026-10-15: biological argument and full author list

Rewrote the title, Abstract, Introduction, Results and Discussion around the biological questions: recurrent unknown families, their concentration in Thermoplasmatota, subsurface distribution, cross-soil sequence sharing, structural relationships and heliorhodopsin diversity. Moved procedural controls, calibration and interpretive assumptions to Methods; consolidated biological limitations and proposed experiments in Discussion. Retained numerical results in the appropriate analysis sections and restored denominators/cross-references after Gemini 4 Argon High proofreading. No new analysis or change to numerical evidence. Review, baseline and root disposition: `reviews/darkmatter-biology-rewrite-20261015/`.

Added the 19-author byline from Life on Air, with Alejandra first, remaining middle authors alphabetical by surname, and Magnus Rueping then Robert Hoehndorf last. The Life on Air manuscript has its own checkout at `~/Documents/papers/empty-quarter-lifeonair/` and Overleaf project `6ad0d1b4338ca7cd7d291188`; its analysis bundle remains in this repository.

### 2026-10-21 — record joined with the 20 October rewrite; cohort residue cleared; PI numbers registered

Record: the September registry rows (SMAG union, cohort 550, relatives, MIMAG), corrections and open items carried
into this repository; the old repository's main merged and archived. Manuscript base = Overleaf main a65f9aa.

Supplement cohort block set to the registered 550-cohort values (S323–329, S504, S821–831): 44 → **45** sites;
photolyase 24.3/12.1 → **25.0/12.3**; carotenoid 70.3/61.0 → **70.7/60.9**; NHEJ 80.0/85.7 → **79.7/85.6**; unmapped
54.3/53.7 p 0.50 → **54.4/53.7 p 0.47**; Halobacteriota ρ 0.451 n 248 → **0.456 n 250**, Cl 0.385 → **0.392**;
salt axes n 248 → **250**; climate axes n 453 → **455**. The 548-cohort paragraph replaced by one cohort sentence.
Not bound to any check in `reproducibility/verify.py`; the harness stays green.

Gradient reruns on the 550 cohort (`resistome_clean/results/v551/`, `osmolyte_clean/results/v551/`; scripts gained
`RES_OUT`): arsenic ~ Ca +0.224 → **+0.226** (q 0.0093, n 250); ~ precipitation −0.243 → **−0.242** (q 2.6e-6, n 455);
climate-axis tests keep 7/7/0 significant rows per axis. Superseded: 248, 453, +0.224, −0.243.

PI October analyses registered with the file that prints each value (new section "PI October analyses"): PGLS,
flanking recurrence, SMAG–GEM 18.54, crops 714/964 (derived), culture DNA support and 22 pairs, AntiFam/shadow 102,
visit-bundle null 80.5 %, Jaccard 0.634, frame null 1.086, plant-root sets 20,136/41,046, contaminant exclusion
6.31 → 6.34 %, chaperone families, heliorhodopsins, CPF, EctC; visit-matched ecology 482/550, 36 of 37 sites,
+0.46 pp (`06_darkmatter/natcomms_20261020`). Superseded in place: 481/549, 40 of 45, 23 pairs, 1.087,
102,005/19,818/40,495. Occupancy wording: 61 location labels, 60 sites.
