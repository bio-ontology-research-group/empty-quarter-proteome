# darkfam_ecology_clean — §2.5 gene–environment rerun on the 106,042 universe

Owner: velazqal. Created 2026-08-16 (session 04-dark-function).

## STATUS (2026-09-20): PERMANOVA hardened against nesting and sequencing depth

Two referee objections to the community block, neither of which moves an R². (1) The 547 samples
sit in 59 sites, up to three positions and three replicates per site visit, so permuting sample
labels freely over samples treats nested samples as independent units: pseudoreplication, and
every p was 0.001. (2) TPM removes library size from each family's value but not from what the
library detected, so an environmental block could be a sequencing-depth gradient. Chain rerun as
job **52137904** (4m04s, eco_0+1+2+3+4, both selftests PASS).

`eco_stats.adonis` and `varpart` take a `strata` argument that permutes sample labels only within
a stratum; `strata=None` keeps the free permutation and the seed, so every registered p is
reproduced bit for bit. `eco_3` runs every block twice, free and within site (the env table's
`site` column, 59 levels, identical to `sample_site_trip.tsv` after its `site` prefix), and adds
`p_within_site` to `permanova_blocks.tsv` and the two within-site p columns to `varpart.tsv`; the
registered columns of both are untouched. `eco_0_depth.py` writes `results/sample_depth.tsv` from
the 573 fastp reports, because no depth table existed anywhere on disk.

**Every environmental block survives the restricted test.** Nesting was not doing the work.

| block | n | R² | p, free | p, within site |
|---|---|---|---|---|
| compartment | 547 | 0.0443 | 0.001 | 0.001 |
| space (lat/lon) | 547 | 0.0788 | 0.001 | **1.000, not testable** |
| climate | 456 | 0.0582 | 0.001 | 0.001 |
| climate \| space | 456 | 0.0296 | 0.001 | 0.001 |
| geochem | 251 | 0.1322 | 0.001 | 0.002 |
| geochem \| space | 251 | 0.0893 | 0.001 | 0.001 |
| geochem \| space+compartment | 251 | 0.0916 | 0.001 | 0.001 |
| geochem v1, no salt | 251 | 0.1231 | 0.001 | 0.001 |
| geochem v1 \| space | 251 | 0.0802 | 0.001 | 0.001 |
| vanadium (declared null) | 251 | 0.0031 | 0.780 | 0.884 |
| compartment [core] | 547 | 0.0462 | 0.001 | 0.001 |
| geochem \| space [core] | 547 | 0.0915 | 0.001 | 0.001 |
| climate \| space [core] | 547 | 0.0305 | 0.001 | 0.001 |
| depth, log10 bases | 547 | 0.0113 | 0.001 | 0.001 |
| geochem \| space+depth | 251 | 0.0889 | 0.001 | 0.001 |
| aridity \| space | 456 | 0.0064 | 0.001 | 0.001 |
| aridity \| space+depth | 456 | 0.0060 | 0.001 | 0.001 |
| climate \| space+depth | 456 | 0.0297 | 0.001 | 0.001 |

Partition, same two schemes: compartment p 0.002 alone and unique, within site 0.002 and 0.002;
geochem 0.002 and 0.002, within site 0.004 and 0.002; space 0.002 and 0.002, within site 1.000
and 0.576.

**Space is the one block a within-site test cannot address, by construction.** lat/lon are site
properties, constant across every sample of a site in all 59 sites, so permuting within site
leaves the design exactly invariant, every permuted F equals the observed F and p is 1.000. That
is the test refusing to answer, not the signal failing: a between-site term needs whole-site
permutation, which the site sizes (1 to 28 samples) make non-trivial, and it is left open. Space
is also the only term that loses its unique variance under the restricted partition (p 0.576).
Every environmental block varies inside a site and is testable there.

**Sequencing depth is not the geochemistry signal.** Depth on its own structures the community
(R² 0.0113, p 0.001), and conditioning on it changes nothing that matters:

| block | R² beyond space | R² beyond space+log10 depth |
|---|---|---|
| geochem, 11 axes, n=251 | 0.0893 | 0.0889 |
| aridity, n=456 | 0.0064 | 0.0060 |
| climate, 3 axes, n=456 | 0.0296 | 0.0297 |

Depth covers all 547 samples (fastp `after_filtering.total_bases`; log10 min 5.26, median 10.33,
max 11.68, i.e. 0.0002 to 476 Gbp, median 21.6 Gbp). Caveat: three libraries are below 2 Gbp
(`F5Dr1` 1,232 reads, `F45PRr3`, `F13Sr1`) and `F5Dr1` and `F13Sr1` are in the XRF cohort, so the
low tail of log10 depth rests on few samples; nothing was dropped.

**Two defects fixed at the root on the way.** (a) `adonis` compared permuted F to observed F
exactly, so an invariant term had its ties split by summation noise of order 1e-10 and returned an
arbitrary p (0.219 in the wiring test) where the honest answer is 1.000; the comparison now allows
a 1e-9 relative tolerance, as vegan's `permutest` does, and no free-permutation p moved. (b)
`eco_3 --selftest` wrote its fabricated tables into `results/`, replacing the registered ones; the
selftest now writes to a temporary directory.

**Registered numbers unchanged, verified by md5.** 47 of the 50 files in `results/` are
byte-identical to the pre-change run, including every `indicators/` and partial output, both PCoA
dumps, `INDICATORS.md` and `PARTIAL.md`. The three that changed are the eco_3 tables that gained
columns or rows, and their registered columns diff clean: geochemistry 0.1322 / 0.0893 / 0.0916,
climate beyond space 0.0296, compartment 0.0443, space 0.0788, vanadium 0.0031 at p 0.780,
partition 0.1058 / 0.0548 / 0.0571, Mantel r +0.341, all 16 indicator counts, all six retentions.

**Open:** the between-site test for the spatial term (whole-site permutation under unequal site
sizes); `code/figures/make_darkfam_env_figures.py` in the paper repo titles fig 2 "all p=0.001
except vanadium", which is now only true of the free permutation, so the caption needs the
restricted p or a footnote (paper repo not touched here).

## STATUS (2026-09-20): full per-family tables persisted, plant-root axis extended

Why: the plant-root contrast was computed but only its top 40 families reached disk, so every
follow-up had to recompute the axis. `eco_1` now persists the full table for each axis, `eco_2`
partials the plant-root contrast on the same footing as the five master axes, and a new step
`eco_4` describes which families the contrast picks out. Chain rerun as job **52135803** (3m13s,
eco_1+2+3+4, selftest PASS).

| new output | content |
|---|---|
| `results/indicators/<axis>_all.tsv` | every tested family, 16 axes: family_id, spearman_r, p, q, n_present (93 MB total) |
| `results/compartment_families_summary.tsv` | counts, effect sizes, phylum shares, over-representation, persistence |
| `results/compartment_families_by_phylum.tsv` | per-phylum cross-tab, 38 phyla |

**Registered numbers unchanged.** Indicators identical (plant-root 60,313 of 102,005 at q<0.05,
67.1% negative; S 54,227 / Ca 47,174 / Cl 46,972 / salinity 37,280 / aridity 23,031; same nulls),
every `*_top.tsv` byte-identical, retention identical (aridity 85.6 / Cl 82.6 / S 82.0 / Ca 71.1 /
salinity 41.4), and the `eco_3` tables checked (community_summary.txt, permanova_blocks.tsv,
varpart.tsv, pcoa_explained.tsv) byte-identical to the 09-17 rerun. Only `partial_summary.tsv` and
`PARTIAL.md` changed, by the added row.

**Plant-root under space partialling:** n=547, marginal 61,182, partial 59,820, **retained 94.1%**,
median abs(dr) 0.011. It is the least spatially confounded axis in the panel. Its n is the 547 samples
carrying lat/lon, not the 551 abundance columns `eco_1` scores, so its marginal count here (61,182)
is the recomputation on that subset, not the registered 60,313.

**Compartment families (`eco_4`, from the persisted axis):**

| set | families | median rho (IQR) | in all 5 expeditions | leading phyla, % of set (% of all attributable) |
|---|---|---|---|---|
| root-enriched | 19,818 | +0.154 (+0.121, +0.200) | 86.8% | Pseudomonadota 19.2 (9.6), Actinomycetota 16.8 (32.8), Myxococcota_A 11.3 (4.9), Thermoplasmatota 11.2 (19.8), Planctomycetota 8.3 (2.4) |
| root-depleted | 40,495 | -0.186 (-0.235, -0.140) | 89.8% | Actinomycetota 52.7, Chloroflexota 13.2 (8.1), Thermoplasmatota 11.3 |
| all families | 106,042 | | 79.4% (baseline) | |

Over-representation in the root-enriched set (share of set / share of all attributable, min 20
families): Bacillota 3.62, Planctomycetota 3.41, Bdellovibrionota 3.27, Cyanobacteriota 2.98,
Verrucomicrobiota 2.94, Myxococcota_A 2.32, Pseudomonadota 2.01. Phylum attribution uses the 77,069
single-phylum families of 78,728 attributable; 13,954 root-enriched and 26,993 root-depleted
families are attributable.

Both sets persist better than the universe, so the compartment signal is not an artefact of
one-off families. Difference from the read-only reproduction that proposed this step: none. Every
count, share, ratio and quartile reproduces; only Planctomycetota's share of all attributable
families reads 2.4% here, from 2.45% at two decimals.

**Open:** not registered yet; the compartment result is §2.3 territory by the shipped division of
labour, so where it lands is Alejandra's call.

## STATUS (2026-09-17) — geochemistry block corrected in `eco_3`, rerun

Two defects (paper repo `notes/dossiers/05_environment_stress.md`): the block omitted S, Cl and Na,
which carry the family-level signal, while including Si, Ti and V, which carry none; and varpart
admitted climate only when complete for every sample (456 of 547) and never offered geochemistry, so
it partitioned compartment and space alone although M.7/M.13 describe the environmental blocks.
Both fixed in `eco_3_community.py` and rerun; `GEO_ELEMS` is the 11-axis block, `GEO_ELEMS_V1` the
old one kept for reproducibility; varpart now runs on the XRF cohort where all three blocks exist.

| block, n=251 | marginal R² | \|space | pseudo-F |
|---|---|---|---|
| v1, no salt (reproduces registered) | 0.1231 | 0.0802 | 4.25 |
| corrected, +S/Cl/Na | 0.1322 | 0.0893 | 3.31 |
| corrected \| space+compartment | — | 0.0916 | 2.53 |

All p=0.001. Reproduction exact to four decimals, so the difference is the design matrix alone.
Headline moves 12.3 → 13.2 %: the omission neither inflated nor deflated the result, but the block
now contains the salt axes and can be described as one. Compartment is not a confound (0.0893 →
0.0916). Partition on the XRF cohort: geochem unique 0.1058, space 0.0548, compartment 0.0571
(p=0.002) — chemistry carries about twice the unique variance of geography.
Coverage caveat: Mn measured in 54/251 samples, V in 163; S 222, Cl 219, Na 241.

**Open:** not registered pending Alejandra's choice of reported block; XRF units/instrument still
blocked on her (M.1).

## STATUS (2026-08-16) — §2.5 SHIPPED (GitHub main `32e24ae`, Overleaf `4cdfcab`)
- **Delivered:** §2.5 Results (three paragraphs) + Methods M.9 extension (retitled "…and
  environmental structure") + Supplementary Table S6 (`tab:darkfam-environment`) +
  fig:environment REBUILT on this basis (`code/make_darkfam_env_figures.py`; supersedes
  `make_darkmatter_figures.py` figs 1–2). Both PDFs compile clean; S6 resolves; snapshot mirror
  synced. Registered: CORRECTIONS 2026-08-16 entry + closed block, CATALOGUE_NUMBERS §2.5 table
  + superseded-map row, PLAN, PAPER_STRUCTURE, DISCUSSION_NOTES deposits.
- **Open (small):** McArdle/Anderson PERMANOVA bibkey = \todo in M.9 (citations session);
  plant-host genus axes unrunnable until the amplicon cache reaches IBEX; plantroot-contrast
  result (60,313/102,005) computed but uncited by design (compartment biology = §2.3);
  **B3-scaffold contradiction now live in the PDF** — the unfenced scaffold paragraph still says
  "all of B3 concerns the same 23,382 dark orthogroups" and cites fig:environment, whose rebuilt
  caption says 106,042; this session created the concrete forcing reason for the fencing decision
  parked with Alejandra; **§2.5 opening re-check when C1 lands** — "the majority of catalogued
  proteins carry no KO" leans on §2.1's (historical 56.7%) claim; re-read the opening sentence
  when the PI's PICRUSt rebuild replaces that number.
- **SETTLED — do not re-litigate:** June ecology dirs untouched (banner-marking = PI item 2b);
  env_matrix.tsv reused as-is (§2.3 precedent); EC-decode scope-out (sibling annotator paper);
  Mo non-reproduction recorded, not a correction (never cited).

## Compute record (job 50575481, 3m53s)
- Purpose: rerun the June 2026 gene–environment analysis (B3 / §2.5) on the rebuilt dark-family
  universe (106,042 families × 551 samples), replacing every number computed on the superseded
  23,382-orthogroup June object.
- **Community level (547 samples with env; 105,977 families at prev≥5; Hellinger+Bray-Curtis,
  adonis 999 perms):** geochem R²=0.123 / beyond-space R²=0.080; climate beyond-space 0.030;
  compartment 0.044; space 0.079 — all p=0.001; vanadium block NULL (p=0.780); distance-decay
  Mantel r=+0.341 p=0.001. Core-stratified (≥25% prev) same picture (geochem|space 0.083).
- **Family level (XRF subset n=251, 94,568 testable at ≥25-present; BH per axis, q<0.05):**
  S 54,227 (57.3% of tested; 91.7% NEGATIVE) / Ca 47,174 (92.5% neg) / Cl 46,972 (90.4% neg) /
  salinity_proxy 37,280 (86.9% neg) / Fe 24,110 (95.0% neg) / K 11,935 / Na 4,603;
  aridity_index (n=456, 100,293 tested) 23,031 (95.9% neg); plantroot contrast 60,313/102,005
  (67.1% neg). NULLS: V, Mo, Ni, Co, Mn, Zn = 0 (Cu = 1).
- **Space-partial retention (control ranked lat+lon, q<0.05 surviving):** aridity 85.6% /
  Cl 82.6% / S 82.0% / Ca 71.1% / salinity_proxy 41.4% (aggregate proxy = most spatially
  confounded, as in June).
- **Direction-stability vs June (superseded basis):** axis ranking, ~90%-negative filter
  pattern, V/Ni/Co/Mn/Zn nulls and community-level ordering all reproduce; every PERMANOVA
  block strengthens (geochem 0.066→0.123, distance-decay 0.161→0.341). CHANGE: June's Mo
  signal (177 families, 100% positive) does NOT reproduce (0 on the rebuilt universe) — was
  never a paper claim; recorded here.
- **Verification:** eco_stats selftest PASS pre-submit; top-Ca-family Spearman independently
  reproduced with scipy.stats.spearmanr (r=-0.525 exact match); eco_1/eco_2 marginal counts
  cross-agree (aridity 23,031, salinity 37,280 identical in both scripts).
- June headlines for reference (superseded 23,382 basis): Ca 4,066 / S 3,964 / Cl 3,594 /
  salinity 3,179, ~90% negative; counts here ~10× larger — universe 4.5× bigger and
  median family prevalence 98/551 puts most families over the ≥25-present power gate.

## Method lineage — reuse, not reinvention
Scripts are copies of the June analysis in the paper repo `code/`, repointed to the rebuilt
basis; statistical method verbatim:
- `eco_0_depth.py` (new 2026-09-20, no June ancestor): per-sample sequencing depth from the
  fastp `*_report.json` files in `EmptyQuarter_Data/soil/wgs_assemblies/processed_reads/`
  (`summary.after_filtering.total_bases`). Nothing in the workbench metadata carried depth, so
  the covariate `eco_3` conditions on is built here; idempotent, rewrites the table each run.
- `eco_1_indicators.py` ← `phase_e_indicators.py` — per-family Spearman vs each axis
  (XRF-gated 13 elements + salinity_proxy; aridity on the climate subset; plant-root
  point-biserial), family tested if nonzero in ≥25 subset samples, BH-FDR per axis.
- `eco_2_partial.py` ← `phase_e_partial_osmolyte.py` part (a) — partial Spearman controlling
  ranked lat+lon for the 5 master axes (salinity_proxy, S, Cl, Ca, aridity_index) and, added
  2026-09-20, the plant-root contrast on the same footing.
- `eco_3_community.py` ← `phase_b_ecology.py` — Hellinger + Bray-Curtis; distance-based
  PERMANOVA (adonis, 999 perms) per block (compartment / space / climate / geochem, the
  latter two also |space); prevalence-stratified core check; Borcard varpart (499 perms);
  distance-decay Mantel. `eco_stats.py` = the June library plus the 2026-09-20 `strata` argument
  and tie tolerance (STATUS block); every free-permutation result is unchanged.
- `eco_4_compartment_families.py` (new 2026-09-20, no June ancestor): joins the persisted
  plant-root axis to phylum (`fam_phylum.tsv`, single-phylum families) and to expedition
  persistence (`fam_presence.tsv` x `sample_site_trip.tsv`, the rule behind the 79.4% baseline).
  It computes no new statistics: it reads `plantroot_vs_other_all.tsv` from `eco_1`.

## Deviations from the June scripts (each forced; no method change)
1. **EC-decode layers dropped** (EC join in eco_1; the whole osmolyte vignette, part (b) of the
   June partial script): no `family_ec.tsv` exists for the 106,042 basis and per-protein
   decoding is deferred to the dark-matter-annotator sibling paper (locked structure).
   Consequence: §2.5 makes no claims about what structured families do.
2. **Plant-host genus axes dropped** (host_Tribulus/Haloxylon/Calligonum): input
   `plants_by_site_trip.tsv` lives in the amplicon repo cache, which is not on IBEX (same
   blocker as §2.2's 16S cross-check). The plant-root compartment contrast IS retained (it
   needs only env_matrix). Panel: 19 June axes → 16 (13 XRF + salinity + aridity + plantroot).
3. **Machine-readable dumps added**: `indicators_summary.tsv`, `partial_summary.tsv`,
   `permanova_blocks.tsv`, `varpart.tsv`, `pcoa_coords.tsv`/`pcoa_explained.tsv` (June wrote
   only stdout/markdown/txt); the PCoA dump feeds the rebuilt fig:environment without
   rerunning permutations.

## Inputs (verified 2026-08-16)
- abundance: `../darkfam_basis_clean/results/fam_abundance.tsv` — 106,042 families × 551
  samples, member-aware CoverM TPM (matrix complete 551/551 after F13PRr1;
  md5 f7cb8c3e939ccb573838704f43604cf1). Key `family_rep` = `sample::gene` tokens.
- environment: repo `results/darkmatter_ecology/env_matrix.tsv` — 547 samples × 45 cols
  (compartment PlantRoot 266 / Deep 165 / Surface 116; has_xrf=1 for 251; climate ≤456).
  SAME table the shipped §2.3 gradients used (osmolyte_clean osmo_5/6/7) — reuse is
  consistency. Built 2026-06-28 from rob's local sources (not rebuildable on IBEX); md5
  4d6eafad897818075a0fc6ae32bdf754.
- joins: `../darkfam_basis_clean/results/fam_recurrence.tsv` (n_members/n_samples),
  `../darkfam_basis_clean/results/selection/dnds.tsv` (median_dNdS; contains dark+ctrl rows,
  keyed join is collision-free).

## Do not touch
- rob's June object: workbench `../darkfam_ecology/` AND repo `results/darkmatter_ecology/`
  outputs (INDICATORS.md etc.) stay as-is — banner-marking is PI-coordination item 2b.
  This folder only READS env_matrix.tsv from there.

## Job log
- eco_run (eco_1+2+3): **50575481** (submitted 2026-08-16; compile + eco_stats selftest passed pre-submit)
- eco_run (eco_1+2+3+4): **52135803** (2026-09-20, 3m13s; full per-axis tables, plant-root partial, eco_4)
- eco_run (eco_0+1+2+3+4): **52137904** (2026-09-20, 4m04s; within-site permutation + depth
  covariate; 47 of 50 result files md5-identical to 52135803)
