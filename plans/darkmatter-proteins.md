> **⚠ SUPERSEDED (2026-08-03).** Pre-pivot / pre-merge plan, kept for history. Current plan → `../notes/PLAN.md`; current narrative → `../notes/PAPER_STRUCTURE.md`; current numbers → `../notes/CATALOGUE_NUMBERS.md`.

# Dark-matter proteins — analysis & paper plan

> **Historical planning record.** Status and headline language below may be superseded. Use
> `manuscript/main.tex`, `handover.md`, and `reproducibility/claims.tsv` for current evidence.

Status as of 2026-06-19. Dedicated plan for the **Pillar 1** dark-matter story (see
`analysis-roadmap.md`). Heavy steps run in `../workbench` on IBEX via `sbatch` (`hohndor_*`,
`${EQ_WORKBENCH}`, no hardcoded home paths). Artifacts land in this repo (`results/`, `paper/`).
Scope: 15,469 soil MAGs (5,438 dRep reps) + 185 isolates + 1,036 enrichment MAGs.

## Context

We can characterise only a small fraction of EQ proteins, and a smaller fraction by function. The
instinct — "make the story by *solving* function assignment" — is a trap. Function assignment has
**plateaued**: the gut version was already cracked (FUGAsseM, *Nat Biotechnol* 2025, with 800
metatranscriptomes + context), and a hard core (>50% of genuine novelty) is **irreducibly
unannotatable** (Vanni/AGNOSTOS *eLife* 2022: "sequence similarity appears to have yielded all it
can"). Our own calibrated C·V·T pipeline confirms it — for truly dark proteins P(correct) ≈
**0.22–0.43**. A "we solved annotation" paper would be untrue and rejected.

**Reframe (chosen headline):** stop trying to *name* the dark proteins; prove the dark proteome is
**ecologically deterministic and under selection** — an environmentally-organised, adaptive
repertoire of hyper-aridity, not a random annotation backlog. This is the one frontier the global
catalogues (Pavlopoulos *Nature* 2023; ESM Atlas *Science* 2023; Durairaj *Nature* 2023) could
**not** touch — they lack fine environmental gradients. We are uniquely positioned: 551 samples ×
60 sites × 5 trips, 59-element XRF (incl. vanadium), SMAP moisture, NDVI, climate + cultures.

**Decisions (PI):**
- **Scope:** write the full computational paper now for a **Nature-family** journal; launch the
  cheapest wet-lab step (metatranscriptomics on the 13 dual-marker strains) in parallel as the
  **upgrade path** to Nature/Science-main. Do **not** gate submission on wet lab.
- **Headline (revised 2026-06-22 after Phase C):** *"The functional dark matter of Earth's largest
  sand desert is real — 99% under purifying selection — and reproducibly structured by environment
  beyond geography, not random novelty."* The "adaptive" claim was DROPPED: structured vs cosmopolitan
  dark families show negligible dN/dS difference (0.389 vs 0.392), so env-responsiveness is not tied to
  stronger selection. Target = Nature-family specialist; Nature-main contingent on Phase F validation.
  The function-assignment pipeline is **supporting** evidence (decodable slice skews osmolyte/stress/C1).

## What already exists (do NOT rebuild)

- **Dark-matter universe** (`workbench/results/darkmatter_function_assignment/`): 5,954,339 dark
  candidate proteins / 2,172,743 families from 5,438 dRep reps; **563,619 sequence-novel** (no
  UniRef50 / <30%); **361,719 triple-negative** (no homolog, no InterPro domain, no AFDB fold).
- **C·V·T function-assignment engine — DONE & calibrated:** gapsmith gaps (1.07M reactions) +
  genomic context + ProstT5/FoldSeek + InterPro + SignalP + antiSMASH + pFBA/FVA. EC-recovery 63%
  exact / 92% class-level on 1.9M held-out genes; `reliability_curve.tsv` maps score→P(correct);
  top-100 reviewed (22 well-supported EC, 35 plausible, 43 pathway-assoc) + opportunity-scored.
  Decoded hits already skew desert-adaptive (osmolyte phosphatases EC 3.1.3.69/3.1.3.70,
  pppGpp-hydrolase, C1/formate-THF ligase).
- **Structural novelty:** ProstT5+FoldSeek vs AFDB done; tier1/tier2 + 9,594 really-dark novel-fold
  families.
- **Ecology axes (ready):** `workbench/data/metadata/sample_metadata.tsv`,
  `workbench/data/geodata/trip{1..5}_geodata.tsv`,
  `workbench/data/geochemistry/xrf_lab_table_filtered.tsv` (59 elements incl. V/V2O5, Na/Cl/S),
  amplicon `../empty-quarter-amplicon/repository/cache/{smap_by_sample,ndvi_by_sample,climate_historical_1995_2024}.parquet`.
- **Cultures (annotated):** 185 isolates (GTDB: 98 novel sp., 16 genera, 2 families, 1 class-level),
  1,036 enrichment MAGs. 13 dual CO+H₂ strains + 166 soluble dark targets in
  `workbench/results/unimatrix_halophilic/experimentalist_targets.tsv`.
- **NOT present:** any metatranscriptomics / proteomics / RNA-seq — must be generated for validation.

## Research Questions (literature-verified open status)

- **RQ1 — Foundational (supporting).** Size/novelty of the hot-sand-desert dark proteome after
  fragment control, and how much is decodable vs irreducibly unknown. *Open:* biome-specific yes;
  the honest **decodability ceiling** is mildly novel. Mostly write-up + a "decodability curve" fig.
- **RQ2 — HEADLINE.** Is the dark proteome ecologically structured along aridity/moisture,
  depth/compartment, geochemistry (V, salinity), temperature; are families endemic to extremes?
  *Open:* **YES — completely unexplored.** Lead result.
- **RQ3 — HEADLINE support ("adaptive").** Are recurrent dark families under purifying selection
  (dN/dS « 1), strongest in the structured/endemic families? *Open:* dN/dS established;
  desert-endemic application novel. Converts "structured" → "adaptive."
- **RQ4 — Positioning.** EQ/desert-endemic vs globally shared (overlap with Pavlopoulos / AGNOSTOS
  GU·EU / AFDB·ESM clusters). *Open:* per-biome aridity endemism unexplored; data-access risk.
- **RQ5 — Decoding↔ecology bridge (supporting).** Is the decodable slice enriched for
  hyper-aridity-adaptive functions (osmolyte, stress, C1, metal/V) mapping onto RQ2 niches? *Open:*
  novel in this framing; doesn't depend on "solving" annotation.
- **RQ6 — Validation / Nature-main upgrade (parallel).** Are dark families expressed; can ≥1
  function be experimentally validated? *Open:* desert dark-matter expression unexplored; one
  validated function pushes to Nature-main. Wet-lab timeline — not a gate.

## Analysis workplan

### Phase A — Dark-family × sample abundance matrix (enables RQ2/RQ4/RQ5)
1. Unit = dark orthogroups (reuse `scores_v3` + `orthogroups_cluster.tsv`); keep multi-member
   families (≥3 members), drop singletons.
2. Per-sample abundance by **read recruitment** (CoverM; `metabolism-env` has CoverM) of each
   sample's reads to dark-family reps — captures rare/unbinned families. Reuse the amplicon CoverM
   approach if one matches; else a 551-sample IBEX array.
3. Output: `darkfam_abundance.tsv` (family × sample) + presence/absence companion.

### Phase B — RQ2 ecological structuring (HEADLINE)
1. Build per-sample environment matrix (aridity from precip/PET or SMAP; V + salinity from XRF;
   compartment; trip; lat/lon; T; NDVI).
2. **Constrained ordination** (db-RDA/CCA) + **variance partitioning** (geochemistry vs climate vs
   compartment vs space) with PERMANOVA/Mantel and a distance-decay control to separate environment
   from geography.
3. **Indicator / differential-abundance** families per niche (extreme- vs less-arid; Surface/Deep/
   PlantRoot; high- vs low-V); per-family **endemism/specialisation index**.
4. Deliverables: ordination + variance-partition figures; ranked niche-indicator table.
   Gate: environment explains composition above the spatial null (PERMANOVA p<0.05).

### Phase C — RQ3 selection signatures ("adaptive" claim)
1. Per-family codon alignments from member CDS (Bakta/pyrodigal); compute **dN/dS** (HyPhy/codeml or
   fast pairwise pN/pS) on families with ≥N members + adequate divergence.
2. Test: structured/endemic families (Phase B) under **stronger purifying selection** than
   cosmopolitan dark families and a matched random-gene null; tie to recurrence.
3. Deliverable: dN/dS by family class + the endemic×selected 2×2.
   Gate: structured/endemic families dN/dS « 1 and below controls. If NOT → **report the negative**
   and soften headline to "ecologically structured" (no overclaim).

### Phase D — RQ4 endemism vs global catalogues
1. Search EQ dark-family reps vs Pavlopoulos/MGnify, AGNOSTOS GU·EU, AFDB/ESM-Atlas (3Di/FoldSeek)
   on IBEX; stage large DB pulls via `ws`.
2. Quantify EQ-specific fraction; overlay endemism onto Phase-B structure.
3. Deliverable: endemism census + shared-vs-EQ-specific figure. Fallback: report whatever subset is
   feasible and state the cap explicitly (no silent truncation).

### Phase E — RQ5 decoding↔ecology bridge (supporting)
1. Reuse C·V·T assignments + reliability curve; functional-enrichment of confident EC/pathway
   assignments vs background for desert-adaptive categories (osmolyte, stress, C1, metal/V).
2. Spatially overlay decoded functions onto Phase-B niches; test V-associated families vs XRF
   vanadium. Deliverable: enrichment table + corroboration figure (bounded by calibrated P).

### Phase F — RQ6 validation (PARALLEL, Nature-main upgrade)
1. In-silico now: map structured/selected dark families into 185 isolates + 1,036 enrichment MAGs →
   culturable, tractable shortlist (intersect with the 166 soluble targets).
2. Commission **metatranscriptomics on the 13 dual-marker strains** under desert stress
   (desiccation/osmotic/oligotrophy); coordinate with unimatrix/culture owners.
3. Optional: heterologous expression + activity assay on 1–2 decoded osmolyte/stress candidates.
   Decoupled from submission timing.

## Honest journal verdict
- **Computational paper (RQ1–RQ5):** first ecological-structuring + selection analysis of a
  hot-sand-desert dark proteome; reframes dark matter as adaptive signal. Target **Nature
  Microbiology / Nature Ecology & Evolution**; fallbacks **ISME J / Microbiome / Genome Biology**.
  This is the explicit, defensible downgrade from Nature-main — a real, high-impact paper.
- **+ RQ6 (≥1 validated function + conditional expression):** competitive for **Nature / Science**
  main, at 6–18-month wet-lab cost — hence the parallel track.
- **Bottom line:** worth doing, worth a strong Nature-family paper now; Nature-main is an upgrade
  contingent on validation, not a precondition.

## Verification / success criteria
- **A:** matrix covers all 551 samples; counts reconcile with the OG catalogue; spot-check 5
  families' recruitment vs MAG occurrence.
- **B (headline gate):** PERMANOVA p<0.05 above spatial null; ≥1 axis (aridity or geochemistry) with
  substantial interpretable variance; reproducible indicator list.
- **C (adaptive gate):** structured/endemic dN/dS significantly <1 and below controls — else soften.
- **D:** defensible EQ-specific fraction with comparison scope stated.
- **E:** desert-adaptive categories significantly enriched vs background.
- **F:** target shortlist produced; metatranscriptomics commissioned (deliverable = expression
  matrix, stress vs control).
- **Reproducibility:** per-step counts/params/job-IDs/DB-versions into `plans/analysis-roadmap.md` +
  results `README`s as each step completes.

## Progress log

**2026-06-21 — Phase B preliminary (545/551 samples; HEADLINE).** `code/phase_b_ecology.py` +
`code/eco_stats.py` (adonis optimised to O(n²)/perm via trace(HGH)=sum(G*H)). Dark-family composition
(21,906 families prev≥5, Hellinger→Bray-Curtis) vs environment, PERMANOVA 999 perms:
- compartment R²=2.0% (p=.001), space(lat/lon) R²=3.8% (p=.001), climate R²=4.7% (p=.001),
  geochem R²=6.6% (p=.001).
- **Beyond geography (the gate): geochem|space R²=5.2% (p=.001), climate|space R²=3.2% (p=.001)** →
  environment structures the dark proteome above the spatial/distance-decay null. **Headline HOLDS.**
- Distance-decay Mantel r=0.163 (p=.001). **Vanadium alone NOT significant** (R²=0.4%, p=0.352) — the
  V-specific hypothesis does not appear in dark-proteome composition.
- **Core-family stratification (prev≥25%): R² essentially unchanged** (geochem|space 5.5%, climate|space
  3.6%, compartment 2.0%) → the modest effect is NOT rare-family noise; it is the genuine magnitude.
- ⚠️ Honest verdict: structuring is highly significant + robust (to geography AND to rare-noise) but
  **modest in magnitude** (~3–5% beyond space). Frame as "significantly & reproducibly structured beyond
  geography," NOT "deterministic" — most dark-proteome variance is environment-independent
  (stochastic/dispersal/unmeasured). This R² is in the normal range for soil community beta-diversity.
  The Nature-family-vs-main calculus now hinges on **Phase C (selection ⇒ adaptive)** and Phase D
  (endemism), not on Phase B effect size alone. Complete 551 matrix won't move these numbers materially.
- Recovery of the 6 straggler samples in flight: 3 pathological (low-complexity reads, stuck in
  minimap2 >14h despite normal file size) → subsample job 47692966; 3 timed-out → rerun 47692989;
  final aggregator 47692990 → complete 551 matrix. Interim matrix pulled to
  `results/darkmatter_ecology/darkfam_abundance.tsv`.

**2026-06-21 — Phase A complete + Phase C selection (pivotal).**
- Phase A: **complete 551/551-sample matrix** (23,382 families; 6 stragglers recovered — 3 pathological
  low-complexity samples via subsample, 3 via 24 h rerun). `darkfam_abundance.tsv` final.
- Phase C (`93/94/96`, scripts member-CDS → mafft+NG86 dN/dS → test): **22,767 dark + 21,579 size-matched
  control families**. **DARK median dN/dS = 0.391, 99.3% < 1 (purifying selection)** → the dark proteome is
  real, functional, selection-constrained genes, NOT spurious ORFs. Control (non-dark) median 0.286, 99.8%<1;
  dark vs control MWU p≈0 (dark slightly LESS constrained — expected, fast evolution is why they're
  homology-dark). Strong validation of the whole enterprise.
- Adaptive test (`code/phase_c_adaptive.py`): structured (9,425) vs cosmopolitan (13,342) dark families
  median dN/dS **0.389 vs 0.392** — significant (MWU p=0.01, large-n) but **negligible difference →
  the "structured = adaptive" lever does NOT fire.** Modest endemism signal: specialist/low-prevalence
  families more constrained (Spearman prevalence~dN/dS r=+0.158, p≈0).
- **Honest verdict:** dark proteome is REAL (99.3% purifying) + ecologically STRUCTURED (modest, beyond
  geography) + specialists slightly more conserved — a strong, defensible story, but the computational
  selection analysis does NOT deliver a Nature-main "adaptive" blockbuster. Target = **Nature-family
  specialist** (Nat Microbiol / Nat Ecol Evol / ISME J / Microbiome / Genome Biology). Nature-main would
  require the Phase F experimental validation (or a standout Phase E decoding discovery) to add a
  demonstrated functional finding.

**2026-06-22 — Phase D (HEAVY: global endemism vs MGnify90).** Scripts `98a`(download+DIAMOND makedb)
/`98b`(blastp)/`code/phase_d_global.py`. Searched all 23,382 dark reps (already UniRef50-dark) vs
**MGnify90** (`databases/mgnify90/`, 74 GB, ~718M clusters, the largest global env protein catalogue;
reference-area `/ibex/reference/EBI` is owner-only, c2014-staged copy was annotation-only → pulled from
EBI FTP). DIAMOND --sensitive, e≤1e-3.
- **36.6% (8,552) absent even from MGnify90 → genuinely desert/EQ-endemic** (no UniRef50 AND no MGnify hit).
- 63.4% (14,830) globally-shared but **divergent** (best-hit median 56.6% id) — related to uncharacterised
  env sequences elsewhere, not known proteins.
- endemic vs shared: dN/dS 0.393 vs 0.390 (negligible); endemic more locally restricted within EQ
  (median 41 vs 62 samples). Fig `fig_darkmatter_global_endemism`. Strongest positioning result.
- Biome breakdown (`98c` → `BIOME_BREAKDOWN.txt`) of the 14,830 shared families: 83.5% aquatic, 48%
  bioreactor, 48% plant-associated, 43% wastewater; **58.5% cosmopolitan (≥3 biomes)** = broadly
  distributed global env dark matter. Only 9.7% in other terrestrial soils — ⚠️ **confounded by MGnify
  biome bias** (aquatic/engineered/human-dominated; soil/desert underrepresented), so do NOT over-read
  the low soil fraction. Honest statement: the shared third connects EQ dark matter to global uncultured
  diversity; ~37% is endemic to this desert.
- **Phase D verdict:** 36.6% desert-endemic (absent from all known env protein space) + 63.4% divergent,
  cosmopolitan global env dark matter. Strongest positioning result; figure `fig_darkmatter_global_endemism`.

**2026-06-22 — Phase D (within-EQ biogeography / endemism; feasible form).** `code/phase_d_endemism.py`
→ `ENDEMISM.md` + `figures/fig_darkmatter_endemism_sites`. Global novelty already established by
construction (UniRef50-dark; really-dark also AFDB-fold-novel); full MGnify/GMGC/AGNOSTOS search deferred
(no catalogues on IBEX; heavy download, low expected payoff given weak amplifiers).
- **EQ-wide stable core: 45.6% of families occupy ≥30/60 sites (median 26 sites, median 5/5 trips)** — a
  pan-desert dark proteome, temporally stable 3+ yr. Only 2.3% single-site, 8.5% ≤3 sites.
- **Endemic families MORE constrained:** ≤3-site median dN/dS 0.332 vs EQ-wide 0.420 (MWU p=5e-69;
  Spearman n_sites~dN/dS r=+0.176, p=3e-157). Cleanest gradient signal — specialists under stronger
  purifying selection. (Adds a real supporting pillar beyond Phase B/E.)

**2026-06-22 — Phase E (decoding↔ecology bridge) + figures + headline reframe.**
- Decode join (`97`, on col1 candidate_protein): **23,382/23,382 families assigned, 23,083 specific EC**
  (best C·V·T per family). Top hits desert-relevant (catalase, dipicolinate synthase, aromatic dioxygenases).
- Bridge (`code/phase_e_bridge.py`): decodable slice contains desert-survival categories (aromatic/oligotrophy
  454, C1 443, metal 394, osmolyte 263, dormancy 221, oxidative-stress 79) — **descriptive positive**. BUT
  desert-adaptive categories are **NOT enriched** among structured families (all p>0.2), and do **NOT** track
  aridity in the expected direction (r≈0.1, mostly opposite). → decoding does not amplify the ecology story.
- Figures (`code/make_darkmatter_figures.py` → `figures/fig_darkmatter_{pcoa,permanova,dnds,endemism}.{pdf,png}`).
- **Cumulative honest verdict:** ONE strong result — dark proteome is REAL (99% purifying selection, Phase C).
  Phase B structuring is solid but MODEST (~5% beyond space). Every amplifier hypothesis came back
  weak/negative: adaptive-selection (no), functional-ecology coupling (no), aridity-tracking (no), vanadium (no).
  Fully-supported headline: *"real + selection-constrained + modestly environment-structured dark proteome
  whose decodable slice encodes plausible desert-survival functions."* Solid Nature-family-specialist paper;
  Nature-main needs Phase F wet-lab. Honest negatives prevent overclaiming.

## Progress log (earlier)

**2026-06-19 — Phase A launched.** Scripts `workbench/scripts/ibex/90_darkfam_catalog.sh`,
`91_darkfam_coverm_array.sh`, `92_darkfam_aggregate.sh`. IBEX workspace
`06_darkmatter/darkfam_ecology/`.
- Analysis unit = orthogroups whose representative is sequence-dark (`reallydark/isdark_ids.txt`,
  563,619 proteins) AND ≥3 members (from `orthogroups/og_recurrence.tsv`) → **23,382 dark
  multi-member families**.
- Nucleotide gene catalogue: CDS extracted from dereplicated genome FASTAs
  (`05_refine/_drep/out/dereplicated_genomes/<mag>.fa`) via pyrodigal GFF coords (genes were
  pyrodigal-called, so Bakta `.ffn` IDs don't match). Build job 47654332 (2 min): **23,382/23,382
  CDS extracted, 0 missing, 14.5 Mbp**.
- Manifest restricted to **551 canonical soil samples** (dropped 19: 16 C-prefix
  separate-experiment + `48PRR1`,`T8DDr2`,`T8DDr3`).
- Recruitment = CoverM `contig` (minimap2-sr), per-gene TPM + covered_fraction; presence at
  covered_fraction ≥ 0.5. Smoke test (10Dr2, job 47654397, 15 min, ~29 GB): 932 present / 3,453
  TPM>0; 0.29% reads recruit to the 14.5 Mbp catalogue.
- Full array 47654618 (550 tasks, %30, 16 CPU/64 GB/8 h) + aggregator 47654619 (`afterany`) →
  `darkfam_abundance.tsv` + `darkfam_presence.tsv` (23,382 × 551). Running.

**2026-06-20 — Canonical placement (completion-gated).** Builder `workbench/scripts/ibex/95_canonical_links.sh`
(job 47669483; symlinks only, non-destructive, idempotent) publishes **only COMPLETE analyses** into
`EmptyQuarter_Data/processed/sequencing-results/catalog/` — the PI's results-pickup point, not a working dir.
Each section is gated behind a completion check and a `MANIFEST.tsv` records complete vs excluded.
- Keystone: `catalog/mags/dereplicated/<mag>.fa.gz` → per-sample bin (5,394/5,438; 44 → dRep genome
  fallback; 1,010 DAS_Tool `_sub` reps → parent SemiBin) + `catalog_membership.tsv`.
- **COMPLETE & linked:** mags, proteins(+orthogroups), taxonomy (GTDB-Tk: 14,752 bac + 717 arc = all
  15,469 MAGs), structure (final novelty figures + AFDB hits only, not working chunks),
  dark_matter/function_assignment, dark_matter/families (419,748 → MAG provenance),
  analyses/atmospheric_guild (10 DONE stages), `../phages/` (geNomad/PhaBox/VIBRANT/vOTU, 550/550 + vOTU).
- **EXCLUDED (in-progress/not-started):** dark_matter/ecology (Phase A running → auto-links on `AGG_DONE`),
  ESMFold Phase 4c (partial), picrust/nitrogen/dormancy/kraken2 (empty).
- Convention saved to memory; per-sample phage anchoring into `wgs/phage/` is a flagged follow-up.

**2026-06-19 — Phase B prep (environment matrix).** `code/build_env_matrix.py` →
`results/darkmatter_ecology/env_matrix.tsv` (547/551 samples; 4 dropped = sites 62/64 absent from
geodata). Coverage by axis: **climate/aridity 456** (precip 6.6–74.7 mm, aridity index 0.17–1.84),
**geochemistry/XRF 251** (V 0–0.071%, salinity proxy Na+Cl+S 0.1–14.2), compartment/trip/space ~547.
⚠️ The workbench-symlinked `data/geochemistry/xrf_lab_table_filtered.tsv` is **trip-5 only** (178
rows, all V-prefix); the real source is amplicon `repository/data/geochemistry/xrf_lab_table_all_trips.tsv`
(725 rows, all trips) → 251 match our samples. Phase B design: primary ordination on the full ~547
(climate+compartment+space), geochemistry/vanadium axis on the 251-sample XRF subset. SMAP/NDVI
(parquet) deferred — needs pyarrow (not installed locally).
