> **⚠ SUPERSEDED (2026-08-03).** Pre-pivot / pre-merge plan, kept for history. Current plan → `../notes/PLAN.md`; current narrative → `../notes/PAPER_STRUCTURE.md`; current numbers → `../notes/CATALOGUE_NUMBERS.md`.

# Manuscript outline — Functional dark matter of the Empty Quarter

> **Historical planning record.** Status and headline language below may be superseded. Use
> `manuscript/main.tex`, `handover.md`, and `reproducibility/claims.tsv` for current evidence.

**Status:** outline + key findings (documentation scaffold). Prose to be written. All numbers below are
final results from the computational study (see `docs/dark-matter-pipeline.md`, `plans/darkmatter-proteins.md`).
**Target:** Nature-family specialist (Nat Microbiol / Nat Ecol Evol / ISME J / Microbiome / Genome Biology).
Nature-main contingent on Phase F wet-lab validation (not yet done).

---

## Working title
*"The functional dark proteome of Earth's largest sand desert is real, environmentally structured, and
largely endemic."*
(Alt: "A real, selection-constrained dark proteome structures the microbiome of the Rub' al-Khali.")

## One-line significance
The uncharacterised protein families of the world's largest sand desert are not random novelty: they are
bona-fide, selection-constrained genes, reproducibly organised by environment, and a third are unknown to
the entire global metagenome protein universe.

## Abstract (key numbers to include)
- 551 soil metagenomes, 5,438 dRep MAGs → **23,382 multi-member dark protein families** (rep sequence-novel
  vs UniRef50, ≥3 members) as the analysis unit; from 563,619 sequence-dark proteins / 361,719 triple-negative.
- **99.3% under purifying selection** (median dN/dS 0.39; control non-dark 0.29) → real functional genes.
- **Environment-structured beyond geography** (geochemistry R²=5.2%, climate 3.3%, PERMANOVA p=0.001).
- **Pan-desert stable core** (45.6% of families across ≥30/60 sites, all 5 trips, 3+ yr) + **endemic
  specialists more constrained** (dN/dS 0.33 vs 0.42, p=5e-69).
- **36.6% absent from UniRef50 AND MGnify90** (718M sequences) → genuinely desert-endemic; 63.4% divergent
  (median 57% id) cosmopolitan environmental dark matter.

## 1. Introduction
- Microbial dark matter is mapped at planetary scale (Pavlopoulos 2023 *Nature*; ESM Atlas 2023 *Science*;
  Durairaj 2023 *Nature*; AGNOSTOS Vanni 2022) but remains **functionally hollow**; function assignment has
  plateaued (FUGAsseM 2025; hard irreducible core).
- **Gap:** is dark matter random annotation backlog, or a real, ecologically-organised feature of a biome?
  Global catalogues lack the fine environmental gradients to test this.
- **System:** Rub' al-Khali — largest sand desert, hyper-arid, never genome-sampled at this scale; 551
  samples × 60 sites × 5 trips with geochemistry (XRF, incl. V), climate, depth/compartment.
- **This study:** treat the dark proteome as an ecological object — quantify its reality (selection),
  environmental structure, biogeography, global endemism, and decodability.

## 2. Results
**R1 — A pan-desert dark-proteome catalogue.** Define dark families (UniRef50-novel + ≥3 members) → 23,382;
abundance across 551 samples by read recruitment. *(Methods/numbers; no headline figure.)*

**R2 — The dark proteome is real: under purifying selection.** 99.3% of families dN/dS<1 (median 0.39),
vs size-matched non-dark control (0.29); dark slightly faster-evolving (expected — why homology-dark).
→ **Fig 1** (`fig_darkmatter_dnds`). *Refutes "spurious ORF" null.*

**R3 — Ecological structuring beyond geography.** PERMANOVA: geochem|space 5.2%, climate|space 3.3%
(p=0.001), robust to rare-family noise (core-family R² unchanged); distance-decay r=0.16. Modest but
reproducible. Vanadium-specific signal null (honest negative). → **Fig 2** (`fig_darkmatter_pcoa`,
`fig_darkmatter_permanova`).

**R4 — Biogeography: stable core + conserved endemic specialists.** 45.6% EQ-wide (≥30 sites, all trips);
site-restricted families more constrained (0.33 vs 0.42, p=5e-69; Spearman n_sites~dN/dS r=+0.18).
→ **Fig 3** (`fig_darkmatter_endemism_sites`, `fig_darkmatter_endemism`).

**R5 — A third is globally endemic.** vs MGnify90 (718M clusters): 36.6% no hit (endemic beyond all known
env protein space); 63.4% divergent (median 57% id), 58.5% cosmopolitan (≥3 biomes, mostly aquatic/engineered).
→ **Fig 4** (`fig_darkmatter_global_endemism`). *(Caveat: MGnify biome sampling bias.)*

**R6 — Decodable slice encodes plausible desert-survival functions (bounded).** C·V·T assigns EC to 23,083
families (calibrated P 0.22–0.43); desert-relevant categories present (oligotrophy, C1, osmolytes, dormancy,
oxidative stress) but NOT enriched among structured families and not aridity-tracking (honest negatives).
Supporting only.

**R7 — Individual dark families track geochemistry and plant host and converge on an osmolyte repertoire.**
Where R3/R6 worked at the community and category level, this resolves single families. (Code:
`phase_e_indicators.py`, `phase_e_partial_osmolyte.py`; tables in `results/darkmatter_ecology/indicators/`,
`INDICATORS.md`, `PARTIAL.md`, `OSMOLYTE.md`.)
- **R7a — per-family geochemical tracking survives spatial control.** Spearman of family abundance vs each
  XRF element, recomputed as partial correlation controlling lat+lon: sulfur 79%, chlorine 84%, calcium 72%,
  aridity 71% of significant families retained at q<0.05 → tracking is environmental, not geography (corroborates
  the R3 geochem|space result at single-family resolution). **Vanadium null reproduced** (0/11,436 families,
  marginal and partial) — the cleanest restatement of the honest V negative.
- **R7b — plant-host specificity.** Families enriched in the rhizosphere of specific dominant plants
  (iNaturalist-nearest host per site×trip): *Tribulus* → glucosylglycerol-phosphate phosphatase (EC 3.1.3.69,
  r=+0.51, q=5e-4) + alkaline phosphatase (P-scavenging); *Calligonum* → glutamate–cysteine ligase
  (glutathione / oxidative-stress, EC 6.3.2.2); *Haloxylon* null. A genome×environment axis the global
  catalogues cannot address.
- **R7c — convergent osmolyte signal (bounded).** 306 dark families decode to osmolyte phosphatases —
  mannosyl-3-phosphoglycerate (EC 3.1.3.70; largest family 127 members / 85 samples, dN/dS 0.49),
  glucosylglycerol (3.1.3.69), trehalose-6-P (3.1.3.12) — all under purifying selection (dN/dS 0.33–0.66) and
  PlantRoot-abundance-skewed; space-partial salinity survivors add glucosyl-3-phosphoglycerate phosphatase
  (3.1.3.85), polyol dehydrogenases and the GlgE pathway. Multiple independent families converging on the same
  compatible-solute chemistry is stronger than any single calibrated-P EC. **Honest negative:** as a *class* the
  306 families are NOT plant-root-enriched beyond the dark background (one-sided MWU p=0.11) — the signal is
  carried by the large individual families, not the category (mirrors the R6 no-category-enrichment result).
- **R7d — who carries it (in flight).** Member→MAG→GTDB resolution for all 23,382 families (IBEX job 47715706,
  `code/ibex/99_darkfam_member_taxa.{sh,py}`) → taxonomic breadth per family + which taxa encode the desert
  osmolyte dark matter. Figure asset pending; rep-MAG taxa already include candidate genera (e.g. JAJPHT01).
→ **Fig 5** (`fig_darkmatter_indicators`: space-partial element bars + host panel; `fig_darkmatter_osmolyte`:
osmolyte-family size/selection/compartment + taxa once R7d lands).

**R8 — Intrinsic resistome and biosafety profile: broad stress tolerance, no pathogen bloom.** A host-genome
trait layer that doubles as reviewer-proofing (contamination control) and as ecology (geochemical stress
adaptation). Worked at the MAG level on the full 15,469-MAG set (+185 isolates / 118 enrichment bins). (Code:
`phase_amr_metal.py`, `pathogen_census.py`, `vfdb_analyze.py`, `code/ibex/100_vfdb_census.sh`; tables in
`results/{amr_metal,pathogens,vfdb}/`.)
- **R8a — broad intrinsic antibiotic + metal resistome in a low-anthropogenic biome.** AMRFinderPlus (Type=AMR):
  β-lactam in **3,816 MAGs (~25%)**, aminoglycoside 2,002, phenicol 1,108, rifamycin 1,086, macrolide 705,
  glycopeptide 622 — a broad ancient/intrinsic resistome with minimal human input. The STRESS layer (metal /
  metalloid) adds **4,159 hits in 3,187 MAGs** (+581 in isolates/enrichment): **arsenic 2,423 MAGs**, **mercury
  1,077**, dominated by Pseudomonadota (As) and Actinomycetota (Hg). Arsenic-resistance load tracks the salinity
  proxy across XRF samples (Spearman r=+0.26, p=3e-5; Hg r=+0.16, p=0.012) — a geochemical-stress signal, not a
  cognate-metal one (the XRF panel has **no As or Hg**, stated as a limit).
- **R8b — pathogen census = clean contamination control.** Of 15,469 MAGs, genera that *contain* pathogens
  account for 126 human/animal (0.8%) and 9 plant (0.1%) MAGs — overwhelmingly environmental relatives
  (*Mycobacterium*, *Nocardia*, *Rhodococcus*, soil *Bacillus*). Actual **named pathogen species: 7 MAGs (0.05%)**
  — 5 *Escherichia coli*, 2 *Acinetobacter baumannii*. The genus→species collapse is exactly the signature of a
  pristine, low-anthropogenic biome with **no kit-ome / human-contamination bloom** (reviewer-proofing).
- **R8c — virulence "arsenal" is a dual-use competition/nutrient toolkit, not pathogenesis.** VFDB setB
  (mmseqs, id≥40/cov≥70/e≤1e-10) hits 541,508 proteins — but the count is inflated by dual-use genes (secretion,
  iron/heme uptake, adhesins, LPS) ubiquitous in soil. The signal collapses with stringency: MAGs with ≥1 VF hit
  fall from 5,435/5,438 (id≥40) to 3,769 (id≥70) to **202 (id≥90)**; the pathogenesis-specific **Exotoxin**
  category collapses from 6,260 proteins to **16 at id≥90** (scattered hemolysin-III / DtxR / phospholipase-D
  homologs). **Honest framing:** what VFDB scores as "virulence" in this biome is competition and nutrient
  acquisition, not human virulence.
→ **Extended Data** (intrinsic-resistome + metal census by phylum/compartment; VFDB identity-tier collapse;
pathogen genus→species gap). No main-text figure required; this layer is safety/control + a stress-ecology hook.

## 3. Discussion
- Reframes dark matter from "annotation backlog" to a **real, selection-maintained, ecologically-organised**
  component of an extreme-biome microbiome.
- Structuring is significant but **modest** — most variance environment-independent (stochastic/dispersal/
  unmeasured); honest about magnitude.
- Endemic third = a genuinely novel desert protein reservoir; conserved specialists suggest specialised
  (still-unknown) functions worth targeted study.
- Limits: capacity not flux; function largely undetermined; "adaptive" not demonstrated (negatives in R6).
- **Next (Phase F):** metatranscriptomics on the 13 dual-marker strains + heterologous expression of the
  166 soluble decoded targets → demonstrate function (Nature-main lever).

## 4. Methods (→ `docs/dark-matter-pipeline.md`, `code/` + `code/ibex/`)
Catalogue + CoverM abundance (90–92); env matrix (build_env_matrix); PERMANOVA/varpart/Mantel/PCoA
(eco_stats + phase_b); member-CDS + NG86 dN/dS (93–96 + phase_c); within-EQ + MGnify90 endemism
(phase_d + 98a–c); C·V·T decode bridge (97 + phase_e); per-family indicators + space-partial (lat/lon)
correlation + osmolyte vignette (phase_e_indicators, phase_e_partial_osmolyte); member→MAG→GTDB taxonomy
(`code/ibex/99_darkfam_member_taxa`). Intrinsic resistome / biosafety (R8): AMRFinderPlus AMR+STRESS census and
XRF link (phase_amr_metal); GTDB genus/species pathogen census (pathogen_census); VFDB setB mmseqs + identity-tier
dual-use analysis (vfdb_analyze, `code/ibex/100_vfdb_census`). Reproducible from this folder.

## 5. Figures (in `figures/`)
1. `fig_darkmatter_dnds` — purifying selection (dark vs control). 2. `fig_darkmatter_pcoa` +
`fig_darkmatter_permanova` — env structuring. 3. `fig_darkmatter_endemism_sites` (+ `_endemism`) —
biogeography/selection. 4. `fig_darkmatter_global_endemism` — MGnify endemism. (Extended data: composition,
variance partitioning.) 5. `fig_darkmatter_indicators` (space-partial element associations + plant-host panel)
+ `fig_darkmatter_osmolyte` (osmolyte-family size/selection/compartment + carrier taxa, pending R7d).

## 6. Honest negatives (state explicitly, prevents overclaiming)
structured ≠ adaptive (dN/dS 0.389 vs 0.392); decode↔ecology uncoupled at the category level; no aridity
tracking; vanadium null (reproduced at single-family resolution, 0 families marginal or space-partial);
osmolyte families not enriched as a *class* (MWU p=0.11 — signal is in large individual families, not the
category); effect sizes modest. The VFDB "virulence" census is dual-use, not pathogenesis — it collapses to ~16
genuine exotoxin homologs at id≥90 and must not be read as a virulence claim; the As/Hg–XRF link is
non-cognate (XRF lacks As and Hg). These bound the claims and set the wet-lab agenda.
