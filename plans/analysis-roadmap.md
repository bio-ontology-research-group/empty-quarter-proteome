> **⚠ SUPERSEDED (2026-08-03).** Pre-pivot / pre-merge plan, kept for history. Current plan → `../notes/PLAN.md`; current narrative → `../notes/PAPER_STRUCTURE.md`; current numbers → `../notes/CATALOGUE_NUMBERS.md`.

# Analysis roadmap — Empty Quarter MAGs paper

> **Historical planning record.** Status and headline language below may be superseded. Use
> `manuscript/main.tex`, `handover.md`, and `reproducibility/claims.tsv` for current evidence.

Status as of 2026-06-15. Operational companion to `../CLAUDE.md`. Heavy steps run in the
workbench (`../workbench`) on IBEX/unimatrix; this file tracks what's done and what's next
for the three MAGs-paper pillars. Scope for every new analysis: **15,469 soil MAGs + 185
isolates + 118 enrichment cultures.**

## Already done (workbench)
- MAG catalog (15,469 HQ+MQ), dRep (3 variants), DAS_Tool refinement.
- GTDB-Tk r226 → **no novel phyla/classes** (the pivot; see CLAUDE.md).
- Bakta + EggNOG on all MAGs; mmseqs reps (23.4M); UniRef50 search (12% no-hit).
- ProstT5 + FoldSeek vs AFDB → ~2.5M structurally-novel reps; novel-fold clustering.
- Dark-matter function assignment (C·V·T, LRO-calibrated). Tier2 FoldSeek rerun = IBEX job
  47540824 (78 missing chunks → 100% coverage).

## Pillar 1 — Functional & structural dark matter
> **Lead story → dedicated plan: `plans/darkmatter-proteins.md`** (2026-06-19). Headline reframe:
> the dark proteome is **ecologically structured & under selection** (adaptive, not random) — NOT
> "we solved function assignment" (that has plateaued; our calibrated P≈0.22–0.43 confirms it).
> Phases A–F: dark-family abundance matrix → gradient structuring (db-RDA/PERMANOVA) → dN/dS
> selection → endemism vs global catalogues → decoding↔ecology bridge → metatranscriptomics
> validation (parallel). Target Nature Microbiology / Nat Ecol Evol; Nature-main contingent on wet lab.
- [ ] Finalize tier2 coverage (job 47540824) → merge tier1+tier2 structural-novelty catalog.
- [ ] Cross-reference novel fold families with EggNOG/Pfam-poor regions → curate the headline
      "novel fold families" table + representative structures for figures.
- [ ] Abundance-weight the dark matter (CoverM) — are novel folds rare or dominant?
- [ ] Global novelty check — are these folds truly absent from AFDB/ESMAtlas/PDB, or just
      undersampled? (already vs AFDB; consider ESMAtlas.)

## Pillar 2 — Novel / extreme metabolism  (NEW — needs envs + scripts)
- [ ] Atmospheric chemosynthesis: curated HMMs for `hhyL` (NiFe group 1h/1l) + `coxL`
      (form-I CODH) over all MAG proteomes; phylogenetic placement to confirm high-affinity
      clades. Tools: METABOLIC (`metabolism-env`) and/or hand-curated HMMs (Greening/HydDB).
- [ ] N-cycle: nitrogenase phylogeny (Nif/Vnf/Anf/Mo-independent) — test the V-nitrogenase ×
      XRF-vanadium hypothesis. amoA/comammox census.
- [ ] Desiccation/UV survival: ectoine/trehalose biosynthesis, photolyases, DNA-repair.
- [ ] Rhodopsins (proteo/xantho/heliorhodopsin) — sunlit dry-soil phototrophy.
- [ ] BGCs / natural products: antiSMASH (`workbench/scripts/ibex/81_antismash_bgc.sh`),
      then BiG-SCAPE → gene-cluster families; highlight novel NRPS/PKS.
- [ ] CAZymes (dbCAN) — carbon acquisition strategy.

## Pillar 3 — Reference bias / PICRUSt2 failure  (NEW)
- [x] PICRUSt2 already computed for the amplicon companion (v2.4.1, KO/EC/pathway + weighted NSTI,
      all 5 trips; `empty-quarter/data/processed/functional/picrust2/`).
- [x] **Compare PICRUSt2 vs MAG-encoded KO content — DONE (2026-06-26).** 436 matched samples,
      15,469 MAGs (IBEX job 47804253; `scripts/ibex/130_mag_ko_metagenome.py`). NSTI typical-soil
      (0.16) so placement is fine; PICRUSt2 recovers the housekeeping core (93% presence-recall) but
      (1) overpredicts 2.4× (median 10,543 vs 4,482 KOs/sample; `hoxH` 99%→0.2%, `anfG` 99%→0%),
      (2) abundance of the signature functions decoupled (Spearman ≈0/neg), (3) **56.7% of MAG proteins
      carry no KO** (hard KO-space ceiling). Results: `results/picrust2_bias/` (SUMMARY.md, README.md).
      Outline §C1.
  - [ ] *(optional extension)* stratify prediction error by GTDB-novelty (species vs genus-level) and
        by phylum; abundance-correlation panel for Fig. F5.

## Cross-cutting rigor  (NEW)
- [ ] Contamination/kit-ome screen vs negative controls (`SNTC`/`SPositive`) + contaminant DBs.
      *(Partly addressed: the pathogen genus→species collapse below is an independent contamination
      control — no human-pathogen bloom; 0.05% named pathogen species.)*
- [ ] Biogeography/temporal: MAG + trait distribution across 5 trips × 60 sites ×
      surface/deep/rhizosphere compartments.
- [x] **AMR / intrinsic resistome — DONE (2026-06-23).** AMRFinderPlus AMR+STRESS census on the
      15,469-MAG set (+185 isolates / 118 enrichment bins). Outline §R8; tables in
      `results/{amr_metal,pathogens,vfdb}/`; staged to IBEX scratch `…/EmptyQuarter_Data/processed/
      sequencing-results/resistome_safety/`. Three deliverables:
  - [x] Intrinsic antibiotic resistome (β-lactam ~25% / 3,816 MAGs, broad) — `results/amr_class_summary.tsv`,
        `amr_long.tsv` (Type=AMR promoted).
  - [x] Metal/metalloid (STRESS) census + XRF link — 4,159 hits / 3,187 MAGs; As 2,423, Hg 1,077;
        As-load ~ salinity r=+0.26 p=3e-5. `phase_amr_metal.py` → `results/amr_metal/`.
  - [x] Pathogen census (contamination control) — 7/15,469 named pathogen species (0.05%).
        `pathogen_census.py` → `results/pathogens/`.
  - [x] VFDB setB virulence census (dual-use) — 202 MAGs at id≥90, ~16 genuine exotoxin homologs.
        `vfdb_analyze.py` + `code/ibex/100_vfdb_census.sh` → `results/vfdb/`.
      *(Optional remaining: CARD/RGI cross-check; isolate-vs-MAG resistome comparison.)*

## Env setup still required (IBEX)
`metabolism-env` exists (has CoverM). Need to create/verify: AMRFinderPlus, antiSMASH (no env
yet despite `81_antismash_bgc.sh`), dbCAN, and curated atmospheric-gas HMM sets. Defense/CRISPR
envs are owned by the phages repo.
