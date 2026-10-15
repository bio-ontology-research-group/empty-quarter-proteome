# resistome_clean — AMRFinderPlus intrinsic resistome on the current 5,438 catalogue (§2.4 recompute)

Owner: velazqal. Created 2026-08-12 (session 01-basis-and-launch).

## STATUS (2026-08-15, session 03-bgc-resistome) — DONE + SHIPPED (GitHub dbe367c / Overleaf 2fef715)
- Array **50413467** complete: 5,438/5,438 TSVs (AMRFinderPlus 4.2.7, DB 2026-05-15.1,
  --plus), 37/37 logs clean; 2,747 species with ≥1 hit (genuine negatives elsewhere —
  every proteome scanned).
- Antibiotic census (Type==AMR, old convention: Class grouping, core+plus pooled
  4,040/142; denominator 5,438 species): **37.40% carry ≥1 determinant**;
  **BETA-LACTAM 1,335 (24.55%)** (old: 24.67% of 15,469 redundant bins — stable across
  the basis change) / AMINOGLYCOSIDE 11.93% / RIFAMYCIN 7.25% / PHENICOL 5.02% /
  MACROLIDE 4.63% / GLYCOPEPTIDE 4.38%. 23 classes.
- STRESS/metal layer RECOVERED (the historical promotion dropped it — METAL.md):
  **ARSENIC 943 (17.34%) / MERCURY 409 (7.52%)**; tabulated separately from AMR.
- Community carriage (548 samples): any_AMR median 15.0% whole / 33.5% mapped;
  BETA-LACTAM 9.4% / 21.5%; ARSENIC 7.2% / 16.3%.
- Gradients (partial Spearman controlling lat/lon, BH-FDR): arsenic occurrence ~ Ca
  ρ=+0.224 q=0.011 n=248; arsenic abundance ~ precipitation ρ=−0.243 q=2.4e-6 n=453
  (drier = more); **old "As tracks salinity proxy r=+0.26" attenuates under geographic
  control (0.24→0.12, q=0.11) — tie the claim to Ca/Na and aridity, not the proxy.**
  Mercury: no positive salt association (negative vs Cl/S).
- Pathogen census on 5,438 (amr_3, PI's panels imported verbatim): **2 named pathogen
  species (0.04%)** — 1 E. coli + 1 A. baumannii rep (old 7 MAGs collapse to 2 species
  after dereplication — RESTATES, does not reproduce, the 0.05%); pathogen-containing
  genera 46+5 reps (0.94%), environmental relatives.
- VFDB: PI's census on the CURRENT basis reused by reference (06_resistome_safety/vfdb/,
  verified by MAG naming 2026-08-15): 5,435/5,438 ≥1 hit at id≥40 → 202 at id≥90;
  surviving categories = secretion/uptake/adhesion (competition, not pathogenesis).
- Written: §2.4 paras 2–3 + Methods M.10 + fig_amr_class rebuilt (two-layer figure,
  AMR + STRESS). Registered: CORRECTIONS 2026-08-15, CATALOGUE_NUMBERS §2.4.
- SETTLED: Type==AMR convention kept for the antibiotic table (like-for-like with the
  historical β-lactam comparator); STRESS reported as its own layer, never pooled in.
- Open: none for this folder.

## Why
The historical AMRFinder run (06_mobilome/amrfinder/, β-lactam 24.67% etc.) ran on the OLD
15,469 single-binning MAG set. CORRECTIONS row A5 = OLD → recompute on 5,438. The resistome
doubles as the biosafety / contamination control in §2.4.

## Inputs
- results/rep_faa_index.tsv — 5,438 rows (mag \t faa.gz), from 06_darkmatter/drep_rep_proteomes
  (the same proteomes as the verified family basis; zcat to $TMPDIR per MAG)
- conda env `amrfinder` (rob's miniconda3); tool + DB versions → results/VERSIONS.txt (task 1)

## Pipeline
- scripts/amr_array.sh — sbatch --array=1-37 (BATCH=150, 4 CPU/8G/6h) →
  results/amrfinder/<mag>.amrfinder.tsv; idempotent
- scripts/amr_1_census.py — Type==AMR + Type==STRESS split, tax_class/amr_class naming
  (the historical amr_long.tsv had two columns named 'class') →
  results/{amr_long,amr_class_summary,stress_long,stress_class_summary,
  amr_census_summary,rep_amr_traits}.tsv
- scripts/amr_2_community.py — CoverM carriage, dual denominators →
  results/community_by_sample.tsv, community_summary.tsv
- scripts/amr_3_pathogens.py — pathogen census, panels imported from
  06_resistome_safety/code/pathogen_census.py → results/pathogen_*.tsv
- scripts/amr_4_gradients.py — partial-Spearman salt/climate blocks →
  results/gradient_tests.tsv

## Job log
- amrfinder array: **50413467** (1-37, submitted 2026-08-13) — COMPLETE.
- aggregation (amr_1..4): run directly 2026-08-15 (TSV parsing, seconds each).
