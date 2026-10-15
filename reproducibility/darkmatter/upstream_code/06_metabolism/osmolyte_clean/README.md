# §2.3 compatible solutes / osmolytes (paper scaffold A3 = RQ D5) — first proper census, 5,438 catalogue

> ## STATUS (2026-08-13): ALL compute DONE (census, compartment, phylum control, salt+climate
> ## gradients, turnover, CLR/detection robustness) + §2.3 DRAFTED + REGISTERED — SHIPPED 2026-08-13: Overleaf b4bdac0, GitHub main ff4b32c.
> Drafted into paper/main.tex §2.3 para 3 + Methods M.8 + Suppl. Tables S3/S4; registered in
> notes/CORRECTIONS.md 2026-08-13 (incl. the phantom-OR removal), PAPER_STRUCTURE.md, RQ.md D5.
> **Headline:** compatible-solute synthesis is widespread (79.5% of species; trehalose 66%,
> ectoine 16.5%) and **flat across compartments** (small trehalose occurrence signal fails the
> phylum control). **The salt gradient goes the "wrong" way for the adaptation hypothesis:**
> trehalose-carrier abundance and osmolyte breadth DECLINE weakly with salinity (rho −0.16…−0.23).
> Resolution (step 6): salt is answered by **taxonomic turnover, not gene content** —
> **Halobacteriota tracks the salt axes strongly** (S rho +0.45 p=7e-14, Cl +0.39 p=3e-10,
> salinity_proxy +0.31 p=6e-07, n=248, lat/lon-partialled). Haloarchaea run the K+ salt-in
> strategy and don't need compatible solutes, which explains the osmolyte decline mechanistically.
> **Section story (three stressors, three response modes):** light → within-community partitioning
> (uv_survival_clean); desiccation → uniform baseline (this folder); salt → community turnover
> (this folder, step 6). Gradient caveat: the 58/224 diffuse weak gradient hits (step 5) read as
> community turnover leaking into every metric, NOT trait-specific ecology — do not cherry-pick
> individual rows from gradient_tests.tsv.
> Glucosylglycerol (3.6%) and mannosylglycerate (2.8%) are RARE, Pseudomonadota-only —
> main.tex:242 overstates them and must be rewritten.

First dedicated osmolyte analysis: the historical panel (`scripts/82_metabolism_marker_census.sh`,
15,469 MAGs, 2026-06-15, hohndor) covered only otsA/otsB + ectC/ectD as a side-table of the
metabolism census; glucosylglycerol/mannosylglycerate/betaine were never censused. The
"OR=1.10 p=0.45" honest-negative attached to A3 in main.tex:245 is a PHANTOM — its source
(`results/darkmatter_ecology/OSMOLYTE.md`) reports MWU p=0.11 on point-biserial r, no OR; and that
test belongs to the dark-family decode (§2.5 / sibling paper), not to this census.

## Panel (synthesis genes only; every KO verified vs `databases/dram/kofam_ko_list.tsv` 2026-08-13)

- trehalose: otsA K00697, otsB K01087, fused K16055; alt routes treS K13057, treY K06044 + treZ K01236
- ectoine: ectC K06720 (committed synthase); ectD K10674 = hydroxyectoine (reported separately)
- glucosylglycerol: ggpS K03692, ggpP K05978
- mannosylglycerate: mpgS K05947 / single-step mgs K18818; mpgP K07026
- EXCLUDED: transporters (opuABCD/proU/BCCT) — dual-use, weak osmotic evidence; glycine-betaine
  synthesis (choline route) not panelled this pass.

Inputs identical to `uv_survival_clean` (rep_ko_v2, GTDB refined summaries, CoverM
profiles_full_recompute 548 samples, sample_metadata). Scripts `osmo_1..4` are the uv_1..4 pipeline
with this panel; osmo_4 additionally catches the all-zero-differences Wilcoxon ValueError (latent,
never triggered, in uv_4).

## Step 1 — census (`osmo_1_census.py` → `{marker,pathway}_prevalence_5438.tsv`, `rep_traits.tsv`)

| pathway (capacity = committed synthase) | species | % of 5,438 |
|---|---|---|
| any osmolyte pathway | 4,322 | **79.5** |
| trehalose any route | 4,099 | **75.4** |
| trehalose (otsAB route) | 3,608 | 66.3 (old 15,469: otsA 67.8) |
| hydroxyectoine (ectD) | 1,028 | 18.9 (old 17.9) |
| ectoine (ectC) | 895 | 16.5 (old 17.0) |
| glucosylglycerol (ggpS) | 196 | **3.6** |
| mannosylglycerate (mpgS/mgs) | 151 | **2.8** |

Strict full-route variants: trehalose 61.8%, glucosylglycerol **0.04% (n=2)**, mannosylglycerate
0.8%. ggpP is nearly absent from the catalogue — GG "capacity" is synthase-only.

## Steps 2–3 — community + site-paired compartment test (44 sites, Wilcoxon, Holm)

Occurrence medians Surface/Deep/root-adjacent: any_osmolyte 87.1/86.4/85.2 — **flat**. By
whole-community abundance every pathway is ns except the two rare ones:

| trait | Δ occ | Δ abund (whole) | phylum control |
|---|---|---|---|
| trehalose | +2.2 (p 2.0e-03) | −1.3 ns | **FAILS** — within-phylum: Actino +0.3 ns, Pseudo −1.1 ns → composition |
| ectoine / hydroxyectoine | +1.1 / +1.2 ns | −0.1 / +0.3 ns | mixed, weak |
| glucosylglycerol | +0.8 (8.1e-04) | +0.9 (1.6e-03) | Pseudomonadota-only trait (+2.8 within, 3.6e-03) |
| mannosylglycerate | −0.5 (2.9e-02) | −0.4 (8.1e-04) | rare (0.3–0.5% within Pseudo), deep-skewed |

**Reading:** osmolyte capacity does NOT partition by 10 cm of depth. The contrast with
`uv_survival_clean` is the result: light-specific defences partition by sun exposure; break repair
partitions the other way; water-retention chemistry is a uniform baseline — consistent with matric
stress being ubiquitous in a hyperarid dune field. Glucosylglycerol/mannosylglycerate movements are
real but tiny, single-phylum, and not headline material.

## Steps 5–6 — gradients and turnover (`osmo_5_gradients.py`, `osmo_6_turnover.py`)

- **Step 5** (`results/gradient_tests.tsv`): 224 tests — every §2.3 trait (both folders, occurrence
  + whole-abundance) × salt axes (salinity_proxy, S, Cl, Na, Ca; n=248 XRF) + climate axes
  (aridity_index, precip, temp; n=453), partial Spearman on ranked lat/lon (PI's PARTIAL.md
  method), BH-FDR per block. 58 hits q<0.05 but ALL weak (|rho| 0.15–0.31), spread across
  near-universal traits, signs inconsistent → diffuse community turnover, not trait ecology.
  Trait-specific reading: **osmolyte capacity does NOT rise with salinity** (trehalose_any
  abundance rho −0.18…−0.22 vs salinity/S/Cl/Ca; osmo breadth −0.16…−0.19 vs Na/Ca). Climate:
  repair abundance rises weakly with dryness (aridity_index/precip rho −0.13…−0.18) — right
  direction, too weak to headline.
- **Step 6** (`results/salt_turnover.tsv`): phylum fraction vs salt axes, same design.
  **Halobacteriota: S +0.451 (7.4e-14) · Cl +0.385 (3.5e-10) · salinity_proxy +0.311 (5.6e-07) ·
  Ca +0.261 (3.3e-05).** Next strongest |rho| ≤0.24 (Chloroflexota/Actinomycetota negative).
  Halobacteriota is near-absent at typical sites (median 0.0%) and rises where salt does — the
  a-priori salt-in (K+) strategists. Salt is answered by turnover toward salt-in taxa, not by
  osmolyte-gene enrichment in the resident community.

## Step 7 — compositionality robustness (`osmo_7_robustness.py`) — PASSES

The Halobacteriota-salt coupling survives both compositionality checks (same partial-Spearman
design): **CLR**-transformed fractions S +0.371 (1.7e-09) / Cl +0.303 / salinity +0.168 / Ca
+0.222; **detection** (present/absent, immune to closure) S +0.441 (3.3e-13) / Cl +0.369 /
salinity +0.299 / Ca +0.244. Halobacteriota detected in **36/248** chemistry samples (rare outside
saline sites — now stated in the §2.3 draft). `results/salt_turnover_robustness.tsv`.

## Files in `results/`
`marker_prevalence_5438.tsv`, `pathway_prevalence_5438.tsv`, `rep_traits.tsv` (step 1) ·
`community_by_sample.tsv`, `community_summary.tsv` (step 2) · `compartment_tests.tsv`,
`site_level_values.tsv` (step 3) · `phylum_composition_test.tsv`, `within_phylum_tests.tsv` (3b) ·
`gradient_tests.tsv` (5) · `salt_turnover.tsv` (6) · `salt_turnover_robustness.tsv` (7)

## v551 STATUS (2026-09-24)

Non-destructive rerun exists at `results/v551/` (`osmo_2_community.py`, `osmo_3_compartment.py`,
`osmo_6_turnover.py`, `osmo_7_robustness.py`, `RES_OUT` env override, default unchanged). Cohort =
551 canonical∩profiled minus 1 degenerate profile (F5Dr1, 0/1,232 reads mapped) = 550 used; XRF salt
cohort 248→250. Halobacteriota vs S rho +0.456 p 3.0e-14 (was +0.451 p 7.4e-14); CLR/detection legs
move the same small amount. Compartment site-paired n_sites goes 44→45. Registered values stay
unchanged until promoted.

## Open / next
1. ✅ Drafted + registered + shipped 2026-08-13 (see STATUS).
2. `verify.py` `a3.osmolyte.or` phantom → real OSMOLYTE.md values (MWU p=0.11, medians +0.001
   vs −0.015, n=212) = pending PI coordination; main.tex:242 sentence already rewritten.
3. Betaine synthesis KOs if a reviewer asks; transporters stay out. Halobacteriota genus-level
   drill-down (which haloarchaea?) = optional armour, not needed for the draft.
