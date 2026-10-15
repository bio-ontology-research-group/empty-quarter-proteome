# §2.3 radiation/desiccation (paper scaffold A2 = RQ D4) — recompute on the 5,438 catalogue

> ## STATUS (2026-08-13): compute DONE + §2.3 DRAFTED + REGISTERED — SHIPPED 2026-08-13: Overleaf b4bdac0, GitHub main ff4b32c.
> Results §2.3 + Methods M.8 + Suppl. Tables S3/S4 written (paper/main.tex, supplementary.tex, both
> compile); registered in notes/CORRECTIONS.md (2026-08-13), PLAN.md, PAPER_STRUCTURE.md, RQ.md
> D4/D5; repo results/uv_survival/UV.md banner-marked superseded. verify.py rebinding = pending PI
> coordination. Ecotype test (step 6) NULL → species sorting, one sentence in draft.
> **Headline:** the light/desiccation dissociation of the old analysis REPRODUCES on the current
> catalogue, on a different statistical unit (communities, not genomes), and — the decisive test —
> **survives within phyla**, so it is not taxonomic composition in disguise.
> **Inputs pinned:** per-species KOs = `04_functional/repko_clean/results/rep_ko_v2.tsv` (5,438 /
> 0 unannotated); abundance = `.../coverm_current/profiles_full_recompute` (548 canonical∩profiled).
> **Runs:** all local (no SLURM); read ladder on 546 BAMs, `samtools idxstats`.
> **Open (small):** bibkeys + numeric NHEJ soil comparator (citations session, see repo PLAN.md);
> verify.py `a2.*` rebinding = PI coordination. PlantRoot framing SETTLED (root zone 5–10 cm, tracks
> light not depth). Osmolytes/salt = sibling folder `osmolyte_clean/`. Panel = 28 KOs (K01859
> excluded 2026-08-13: chalcone isomerase mislabelled as phrA; 1 carrier, also has K01669, zero
> count changes).

Old analysis: `06_metabolism/uv_survival/` (15,469 single-sample MAGs, plain `__SemiBin` naming,
job 47736827; stats in the paper repo `code/uv_survival_analysis.py`). That folder is HISTORICAL —
left untouched, read-only. Nothing here is re-annotated; only the census and the statistics are new.

## Why the design changed

The old compartment test tagged each of 15,469 MAGs with the compartment of its source sample and
compared percentages of MAGs. Three problems, all verified in this catalogue:

1. **Pseudoreplication.** MAGs from one sample are not independent observations. The old
   `p ≈ 1e-14` came from treating 15,469 genomes as 15,469 independent data points.
2. **Recovery skew.** MAG recovery per sample differs strongly by compartment — PlantRoot 27.6
   bins/sample, Deep 21.1, Surface 15.3 — so "% of recovered genomes" partly measures assembly
   success, not biology.
3. **Dereplication artefact.** The 5,438 representatives come from only 473 of 551 samples and are
   PlantRoot-skewed, so a rep-level compartment split would measure which sample won dereplication.

Replacement: the **sample** is the unit. Per-species traits × CoverM abundance → one value per
sample; replicates collapsed to one value per site×compartment; **Surface vs Deep paired within
site** (44 sites have both profiled), which holds climate and geochemistry constant.

## Verified inputs

| what | path | note |
|---|---|---|
| 29-KO damage-response panel | `06_metabolism/uv_survival/uv_markers.tsv` | reused verbatim from `scripts/ibex/124_uv_survival_census.sh` |
| per-species KOs | `04_functional/repko_clean/results/rep_ko_v2.tsv` | 5,438 rows / 0 unannotated / 690 reps had spurious KOs shed |
| species ↔ members | `05_refine/_drep/out/data_tables/{Cdb,Widb}.csv` | 12,433 HQ+MQ bins → 5,438 reps |
| GTDB taxonomy | `05_refine/_gtdbtk_refined/gtdbtk.{bac120,ar53}.summary.tsv` | 4,847 + 591 = 5,438 exactly |
| abundance | `rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute` | 548 canonical∩profiled |
| design | `metadata/sample_metadata.tsv` | 551 samples, 61 sites, compartment |

## Step 1 — census on the current catalogue (`uv_1_census.py`)

`results/group_prevalence_5438.tsv`, `marker_prevalence_5438.tsv`, `rep_traits.tsv`.
Denominator 5,438 (all annotated after `repko_v2`).

| functional group | species | % of 5,438 | old % (15,469) |
|---|---|---|---|
| base excision repair | 5,415 | **99.6** | 99.5 |
| nucleotide excision repair | 5,411 | **99.5** | 99.4 |
| DNA protection / oxidative | 5,083 | **93.5** | 91.9 |
| SOS / translesion | 4,899 | **90.1** | 94.2 |
| double-strand break / NHEJ | 4,051 | **74.5** | 78.4 |
| mismatch repair | 3,136 | **57.7** | 48.8 |
| photoprotective pigment | 2,959 | **54.4** | 58.5 |
| photoreactivation (photolyase) | 975 | **17.9** | 17.5 |

Every representative but one carries at least one panel marker (5,437/5,438).

## Step 2 — community metrics (`uv_2_community.py`)

`results/community_by_sample.tsv`, `community_summary.tsv`. 548 samples × 11 traits. Two views:
**occurrence** (share of species detected in the sample that carry the trait) and **abundance**
(share of the community by relative abundance), the latter under two denominators — whole community
including the ~51% unmapped (a hard lower bound) and the mapped community (extrapolation).

Median occurrence, by compartment:

| trait | Surface | Deep | PlantRoot |
|---|---|---|---|
| photolyase | **17.0** | **12.7** | 17.3 |
| carotenoid | **63.2** | **57.3** | 63.5 |
| NHEJ | **83.1** | **85.7** | 84.6 |
| NER | 99.8 | 99.7 | 99.8 |
| BER | 99.6 | 99.7 | 99.7 |

## Step 3 — site-paired compartment test (`uv_3_compartment.py`)

`results/compartment_tests.tsv`, `site_level_values.tsv`. Wilcoxon signed-rank on the per-site
Surface − Deep difference, 44 paired sites, Holm-corrected across traits within each metric.

| trait | Δ occurrence | sites S>D | p (Holm) | Δ abundance (mapped) | p (Holm) |
|---|---|---|---|---|---|
| photoreactivation | **+3.3** | 36/44 | 1.0e-03 | **+10.7** | 1.9e-03 |
| photoprotective pigment | **+4.8** | 37/44 | 3.6e-05 | **+10.7** | 3.6e-05 |
| dsb / NHEJ | **−2.6** |  9/44 | 3.9e-03 | **−5.6** | 2.8e-02 |
| mismatch repair | +4.3 | 32/44 | 1.1e-02 | +6.9 | 2.8e-02 |
| SOS / translesion | +0.7 | 36/44 | 1.0e-03 | +0.9 | 4.0e-02 |
| oxidative | +0.4 | 29/44 | 2.4e-02 | +0.7 | 7.3e-02 |
| NER | +0.1 | 32/44 | 4.7e-03 | +0.1 | 5.7e-03 |
| BER | −0.0 | 16/44 | 7.0e-02 | −0.1 | 1.5e-01 |

**Read this by effect size, not by p-value.** With 44 paired sites and a trait sitting at 99.6%
prevalence, a 0.1-point shift is "significant" and meaningless. The near-universal groups move by
0.1–0.9 points; the light-dependent traits move by 3–11. On the most conservative metric
(whole-community abundance, unmapped counted as absent) **only photoreactivation (+4.2, p=1.4e-04)
and photoprotective pigment (+2.9, p=6.6e-03) survive at all** — nothing else does.

## Step 3b — the composition control (`uv_4_phylum.py`)

`results/phylum_composition_test.tsv`, `within_phylum_tests.tsv`.

**Composition really does shift with compartment**, so the naive test was genuinely at risk:
Pseudomonadota **+10.0 points** at the surface (p=1.4e-03), Bacteroidota +2.6 (p=1.8e-05),
Chloroflexota −1.7 (p=9.9e-04), Thermoplasmatota −0.9 (p=1.6e-06), Actinomycetota −2.0 (p=2.9e-02).

**But the effects survive within phyla** (occurrence recomputed using only species of that phylum):

| trait | Actinomycetota | Pseudomonadota | Bacillota | Bacteroidota | Thermoplasmatota |
|---|---|---|---|---|---|
| photoreactivation | **+3.5** (1.4e-03) | **+2.6** (2.0e-02) | **+2.2** (2.1e-02) | **+1.4** (4.1e-02) | 0.0 ns |
| photoprotective pigment | **+2.9** (1.5e-03) | +0.6 ns¹ | **+5.0** (1.8e-04) | −0.6 ns | −0.9 ns |
| dsb / NHEJ | **−0.8** (2.5e-03) | −2.2 ns | **−5.9** (1.6e-02) | **−3.5** (1.6e-03) | **−5.2** (3.4e-03) |
| mismatch repair | 0.0 ns | −0.6 (5.9e-04) | +0.5 ns | +0.2 ns | +0.3 ns |

¹ already at 85% within Pseudomonadota — ceiling.

- **Photolyase rises at the surface inside four independent phyla**; carotenoid inside two.
- **NHEJ falls at the surface inside four independent phyla** — opposite direction, same design.
- **Mismatch repair fails the control**: significant in opposite directions in different phyla
  (Pseudomonadota −0.6, Chloroflexota +1.3), null elsewhere. Its whole-community effect was
  composition, not biology. **Drop it** — and note the control has teeth, since it kills a candidate.

## Step 6 — within-species ecotype test (`uv_6_ecotype.py`) — NULL, informative

For the 327 species with member genomes in both Surface and Deep (324 usable after annotation
resolution; 23/2,484 member bins unresolvable), conspecific genomes almost always share the same
trait state (photolyase state differs in only 37/324 species), and surface members do not differ
from subsurface members for any focal trait (Wilcoxon p ≥ 0.18; HQ-only variant n=77, also null,
underpowered at 2–8 informative pairs). **Meaning:** the compartment partitioning is species
SORTING, not intraspecific gene gain/loss — drafted as one sentence in §2.3.
`results/ecotype_{per_species,test}.tsv`.

## NHEJ definition breakdown (2026-08-13, for the Sharda 2020 comparator)

On 5,438 / rep_ko_v2: ku (K10979) 3,181 = 58.5% · ligD (K01971) 3,891 = 71.6% · **BOTH 3,021 =
55.6%** · either 4,051 = 74.5% (= the paper's dsb_nhej group flag). Sharda 2020 (GBE, ~6,000 mostly
cultured complete genomes) reports **22% with the complete Ku+LigD system** (Ku-only 2.3%, LigD-only
10.3% separate). **Definition-matched comparison = our 55.6% vs their 22%** — the contrast survives
alignment (2.5×), but cite 55.6, not 74.5, next to their number; note their cultured/complete-genome
denominator and the genome-size/slow-growth confound (soil Actinomycetota carry NHEJ anyway).

## Step 5 — read ladder (`uv_5_ladder.py`)

`samtools idxstats` over each sample's own-assembly BAM × `contig_to_bin_refined.tsv.gz`, splitting
reads into never-assembled / assembled-but-unbinned / binned. Median per compartment, 546 samples:

| compartment | never assembled | assembled, unbinned | binned |
|---|---|---|---|
| ALL | 44.2% | 40.1% | 11.7% |
| Surface | 50.3% | 41.1% | 7.4% |
| Deep | 48.9% | 37.4% | 11.3% |
| PlantRoot | 37.6% | 42.1% | 14.7% |

Assembly loss and binning loss are comparable in size, with assembly slightly larger. PlantRoot
recovers roughly twice as much of its community into the binned tier as Surface. This is a
within-sample ladder; CoverM's ~49% mapped is higher because the catalogue pools species across
551 samples. Surface vs Deep are comparable at every step, so the paired photolyase/NHEJ result is
not driven by differential assembly loss.

## Caveat regime

- **Capacity, not activity.** Gene presence, DNA-based abundance. No metatranscriptomics.
- **Genes are called from one representative genome per species** (median completeness 89.8%), so
  presence is under-called by roughly 10%. This biases toward the null: positive results are safe,
  a null result would be ambiguous. Member-pooled sensitivity (12,433 bins) is deferred — 66% of
  species are singletons, so pooling would favour widely-recovered species.
- **~51% of reads map to no representative.** Surface and Deep are missing equally (54.3% vs 53.7%,
  p=0.50), so the paired contrast is not biased by it; **PlantRoot maps significantly better**
  (47.3%, Kruskal p=3.9e-10), so any PlantRoot comparison must use the mapped-only denominator.
  Equal missing *magnitude* is not proof of equal missing *composition* — that residual is what an
  assembly-level gene search would close.
- **PlantRoot photolyase (17.3%) is close to Surface (17.0%) and well above Deep (12.7%).**
  Framing (sun exposure vs surface-vs-subsurface) depends on the root-sampling depth — open.

## v551 STATUS (2026-09-24)

Non-destructive rerun exists at `results/v551/` (`uv_2_community.py`, `uv_3_compartment.py`,
`uv_4_phylum.py`, `RES_OUT` env override, default unchanged). Cohort = 551 canonical∩profiled minus
1 degenerate profile (F5Dr1, 0/1,232 reads mapped) = 550 used; site-paired n_sites goes 44→45 (site
51 newly paired via the repaired 51Sr1 Surface profile). Compartment (abund_annotated_pct): photolyase
25.0/12.3 (was 24.3/12.1), carotenoid 70.7/60.9 (was 70.3/61.0), NHEJ 79.7/85.6 (was 80.0/85.7);
unmapped Surface/Deep median 54.4/53.7 p=0.47 (was 54.3/53.7 p=0.50). Within-phylum control unchanged
in direction/significance pattern. Registered values stay unchanged until promoted.

## Not done / next

1. **Draft the §2.3 repair half in `paper/main.tex`** on these numbers (line 170, currently a
   `\todo` stub); write the corresponding piece of Methods M.4/M.5.
2. **Register old→new in `notes/CORRECTIONS.md`** with basis and provenance.
3. **Osmolytes (RQ D5)** — separate panel verified against `databases/dram/kofam_ko_list.tsv`, then
   steps 2–3b with the same design.
4. **Optional armour** (not required now): member-pooled sensitivity across 12,433 bins;
   assembly-level KOfam search for the unbinned half using the existing `04_protein_catalog`.

## Files in `results/`

- `group_prevalence_5438.tsv`, `marker_prevalence_5438.tsv`, `rep_traits.tsv` — step 1
- `community_by_sample.tsv`, `community_summary.tsv` — step 2
- `compartment_tests.tsv`, `site_level_values.tsv` — step 3
- `phylum_composition_test.tsv`, `within_phylum_tests.tsv` — step 3b
- `read_ladder_by_sample.tsv`, `read_ladder_summary.tsv` — step 5
