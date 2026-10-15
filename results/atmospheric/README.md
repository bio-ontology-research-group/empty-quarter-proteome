# Atmospheric trace-gas metabolism — phylogenetically confirmed ("life on air")

> **SUPERSEDED (2026-08-11) — HISTORICAL, 15,469-MAG basis.** These are the June 2026 outputs on the old
> single-sample catalogue (MAG ids `<sample>__SemiBin_N`). The paper's §2.2 now uses the **5,438** catalogue:
> 946 high-affinity hydB / 1,147 form-I coxL / 1,822 guild / 271 both; abundance 18.0% whole-community
> (42.9% of mapped). Current tables: `../atmospheric_5438/`; provenance:
> `M/06_atmospheric/lifeonair_clean/README.md`. Kept for audit/provenance only — do not cite as a paper number.


Generated 2026-06-16 from `06_atmospheric/` on IBEX (jobs 47554112 / 47554035 / 47554113).
Confirms, by clade/form assignment, that a large fraction of the 15,469 EQ soil MAGs encode
**high-affinity, atmospheric** trace-gas-oxidising enzymes — not merely a coxL/hydB KO. This is the
methodological upgrade from the first-pass marker census (`../metabolism_*`).

## Methods (field-standard, citable)
- **H₂ — [NiFe]-hydrogenase class**: GreeningLab **HydDB** (Søndergaard et al. 2016, *Sci Rep* 6:34212).
  DIAMOND `--ultra-sensitive --max-target-seqs 1` vs `HydDB_all_hydrogenases.faa` (3,260 group-labeled
  refs); class = best-hit group. **Groups 1h + 1l = high-affinity, O₂-tolerant atmospheric H₂ uptake.**
- **CO — Mo-CODH form I**: Cordero et al. 2019, *ISME J* (PMC6794299). DIAMOND vs three RefSeq form-I
  references (*Oligotropha carboxidovorans* WP_013913730.1, *Mycobacterium smegmatis* WP_003892166.1,
  *Natronorubrum bangense* WP_006067999.1); **form I = >55% identity & ≥90% query coverage**. Independent
  2nd axis: active-site motif (form I **AYRGAGR**). A protein is called form I on concordance of the two.

## Files
| File | What |
|---|---|
| `atmo_classification.tsv` | per-protein call: mag_id, protein_id, marker, assignment, pident, qcov, **is_high_affinity** |
| `hydB_class_summary.tsv` | HydDB class → n_proteins, n_MAGs |
| `coxL_form_summary.tsv` | CoxL form → n_proteins, n_MAGs |
| `guild_core_mags.txt` | the 1,204 MAGs encoding **both** confirmed atmospheric H₂ and CO oxidation |
| `coxL.nwk`, `hydB.nwk` | FastTree of dereplicated queries + references (figure input) |

## Headline numbers (of 15,469 MAGs)
- **Atmospheric H₂ oxidation**: 6,602 hydB → **6,278 high-affinity** (95%): **Group 1l 5,292 / Group 1h 986**;
  minor tail across 10 other classes (2a 246, 1d 32, 1c 14 …) → not a one-reference artifact; median id 58%.
  **3,386 MAGs (22%)** carry a high-affinity atmospheric H₂-uptake hydrogenase.
- **CO oxidation**: 13,675 coxL → **5,645 form I** (identity + AYRGAGR motif concordant), 8,003 ambiguous,
  27 form II. **3,828 MAGs (25%)** carry a bona fide form-I CO dehydrogenase.
- **Guild core**: **1,204 MAGs (~8%) encode BOTH** confirmed atmospheric H₂ *and* CO oxidation.

The **Group-1l dominance** is the recently-described arid/desert-soil high-affinity signal (vs the
classic Group-1h "Huc" of mesic soils) — observing it community-wide in the most carbon-starved biome on
Earth is the headline: the Empty Quarter microbiome is, to first order, **powered by atmospheric trace gases**.

## Status / caveats
- Sequence-based clade assignment confirms **genetic capacity**, not in-situ rate. Strengthened (not proven)
  by the Stage-D layers: gapseq/gapsmith FBA feasibility, operon completeness, and 16S abundance weighting.
- The 8,003 "ambiguous" coxL are CoxL-family but below the form-I threshold (likely form II / coxL-like);
  reported separately, not counted as CO oxidisers.
- Fragmented ORFs flagged by length in `atmo_classification.tsv` (pident/qcov columns).

See `CLIMATE_IMPLICATIONS.md` for why this matters beyond the desert.

## Stage D — integration (`integrate/`)
- **`atmo_guild_per_mag.tsv`** — per-MAG guild call + GTDB phylum/class + trip/compartment + gapsmith model
  corroboration (`model_has_H2ase`, `model_has_CODH`).
- **`atmo_guild_by_phylum.tsv`** — guild by phylum: **Actinomycetota** dominant (2,175 CO / 2,198 H₂ / **799 both**),
  then **Chloroflexota** (332 / 305 / **224** — highest core fraction), **Gemmatimonadota** (132 / 257 / 114),
  Acidobacteriota (50 / 118 / 35). **Pseudomonadota** are CO-rich but H₂-poor (1,007 / 24) — i.e. the
  atmospheric-H₂ guild is specifically the actinobacteria/chloroflexi/acidobacteria, exactly the clades
  Greening-lab work flags as the dominant soil H₂ scavengers.
- **gapsmith FBA cross-check**: of the 648 guild MAGs that have a gap-filled metabolic model, **616 (95%)
  independently contain a hydrogenase reaction** and 87 a CO-dehydrogenase reaction — metabolic-reconstruction
  corroboration of the sequence call. (Only 648/6,010 guild MAGs currently have a gapsmith model; corroboration
  is on that subset.)
- **`atmo_operon_context.tsv`** — operon neighbourhood (±5 genes) of every confirmed coxL/hydB. 9,007 of
  11,923 (76%) sit next to a Bakta **"hypothetical protein"** — but ⚠️ that is an **annotation-gap proxy, not
  structural novelty**. The proper crosswalk to the dark-matter catalog (`OPERON_DARKMATTER.md`) shows the
  real figure is **1.5%**: gas operons are NOT enriched for structurally-novel proteins. Use 1.5%, not 76%.
- **`OPERON_DARKMATTER.md` + `operon_darkmatter_context.tsv` + `dark_gas_families.tsv`** — the corrected,
  dark-matter-catalog-based analysis (job 82e). **Clean negative: 177/11,923 (1.5%) of gas operons have a
  genuinely structurally-novel neighbour; no dark orthogroup recurs near gas genes in >7 MAGs.** The
  novel-fold-in-gas-operon hypothesis does not hold — gas machinery is characterizable; the paper's strength
  is the guild itself, not a novel-protein-in-operon story.

## Cultures / isolates — experimental targets (`cultures/`)
Same classifier on the 1,221 EQ culture genomes (1,036 enrichment + 185 isolates; all EggNOG-annotated).
- **CO (form-I CODH)**: 144 proteins → **86 genomes** (84 enrichment + 2 isolates).
- **High-affinity H₂**: only **16 Group-1h** proteins, **0 Group-1l** — the environmentally dominant 1l clade
  is **absent from current cultures** (a clean "dominant atmospheric scavenger still uncultured" result;
  motivates targeted H₂-amended oligotrophic isolation).
- **`wetlab_targets.tsv`**: strains carrying high-affinity enzymes (gene IDs included). **13 enrichment strains
  encode BOTH** CO + H₂ — prime targets for gas-uptake kinetics + metatranscriptomics:
  `C-16.1, C-24.2, D16_ont.8, D16_pb.5, D23_ont.3, D23_pb.9, MO-1.21, MO-2.9, MO-3.2, MO-5.3, MO-P-5.1, NaB-1.4, NaB-6.6`.

## 16S abundance weighting (`abundance/`)
Mapping the guild onto the 1,227-sample amplicon community (GTDB-linked at genus + phylum):
**~14% (genus, conservative) – 29% (phylum) of the EQ soil community belongs to taxa with the genomic capacity
for atmospheric trace-gas oxidation** (up to ~74% presence-based upper bound); CO capacity (13–20%) more
widespread than H₂ (2–16%); **uniform across Surface/Deep/PlantRoot** = a pervasive baseline, not a niche.
Secondary: guild fraction is lower in the 'alive' (non-relic) subset — possible dormancy/persistence signal,
flagged cautiously. See `abundance/README.md`. Script: `../../code/make_16s_weighting.py`.

## Prior literature
A curated survey of known trace-gas scavenging in hyper-arid deserts (with PDFs) is in
`../../resources/literature/trace-gas-deserts/` — see its `README.md`.
