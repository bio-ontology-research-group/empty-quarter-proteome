# Carbon-fixation mechanism behind atmospheric scavenging — KNOWN (Calvin cycle), not novel

Generated 2026-06-16 (`scripts/ibex/82j_carbon_fixation_census.sh` gene-level + `82i_atmo_autotrophy_fba.sh`
FBA). The trace-gas *enzymes* (coxL/hydB) supply energy/electrons; this asks what **carbon-fixation pathway**
they power — the actual "atmospheric chemosynthesis" mechanism.

## Answer (gene/KO-based, model-independent — the authoritative result)
Among the 6,010 guild MAGs:
| pathway | diagnostic genes | guild MAGs |
|---|---|---|
| **Calvin–Benson–Bassham (CBB)**, RuBisCO-form-validated | carboxylating `rbcL` (form I/II/III) **+** `prk` (K00855) | **974 (16.2%)** |
| &nbsp;&nbsp;(naive, before form validation) | `rbcL` any (K01601) + `prk` | 1,141 (19.0%) |
| reductive TCA | ATP-citrate lyase `aclA/B` (K15230/1) | **7** |
| Wood–Ljungdahl | `acsB/cooS` (K00198) | **1** |
| no complete autotrophic pathway | — | ~4,900 (82%) |

(`rbcL` alone 1,875; `prk` alone 1,972.) **The mechanism is the canonical Calvin cycle** — rTCA and
Wood–Ljungdahl are essentially absent. **There is no evidence for a novel carbon-fixation pathway.**

### RuBisCO-form validation (`82k`/`82l`) — excludes non-fixing form IV RLPs
`rbcL` (KEGG K01601 / EggNOG) lumps carboxylating RuBisCO with **form IV "RuBisCO-like proteins" (RLPs)** that
do **not** fix CO₂ (Tabita et al. 2008, *J Exp Bot*). Classifying the 2,292 guild `rbcL` proteins by DIAMOND
best-hit against a form I/II/III/IV reference panel:
- **form I 1,172 · form IV (RLP) 1,062 · form III 53 · form II 0 · no-hit 5** — i.e. **~46% of the raw `rbcL`
  hits are non-carboxylating RLPs**, which inflated the naive count.
- **1,122 guild MAGs carry a bona-fide carboxylating `rbcL`**; of those **974 also encode `prk`** = the
  form-validated CBB chemoautotroph set (**16.2%** of the guild). Where carbon is fixed it is **form-I CBB**.
- Method discriminated by similarity, not reference abundance (0 form II assigned despite 25 form II refs).
- Caveat: the form-I reference panel was thin (2 seqs; the Actinobacteria UniProt query returned 0 reviewed),
  and assignment is best-hit (not residue-level). Residue-level confirmation of the 19 catalytic positions, or
  a form-I phylogeny analogous to CoxL, would further tighten the 974. Files: `rubisco_summary.txt`,
  `rbcL_form_per_protein.tsv`, `rubisco_refs.faa`.

### CBB by phylum (the candidate atmospheric chemoautotrophs)
Actinomycetota **916** (25.6% of guild Actinomycetota) ≫ Pseudomonadota 142 (13.9%) > Chloroflexota 76 (18.4%);
negligible elsewhere. CBB-powered atmospheric chemosynthesis is an **Actinomycetota** story, the desert-dominant phylum.

## Two-tier interpretation
1. **~19% (1,141 MAGs) are candidate CBB chemolithoautotrophs** — couple H₂/CO oxidation to **Calvin-cycle**
   CO₂ fixation = atmospheric chemosynthesis (Ji et al. 2017, *Nature*, shown in Antarctic deserts; here in the
   largest **hot sand** desert, via the *known* mechanism).
2. **~81% lack a complete fixation pathway** — for the majority, trace-gas oxidation is **supplementary
   maintenance/survival energy** (mixotrophy), not net carbon fixation — the Greening-lab survival-energy paradigm.

## Why the gapsmith FBA could not be used to call the mechanism (honest caveat)
`82i` ran the FBA with **scavenging as the objective and biomass ≥ 0** (per Robert — gapsmith gap-filled with a
biomass objective, so plain biomass-FBA bypasses the gas reactions). Results: **0% grow autotrophically** on a
trace-gas-only medium, and the apparent "rTCA fixation flux" was driven by **2-oxoglutarate synthase** —
a non-diagnostic central-metabolism enzyme, **not** real carbon fixation. The models showed **0 RuBisCO** because
the gapseq/ModelSEED reconstruction **omits the RuBisCO reaction even when `rbcL` is in the genome** (a known
reconstruction gap), so the models cannot represent CBB autotrophy. → **The FBA confirms these are not obligate
autotrophs, but the carbon-fixation mechanism must be read from genes, not the models.** (`autotrophy_fba_*`.)

## Caveats / refinement
- KO-based: `rbcL`+`prk` co-occurrence strongly implies functional CBB, but **RuBisCO form** (carboxylating
  form I/II vs non-carboxylating RuBisCO-like form IV) is not resolved here — a phylogenetic check (as done for
  coxL form I) on the 1,875 `rbcL` would tighten the 1,141 CBB count. Flagged.
- Capacity, not in-situ flux (resolve with metatranscriptomics / ¹³CO₂ labelling on the cultured CBB+guild strains).

## Files
`carbon_fixation_per_mag.tsv` (per guild MAG: rbcL, prk, CBB_complete, pathway, phylum, gas markers),
`carbon_fixation_by_phylum.tsv`, `carbon_fixation_summary.txt`; FBA: `autotrophy_fba_per_mag.tsv`, `autotrophy_fba_summary.txt`.
