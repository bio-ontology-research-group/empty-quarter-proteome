# RQ-D3 — Result

## Answer
The Empty Quarter soil community carries a **large, pervasive genomic capacity to oxidise
atmospheric H₂ and CO**. By the conservative genus-level estimate, the **median sample has ~13.5%
of its community abundance in high-affinity trace-gas oxidisers** (IQR 9.6–18.0%; full-coverage
phylum-level median ~29%, up to ~30% as a presence upper bound). This capacity is a **uniform
baseline across every compartment**, not a niche — consistent with trace-gas scavenging as a
community-wide survival strategy. **This is CAPACITY (genomic potential × abundance), NOT a measured
flux** — but it places the world's largest sand desert, a biome essentially absent from current
global trace-gas budgets, squarely on the climate-relevant side.

## Numbers

**1. Capacity index across 1,227 samples** (full community, genus rank, prevalence-weighted "any" =
H₂ or CO high-affinity guild):
- **median 13.5%** of community abundance; IQR **9.6–18.0%**; mean ~14.3%; range 0–99.9%
  (the 99.9% sample is a low-diversity outlier dominated by one guild genus).
- This index **is** the answer to "what fraction of the community carries the capacity": **~1 in 7
  cells (by 16S abundance), conservatively, belongs to a high-affinity trace-gas oxidiser.**
- Full-coverage **phylum-level median 28.9%**; presence-based **upper bound 29.9%**. So **~14–29%**
  of the community can live (in part) on air. The genus estimate is a **lower bound**: only ~65% of
  community abundance (`covered_fraction`) sits in genera with a MAG; the uncovered remainder is
  scored 0.

**2. By compartment (uniform baseline, n=1,227):**
| compartment | n | median capacity |
|---|---|---|
| Surface | 372 | 13.7% |
| Deep | 381 | 14.2% |
| PlantRoot | 474 | 12.8% |

The capacity is **flat across compartments** — a pervasive maintenance strategy, not a specialised
guild. (Of 1,227 samples, 445 mapped to a survey site via `env_matrix.tsv`; site-level breakdown in
`table_capacity_by_site.tsv`.)

**3. Genome-resolved guild backing the index** (`table_guild_genomics.tsv`):
- **3,828 MAGs** carry form-I `coxL` (atmospheric CO); high-affinity H₂ [NiFe] hydrogenases are
  dominated by the desert-clade **Group 1l (2,609 MAGs)** over the classic Group-1h "Huc" (893).
- **1,204 "core" MAGs co-encode both H₂ and CO** oxidation — built to scavenge multiple atmospheric
  reductants for survival energy.

## This is CAPACITY, not flux (explicit caveat)
The index multiplies **gene presence** (high-affinity hydB/coxL in MAGs) by **16S abundance**. It
says *who* and *how many*, **not how fast**. Realised atmospheric H₂/CO uptake — and its moisture and
temperature dependence — needs gas-chromatography soil incubations, metatranscriptomics for in-situ
expression, and kinetic characterisation of representative Group-1l enzymes
(`results/atmospheric/CLIMATE_IMPLICATIONS.md`). High-affinity oxidation also needs thin water films,
so extreme desiccation may suppress the realised flux below this genomic ceiling.

## Why it is climate-relevant
- Soil microbes are the **dominant biological sink for atmospheric H₂** (~60–80 Tg yr⁻¹, ~70–80% of
  total removal) through exactly the high-affinity [NiFe]-hydrogenases (Group 1h/1l) enumerated here
  — and H₂ is an **indirect greenhouse gas** (GWP₁₀₀ ≈ 11) via OH/methane chemistry, with a coming
  hydrogen-economy H₂ leak making the soil buffer urgent
  (`2025_soil_moisture_threshold_H2_uptake_NatComm.pdf`).
- Soils are likewise a **widespread net CO sink** (~250 Tg yr⁻¹) via form-I CODH
  (`Cordero_2019_atmospheric_CO_widespread_ISME.pdf`); trace-gas scavenging is near-universal in
  soil bacteria (>70% of taxa, Bay et al. 2021), and in desert soils it can power CO₂ fixation
  ("atmospheric chemosynthesis", Ji et al. 2017, *Nature*).
- **Drylands (~35–45% of land area, expanding) are under-represented in global H₂/CO budgets.** A
  dense, abundance-weighted high-affinity guild across the largest sand desert implies the **dryland
  contribution to these sinks is larger and more taxonomically distinct than current parameterisations
  assume** — this dataset replaces a modelled guess with a genome-resolved, abundance-weighted census.

## Files
- `table_capacity_distribution.tsv` — capacity index distribution + CO/H₂/presence/phylum variants + covered fraction.
- `table_capacity_by_compartment.tsv` — capacity by Surface/Deep/PlantRoot.
- `table_capacity_by_site.tsv` — capacity by survey site (env_matrix-mapped subset).
- `table_capacity_per_sample.tsv` — per-sample index + compartment + site.
- `table_guild_genomics.tsv` — hydB 1l/1h and coxL form-I MAG counts + 1,204 guild-core MAGs.
- `fig_tracegas_capacity.{pdf,png}` — distribution across samples + by compartment.

## References
- 2025, *Nature Communications* — soil-moisture threshold for high-affinity H₂ uptake
  (`2025_soil_moisture_threshold_H2_uptake_NatComm.pdf`).
- Cordero et al. 2019, *ISME J* — atmospheric CO oxidation is widespread in soils
  (`Cordero_2019_atmospheric_CO_widespread_ISME.pdf`).
- Bay et al. 2021, *Nat Microbiol* — trace-gas oxidation in >70% of soil taxa.
- Ji et al. 2017, *Nature* 552:400 — atmospheric chemosynthesis (H₂/CO-powered CO₂ fixation) in
  desert soils. HydDB: Søndergaard et al. 2016 / Greening et al. 2016.
