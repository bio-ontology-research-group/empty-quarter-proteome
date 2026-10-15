# RQ-D3 — Methods

**Question.** What is the community's atmospheric trace-gas (H₂/CO) sink **CAPACITY**, and is it
climate-relevant?

**Capacity, not flux.** Every number here is **genomic potential × community abundance** — the
fraction of the desert microbial community that *can* oxidise atmospheric H₂/CO — **NOT a measured
gas-uptake rate**. Realised flux requires gas-chromatography incubations + metatranscriptomics
(`results/atmospheric/CLIMATE_IMPLICATIONS.md`).

**Inputs.**
- `results/atmospheric/abundance/guild_abundance_by_sample.tsv` — per-sample, prevalence-weighted
  community guild fraction. Columns: `table` (full | alive-non-relic), `rank` (genus | phylum),
  `sample`, `compartment`, `metric` (any | CO | H2 | any_presence), `guild_fraction`,
  `covered_fraction`. 1,227 samples (16S amplicon companion).
- `results/atmospheric/abundance/guild_abundance_summary.tsv` — means by table × rank × metric × scope.
- `results/atmospheric/{guild_core_mags.txt (1,204 MAGs), hydB_class_summary.tsv, coxL_form_summary.tsv}`
  — genome-resolved guild: high-affinity [NiFe] hydrogenase classes and form-I CO dehydrogenase.
- `results/darkmatter_ecology/env_matrix.tsv` — sample → site / compartment (547 samples; id mapped
  by stripping the run prefix, `<run>_<envsample>` → `<envsample>`).

**Guild definition (upstream, in workbench).** A MAG joins the high-affinity atmospheric trace-gas
guild if it encodes a high-affinity [NiFe] hydrogenase (Group 1h "Huc" or the desert-dominant Group
1l) for H₂ and/or a form-I `coxL` CO dehydrogenase for CO — markers classified by HydDB-style
phylogenetics (`hydB_class.tsv`, `coxL_class.tsv`). Prevalence per taxon = guild MAGs / total MAGs
(denominator = all 15,469 GTDB-classified MAGs).

**Capacity index (this analysis, `run.py`, stdlib + matplotlib, run from repo root).**
Per sample, `capacity_index = Σ_taxon relative_abundance(taxon) × guild_prevalence(taxon)`, i.e. the
community-abundance fraction carried by trace-gas-capable taxa. We take the **conservative, primary**
definition = **full community, genus rank, metric "any"** (H₂ *or* CO high-affinity guild); genus-level
linkage is precise but only `covered_fraction` (~65%) of community abundance sits in genera with ≥1
MAG (the rest counted as 0 → lower bound). We also report the **phylum-level full-coverage** estimate
(~98% covered, treats within-phylum prevalence as uniform → upper-ish) and the presence-based upper
bound (`any_presence`). Distribution summarised as median / IQR / range across the 1,227 samples and
broken out by compartment (Surface / Deep / PlantRoot) and, where mappable, by site.

**Figure.** `fig_tracegas_capacity.{pdf,png}` — (A) histogram of the capacity index over 1,227
samples with median; (B) boxplot by compartment.

**Caveats.** Capacity ≠ flux. Prevalence-weighting assumes MAG-estimated within-genus guild
prevalence transfers to amplicon members of that genus. Genus-level underestimates (community in
MAG-less genera scored 0); phylum-level overestimates (within-phylum uniformity). The truth lies
between; both are reported. The single 99.9% sample is a low-diversity outlier dominated by a guild genus.

**Climate context / references.**
- Soil microbes are the dominant biological sink for atmospheric H₂ (~60–80 Tg yr⁻¹, ~70–80% of
  removal) via high-affinity [NiFe]-hydrogenases:
  `resources/literature/trace-gas-deserts/2025_soil_moisture_threshold_H2_uptake_NatComm.pdf`.
- Soils are a widespread net sink for atmospheric CO (~250 Tg yr⁻¹) via form-I CODH:
  `resources/literature/trace-gas-deserts/Cordero_2019_atmospheric_CO_widespread_ISME.pdf`.
- Trace-gas scavenging is a near-universal soil survival trait — >70% of soil taxa (Bay et al. 2021,
  *Nat Microbiol*); atmospheric chemosynthesis powering carbon fixation in desert soils (Ji et al.
  2017, *Nature* 552:400). HydDB classification: Søndergaard et al. 2016 / Greening et al. 2016.
