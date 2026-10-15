# Nature headline options for the Empty Quarter metagenomics project

**Bottom line:** no result now in hand clears the *Nature* bar without either a
decisive new measurement or a much stronger within-species comparative result.
The ranked options below are the least implausible routes, not promises.
Path shorthand: `$M=/ibex/scratch/projects/c2014/rob/metagenomics`,
`$A=/home/leechuck/Documents/papers/empty-quarter-amplicon/repository`, and
`$W=/home/leechuck/Public/software/empty-quarter-workbench`.

## 1. Relic DNA creates the apparent desert climate–diversity gradient

- **Claim:** Across the world's largest sand sea, dead/extracellular DNA forms an ecological archive that generates the apparent climate and temporal structure of microbial diversity, whereas the living community is smaller, climate-insensitive and stable.
- **What is new:** Carini *et al.* (*Nature Microbiology*, 2017) already showed that relic DNA can obscure soil diversity, while Lennon *et al.* (*mBio*, 2018) found that even large relic pools often have little effect; showing that relic DNA creates an entire biome-scale climate–diversity relationship would resolve that controversy in an extreme system and is not already done.
- **Analysis and inputs:** Reanalyse the nine pairs in `/home/leechuck/Public/software/empty-quarter/relic-dna/ASV_table_rel_abundance.tsv`, then propagate only cross-site-calibrated effects into `$A/cache/feature_table.parquet`, `$A/cache/relic_priors/relic_score_with_mag_prior.tsv`, `results/darkmatter_ecology/env_matrix.tsv` and the 551 metagenomes; paired untreated/PMA/DNase shotgun plus spike-in absolute 16S qPCR at ≥20 climate-stratified sites × all three compartments × ≥2 trips is unavoidable, although archived soils could avoid new fieldwork.
- **Verdict:** **Nature-plausible only after that measurement; specialist now.** Strongest objection: nine pairs from two sites, with no Deep samples and held-out-site AUC about 0.69, cannot support extrapolation to 1,271 samples or an absolute “60–80% relic” claim from compositional libraries.

## 2. Desert adaptation occurs through parallel within-species pangenome remodelling

- **Claim:** Across hundreds of bacterial species, the same salinity and mineral gradients repeatedly select retention or gain of trace-gas, osmolyte and repair modules within species, making accessory-gene turnover—not species replacement—the dominant genomic response to desert heterogeneity.
- **What is new:** Open pangenomes themselves are old news (Tettelin *et al.*, *PNAS*, 2005; McInerney *et al.*, *Nature Microbiology*, 2017), and Conrad *et al.* (*ISME Journal*, 2022) already showed salinity-responsive accessory genes in one *Salinibacter* population; a replicated result across 345 co-occurring species would directly address the Vos–Eyre-Walker “are pangenomes adaptive?” controversy and would be new.
- **Analysis and inputs:** Build completeness-aware metapangenomes from `$M/06_pangenome/out/{family_table.tsv,species_members.tsv}`, `$M/04_drep_self_vs_single/out/data_tables/Cdb.csv`, `RQ/C2/species_openness.tsv`, `results/darkmatter_ecology/env_matrix.tsv` and `results/atmospheric/integrate/atmo_guild_per_mag.tsv`; test environment–gene associations within species, meta-analyse direction across lineages, validate presence/absence by raw-read mapping, and require gene-tree/synteny evidence for repeated gain; no new wet lab is unavoidable.
- **Verdict:** **Nature-plausible only if the same modules move in parallel across many lineages and explain more response than taxon turnover; otherwise specialist.** Strongest objection: incomplete MAGs plus lineage–environment population structure can manufacture both gene loss and apparent adaptation.

## 3. Hot-desert microbes live on air without growing on it

- **Claim:** Atmospheric H₂ and CO sustain maintenance respiration but little primary production across the Empty Quarter, revealing a hot-desert energetic regime distinct from Antarctic atmospheric chemosynthesis.
- **What is new:** Ji *et al.* (*Nature*, 2017) already demonstrated trace-gas-powered CO₂ fixation in Antarctic desert soil, and Bay *et al.* (*Nature Microbiology*, 2021) already showed widespread active trace-gas oxidation sufficient for maintenance in diverse soils; therefore “trace gases are common” and “they support persistence” are already done—the new result must be a measured hot-versus-cold partition of energy into maintenance versus fixation.
- **Analysis and inputs:** Correct pathway absence for MAG completeness using `results/atmospheric/atmo_classification.tsv`, `results/atmospheric/integrate/{carbon_fixation_per_mag.tsv,rbcL_form_per_protein.tsv}`, `results/atmospheric/abundance/guild_abundance_by_sample.tsv` and matched cross-desert genomes; fresh-soil ambient H₂/CO uptake plus quantitative ¹³CO₂-SIP across ≥12 climate-stratified sites × three compartments at measured field water potentials, ideally with metatranscriptomics, is unavoidable.
- **Verdict:** **Nature-plausible with flux and isotope evidence; specialist genomics now.** Strongest objection: a genome census cannot show that trace gases are consumed, power maintenance, or fail to support carbon fixation in situ.

## 4. Deserts converge on an energy strategy, not on its molecular implementation

- **Claim:** Hot and cold deserts repeatedly converge on atmospheric trace-gas scavenging while using different taxa, hydrogenase classes and carbon-acquisition modes, showing that “desert adaptation” is functional convergence without a universal desert genome.
- **What is new:** Fierer *et al.* (*PNAS*, 2012) already found desert metagenomes enriched for osmoregulation and dormancy, and Bay *et al.* (2021) established broad trace-gas use; the molecular non-convergence would be new only if it extends far beyond the already-negative Group-1l aridity test.
- **Analysis and inputs:** Run raw comparator metagenomes through one workflow using `$W/data/comparator_metagenomes/comparator_samples.tsv`, `$W/docs/COMPARATOR_PIPELINE.md`, the proteomes under `$M/databases/deserts/`, `results/atmospheric/aridity/cross_biome_hydrogenase.tsv` and the EQ carbon-fixation/pangenome tables; no new wet lab is required for a genomic claim, but every accession must first be verified and a real matched mesic metagenome reference is required.
- **Verdict:** **Specialist**, with an outside chance at *Nature Ecology & Evolution*, not *Nature*. Strongest objection: study-of-origin, assembly and climate covary perfectly, so apparent molecular divergence can be a pipeline/batch artefact.

## 5. The Empty Quarter yields a chemically new natural product

- **Claim:** A conserved, desert-specific biosynthetic gene cluster from the living Empty Quarter community produces a previously unknown molecule with a new biological activity or mechanism.
- **What is new:** Computational BGC novelty is already saturated: Nayfach *et al.* (*Nature Biotechnology*, 2021) reported 87,187 putatively novel BGCs in the GEM catalogue and the 2023 SMAG soil catalogue reported 70,081; the present 98.8% complete-BGC GCF novelty is therefore not itself new.
- **Analysis and inputs:** Prioritise complete cultured-host clusters from `RQ/F2/results/{result_summary.tsv,gcf_composition.tsv,novelty_by_product.tsv}`, `$M/06_metabolism/bgc_novelty/` and `results/atmospheric/cultures/wetlab_targets.tsv`, then require cross-desert absence; long-read closure, expression or cultivation, LC–MS/NMR structure elucidation and activity/mechanism assays are unavoidable.
- **Verdict:** **Not worth it as the paper headline without a molecule; Nature-plausible only if the molecule has exceptional biology.** Strongest objection: distance from MIBiG predicts sequence novelty, not a new compound, activity or mechanism.
