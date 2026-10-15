# Reading notes — mined facts per paper

Distilled from reading the actual PDFs (via PyMuPDF). Each entry: exact numbers + how it bears on our
manuscript. Read status tracked so we don't re-read. Page-level quotes are paraphrased unless in "quotes".

---

## AFESM 2025 — Yeo, Han, Bordin, Lau, ... Jones, Orengo, Steinegger. bioRxiv 2025.04.23.650224
`structure-darkmatter/Yeo_2025_AFESM_structure_universe_bioRxiv.pdf` · 44 pp · **read: abstract+intro (2026-08-09)**

**What it is.** Combines AFDB (>200M / 214M AF2 models) + ESM Metagenomic Atlas (>600M ESMFold models from
MGnify metagenomes) into **AFESM = 821M entries**, clustered by sequence+structure into **5.12M non-singleton
structural clusters**.

**The load-bearing result for our §3 (the consensus we must not contradict):**
- "Initial ESMfold-based predictions revealed **no novel domain folds**; re-predicting **2.3 million** proteins
  with AlphaFold2 yielded **only one new fold**, suggesting **near-saturation** of the domain space and
  limitations of predictors."
- Novelty that *does* survive is in **combinations**, not new folds: **11,941 multi-domain architectures not
  previously observed**.
- TED context they cite: ~365M domains; **77%** of non-redundant domains match CATH; **>7,000 putative novel
  folds** across the tree of life (the TED/Lau benchmark).

**How we use it.**
- Confirms our framing is correct: **§3 must be calibration-first ("bounded and calibrated"), NOT a novel-fold
  count.** A "we found N new folds" headline is refuted on sight by AFESM+TED.
- Our differentiator is legitimate and non-overlapping: **single-biome attribution**, both TM-score
  normalisations, and a **blinded known-fold error floor** — none of which AFESM does (they work globally).
- Position residual candidates as "a defined minority surviving explicit controls," explicitly consistent with
  near-saturation, not against it.
- Cite as `yeo2025afesm` (already in references.bib).

---

## Qaidam Basin 2025 — Liu, Wang, Zhang, Lin. *Scientific Data* 12:1795 · doi:10.1038/s41597-025-06085-3
`desert-mags-comparators/Qaidam_2025_1773MAGs_Marsanalog_SciData.pdf` · 7 pp · **read: abstract (2026-08-09)**
- 58 soil samples, Mars-analog Qaidam Basin (hyperarid <45 mm/yr, ~2800 m, high UV, saline). **1,773 mMAGs +
  2,060 vMAGs; >94% novel taxa.** HQ = 327 mMAGs (≥90%/≤5%), 325 vMAGs. Top phyla: Actinomycetota (565),
  Halobacteriota (111, archaea), Uroviricota (836, viral).
- **Use (§1 comparator):** the mirror image of our result. They report >94% *taxonomic* novelty from an
  under-sampled desert; we report taxonomic *conventionality* (0 novel phyla; Actinomycetota/Pseudomonadota).
  Both deserts are Actinomycetota-dominated. Sharpens "known organisms, unknown biology": our novelty is
  functional/structural, not phylogenetic. (Caveat we can note: their ">94% novel" is genome-level from a
  reference-blank region — partly a database-gap signal; ours is HQ-adjudicated conventional placement.)
  New bib key: `liu2025qaidam`.

## GMSC 2024 — Duan … Bork, Coelho. *Nat Commun* · doi:10.1038/s41467-024-51894-6  (bib: `duan2024gmsc`)
`smorfs-smallproteins/GMSC_2024_global_smallproteins_catalog_NatComm.pdf` · 11 pp · **read: abstract (2026-08-09)**
- Global smORF catalog from 63,410 metagenomes (75 habitats) + 87,920 isolate genomes → **965 million
  non-redundant smORFs**. smORFs (<100 codons) ≈ 5–10% of annotated genes; archaea proportionally richer.
- **Use (§2):** anchors that the short-protein space is enormous and underexplored — consistent with our finding
  that **length (not completeness) discriminates the dark proteome** (dark skews <100 aa). Also flags Prodigal
  under-calling of smORFs (methods caveat).

## Low-biomass contamination guidelines 2025 — *Nat Microbiol* Consensus Statement 10:1570–1580 · doi:10.1038/s41564-025-02035-2
`contamination-armor/LowBiomass_2025_contamination_guidelines_NatMicro.pdf` · 11 pp · **read: abstract (2026-08-09)**
- Consensus statement; **explicitly names "hyper-arid soils"** among low-biomass environments where contamination
  becomes critical near the detection limit. Provides minimal reporting standards + removal workflows for every
  study stage.
- **Use (Methods/armor):** the citation that justifies our contamination controls (extraction blanks, decontam,
  index-hopping estimate) and pre-empts the kit-ome/relic-DNA desk-reject. New bib key: `contamination2025guidelines`.

## Delgado-Baquerizo 2019 — "Increases in aridity … drastic shifts in dryland microbial network assembly"
`dryland-ecology/DelgadoBaquerizo_2019_aridity_shifts_assembly.pdf` · 10 pp · **read: abstract (2026-08-09)**
- Aridity increases drive drastic reorganisation of dryland microbial co-occurrence networks (global drylands).
- **Use (§4/ecology):** supports aridity as the master variable structuring the community, backing the
  environmental-structuring of the dark families. New bib key: `delgadobaquerizo2019networks` (confirm DOI).

## Zhang et al. 2022 — environmental selection overturns distance-decay. *eLife* 11:e70164 · doi:10.7554/eLife.70164
`biogeography-evolution/EnvSelection_2021_distancedecay_soil_eLife.pdf` · 19 pp · **read: abstract (2026-08-09)**
- Soil prokaryotes across Chinese grasslands to 4000 km: a **U-shaped** distance-decay — similarity *increases*
  past tipping points (~1760–1920 km) because environmental heterogeneity decreases across biotas. Environmental
  selection, not dispersal limitation, drives large-scale pattern.
- **Use (§4 framing discipline):** supports stating our result as **environmental selection / database-relative
  specificity under an explicit search rule, NOT geographic endemism** (matches CLAUDE.md caution). New bib key:
  `zhang2022envselection`.

## Waschulin et al. 2021 — long-read BGCs, Antarctic soil. *ISME J* 16:101–111 · doi:10.1038/s41396-021-01052-3
`metabolism-bgc/Antarctic_2021_longread_novel_BGCs_ISME.pdf` · 11 pp · **read: abstract (2026-08-09)**
- Long-read metagenomics (Mars Oasis, Antarctica) recovers **>1,400 mostly full-length BGCs** from uncultivated
  lineages; highly divergent BGCs across Acidobacteriota/Verrucomicrobiota/Gemmatimonadota + a candidate novel
  RiPP family.
- **Use (BGC subsection + Methods caveat):** short reads fragment BGCs (our antiSMASH calls need contig-edge
  flags; no novelty fraction without long-read validation). Comparator for cold-desert biosynthetic dark matter;
  BGC composition tracks lineage poorly (supports decoupling). New bib key: `waschulin2021bgc`.

---

# Batch 2 — ingested 2026-08-12 (17 new PDFs + 5 verified from the June trace-gas folder)

All metadata below verified this session against the PDF page 1, Crossref (`api.crossref.org/works/<doi>`)
or the Europe PMC REST API. Bib keys marked **added** are now in `paper/references.bib` and
`manuscript/references.bib` (the two files are byte-identical).

## §1.2–1.3 — dark matter, domains, structural novelty

### Lau et al. 2024 — The Encyclopedia of Domains (TED). *Science* 386:eadq4946 · doi:10.1126/science.adq4946
`structure-darkmatter/Lau_2024_TED_encyclopedia_domains_Science.pdf` · 15 pp · **read: abstract + results + methods (2026-08-12)** · bib: `lau2024ted` (**added**)
- Deep-learning domain parsing of all 214M AFDB models: **~365 million domains** across **>1 million taxa**,
  **>100 million** of them undetectable by sequence methods. **77% of non-redundant domains** match an existing
  CATH superfamily (p2, p4). **>10,000 new superfamily–superfamily structural interactions.**
- **The methods precedent that matters for §1.3.** Their novel-fold pipeline is a filter cascade with published
  thresholds, and the attrition is the point (p11): 41,879,858 non-CATH-labelled clusters →
  **13,820,550** pass globularity (normalised radius of gyration <0.356, packing density >10.333) + ≥6 secondary
  structure elements (STRIDE) → **8,612,318** pass pLDDT80 ≥90 → Foldseek TM-align exclusion at **>0.56** with
  60% coverage against PDB/CATH/ECOD/SCOPe → **240,674** candidate clusters → final **7,427 putative novel folds**
  (down from 24,653 before chopping checks).
- **Use (§1.3):** this is the strongest available justification for our own gate. It licenses (a) reporting a
  *filtered minority* rather than the raw no-hit count, (b) an explicit confidence gate (they use pLDDT80 ≥90;
  we use mean pLDDT ≥0.70 — ours is looser, worth stating), and (c) a TM threshold above 0.5 to absorb Foldseek's
  approximation of true TM-align. It also sets the denominator: 7,427 candidate novel folds for the *entire*
  AFDB. Any single-biome number of comparable magnitude would be implausible and should be presented as such.
- Cross-check with AFESM (Yeo 2025, above): AFESM cites TED's ">7,000 putative novel folds" and reports that
  re-predicting 2.3M proteins with AF2 yielded only one genuinely new fold. The two together define the
  near-saturation consensus our §1.3 must sit inside.

### Waman et al. 2025 — CATH v4.4. *Nucleic Acids Res* 53:D348–D355 · doi:10.1093/nar/gkae1087
`structure-darkmatter/Waman_2025_CATHv44_expansion_NAR.pdf` · 8 pp · **read: abstract (2026-08-12)** · bib: `waman2025cath` (**added**)
- CATH v4.4 absorbs ~90 million predicted TED domains plus ~64,844 new experimental PDB domains. Superfamilies
  **5,841 → 6,573**, folds **1,349 → 2,078**, architectures **41 → 77**. FunFam sequence coverage +276%.
- **Their own caveat, quoted:** "TED data comprises predicted structures, so these new folds and architectures
  remain hypothetical until experimentally confirmed."
- **Use (§1.3, one sentence):** the reference-database version statement and the precedent for calling
  predicted-structure novelty *hypothetical*. Our "reference-relative candidate folds, not experimentally
  validated structures" line is exactly the language CATH itself uses — worth citing so it reads as field
  convention rather than hedging.

### Harteveld et al. 2024 — Genesis VAE, dark-matter folds. *Cell Systems* 15:898–910.e5 · doi:10.1016/j.cels.2024.09.006
`structure-darkmatter/Harteveld_2024_darkmatter_folds_deeplearning_CellSystems.pdf` · 19 pp · **read: abstract (2026-08-12)** · bib: `harteveld2024darkfolds` (**added**)
- A *design* paper: a convolutional VAE (Genesis) coupled to trRosetta generates de novo proteins with
  specified architectures, including "dark-matter" folds not seen in nature, validated by high-throughput
  protease-resistance assays.
- **Use (cite-only, §1.3 discussion at most).** Note the term collision: their "dark matter" means
  *designable-but-unobserved* fold space, ours means *unannotated natural sequence*. Do not cite it as evidence
  about natural novelty. Its only real use is a half-sentence that unoccupied fold space exists and is
  physically buildable, which is orthogonal to whether desert soil occupies any of it. Low priority.

### Rodríguez del Río et al. 2024 — FESNov. *Nature* 626:377–384 · doi:10.1038/s41586-023-06955-z
`structure-darkmatter/RodriguezdelRio_2024_FESNov_unknown_genes_Nature.pdf` · **read: abstract (2026-08-12)** · bib: `rodriguezdelrio2024fesnov` (already present)
- **149,842 environmental genomes** → **404,085** functionally and evolutionarily significant novel (FESNov)
  gene families exclusive to uncultivated taxa, **nearly tripling** the number of described prokaryotic gene
  families. All families span multiple species and show **strong signals of purifying selection**.
- **1,034 novel families distinguish entire uncultivated phyla, classes and orders** (candidate synapomorphies).
  Functional associations predicted for **32.4%** of families, of which **4,349** are high-confidence.
- **Use (§1.2, load-bearing):** FESNov is the closest prior art to our dark-family ladder and it sets the
  criterion we already meet — *multi-species occurrence + purifying selection* is what upgrades a dark cluster
  from "possible spurious ORF" to "gene family". Our 23,382 recurrent dark families with filtered pairwise
  dN/dS 0.391 (0.499 at the ≥30-member floor) is the same argument by a different estimator, so §1.2 can say
  our families qualify under the FESNov standard rather than inventing a bar.
- **Also the honest comparator for §3.1:** FESNov's novelty is *taxon*-exclusive (uncultivated clades); ours is
  *biome*-attributed under an explicit search rule. Different axes — say so, or a reviewer will read our
  21.4% as competing with their 404,085.

## §2.1 — reference bias (16S → function)

### Sun, Jones & Fodor 2020. *Microbiome* 8:46 · doi:10.1186/s40168-020-00815-y
`reference-bias/Sun_2020_PICRUSt_accuracy_sampletypes_Microbiome.pdf` · **read: abstract + results (2026-08-12)** · bib: `sun2020inference` (already present)
- 7 matched 16S/shotgun datasets (human, non-human animal, **soil**). Predicted-vs-shotgun gene composition
  correlates strongly, **Spearman ρ 0.53–0.87** — but the correlation **survives permuting gene abundances
  across samples**: one soil sample gives ρ=0.85 unpermuted and **ρ=0.84 permuted** (p2, Fig. 1b,c).
- **For the two soil datasets, unpermuted correlations were not significantly different from permuted ones for
  all three tools** (PICRUSt, PICRUSt2, Tax4Fun) (p2). The correlation is carried by the average gene-abundance
  profile, not by any sample-specific signal.
- **Use (§2.1, load-bearing):** this is the citation that makes our replacement of PICRUSt2 a methodological
  necessity rather than a preference. The right framing: 16S-inferred function reproduces *what is typical*,
  and soil is the sample type where that failure is most complete. Pairs with our 56.7%-no-KO result — the
  inference cannot represent proteins that are absent from the reference genomes it interpolates over.

### Matchado et al. 2024. *Microbial Genomics* 10:001203 · doi:10.1099/mgen.0.001203
`reference-bias/Matchado_2024_16S_prediction_limits_MicrobGenom.pdf` · **read: abstract + results (2026-08-12)** · bib: `matchado2024limits` (already present)
- Systematic evaluation of PICRUSt2, Tax4Fun2, PanFP and MetGEM against matched shotgun and simulated data,
  plus a custom rrnDB 16S copy-number normalisation. Conclusion: these tools "generally do not have the
  necessary sensitivity to delineate health-related functional changes in the microbiome and should thus be
  used with care."
- **Confirms Sun 2020 explicitly** (p7): Spearman correlations are not affected by permutation; they adopt
  Sun's recommendation that **differential analysis, not correlation, is the correct performance measure** (p5).
  They also note prior reviews flag the **lack of reference genomes specifically for soil ecosystems** (p3).
- **Use (§2.1):** the independent replication. Cite Sun 2020 for the permutation result and Matchado 2024 for
  "and this was confirmed four years later with more tools"; the soil-reference-gap sentence is the bridge to
  our catalogue's justification.

### Wright & Langille 2025 — PICRUSt2-SC. *Bioinformatics* 41:btaf269 · doi:10.1093/bioinformatics/btaf269
`reference-bias/Wright_2025_PICRUSt2SC_reference_update_Bioinformatics.pdf` · 5 pp · **read: abstract + intro (2026-08-12)** · bib: `wright2025picrust2sc` (**added**)
- Updated PICRUSt2 reference database: bacteria **19,493 → 26,868** genomes, archaea **406 → 1,002**; ships with
  PICRUSt2 ≥2.6.0. Genomes more closely match study sequences across environments, improving predictions.
- **Use (§2.1, defensive and important):** this is the paper a reviewer will raise against us — "the reference
  database has been fixed." It has not, for our purposes: 26,868 genomes is still culture- and human-weighted,
  and a *bigger* reference does not fix Sun's permutation result, which is about what correlation measures, not
  about database size. Best use: acknowledge PICRUSt2-SC as the current state of the art, then note that our
  catalogue contributes 5,438 species representatives from a biome with essentially no reference genomes —
  i.e. we are supplying the input such updates depend on, not competing with the tool.

## §2.1 / pangenome decoupling

### Whelan, Hall & McInerney 2021. *Mol Biol Evol* 38:3697–3708 · doi:10.1093/molbev/msab139
`pangenome-decoupling/Whelan_2021_accessory_selection_pangenome_MBE.pdf` · 12 pp · **read: abstract (2026-08-12)** · bib: `whelan2021accessory` (**added**)
- 40 *Pseudomonas* species tested for statistically significant coincident (co-occurring/mutually avoiding)
  gene patterns. **86.7% of common accessory genes participate in at least one coincident relationship.**
  Coincident, non-vertically-inherited gene pairs are more likely to share functional categories, be
  co-transcribed, and produce interacting proteins than chance — and this is not an artefact of chromosomal
  adjacency.
- **Use (§2.1):** upgrades the Tettelin open-pangenome citation. Tettelin says gene content is not predictable
  from taxonomy; Whelan says the accessory genome is nevertheless *structured by selection into functional
  modules*. That combination is precisely our claim — "taxonomy places these organisms confidently, yet
  placement does not make their gene content predictable" — and it forecloses the reviewer response that
  accessory content is just neutral HGT noise not worth modelling.

### Guerra 2026 — the pangenome as a statistical model. *Bioinformatics Advances* 6:vbag069 · doi:10.1093/bioadv/vbag069
`pangenome-decoupling/Guerra_2026_pangenome_statistical_model_BioinfAdv.pdf` · 20 pp · **read: abstract + intro (2026-08-12)** · bib: `guerra2026pangenome` (**added**)
- Single-author review (RICAI, Colombia). Core argument: the pangenome is **not a fixed biological entity but a
  statistical model** whose open/closed verdict depends on dataset quality, curation and taxonomic resolution;
  advocates moving past the open/closed dichotomy toward pangenome fluidity and the "ecogenome".
- **Use (cite-only, §2.1 caveat).** A review in a low-profile venue — not load-bearing, and not a substitute
  for Tettelin 2005 or McInerney 2017. Its single value is one hedging sentence if we ever quantify pangenome
  openness from MAGs, since MAG incompleteness is exactly the "dataset quality" dependence it warns about.
  Do not build an argument on it.

## §2.2 — trace gas / life on air

### Bay et al. 2021 — the soil-wide benchmark. *Nat Microbiol* 6:246–256 · doi:10.1038/s41564-020-00811-w
**No PDF (paywalled).** Abstract verified via Europe PMC REST, 2026-08-12 · bib: `bay2021trace` / `bay2021desert` (both present — see duplicate warning below)
- Verbatim from the abstract: "Analysis of **40 metagenomes and 757 derived genomes** indicated that **over 70%
  of soil bacterial taxa encode enzymes to consume inorganic energy sources**. **Bacteria from 19 phyla**
  encoded enzymes to use the trace gases hydrogen and carbon monoxide as supplemental electron donors for
  aerobic respiration." A fourth phylum (Gemmatimonadota) is flagged as potentially aerobically methanotrophic.
- Thermodynamic modelling in the same paper indicates the power from H₂/CO/CH₄ oxidation is **sufficient to meet
  the maintenance needs** of the cells capable of consuming them — the same argument our §2.2 energetic balance
  makes independently for the Rub' al-Khali.
- **SCOPE TRAP — read before citing.** The ">70%" is for **inorganic energy sources generally**, not for
  trace-gas (H₂/CO) enzymes specifically; the trace-gas-specific statement is the "19 phyla" clause. Writing
  "Bay et al. found >70% of soil taxa encode trace-gas enzymes" would be a mis-citation. Our comparable
  quantity is the 34% guild union, which is a *narrower* denominator than Bay's >70%. If the §2.2 sentence
  needs a taxon-level soil benchmark, the defensible phrasing is: "trace-gas oxidation capacity is widespread
  in soils, spanning 19 bacterial phyla, and inorganic energy sources are encoded by over 70% of soil bacterial
  taxa \citep{bay2021trace}" — two claims, correctly separated.
- **Duplicate-key defect:** `bay2021trace` and `bay2021desert` are the same paper and **both are cited in
  main.tex**. Crossref confirms the seventh author is **Stefan K. Arndt**; `bay2021trace` has "Sandra K. Arndt"
  and is wrong. Not fixed here (would require a main.tex edit) — flagged for Alejandra.

### Hutchinson et al. 2025 — Antarctic soil MAGs, thermal resilience. *ISME J* (bioRxiv 2025.09.23.677255)
`trace-gas-deserts/Hutchinson_2025_antarctic_trace_gases.pdf` · **read: results (2026-08-12)** — June download, never distilled
- **554 medium-to-high-quality MAGs** from Antarctic soils, 15 bacterial + 1 archaeal phyla; Actinomycetota
  dominant (283 MAGs, 29.5 ± 11.5% of reads).
- **Group 1h/1l [NiFe]-hydrogenase in n=234 and n=119 MAGs respectively, 31 encoding both** → **322/554 = 58%**,
  spanning **11 of 15 phyla** (p7). **coxL in 60 MAGs from 5 phyla (11%)**, average 0.49 copies per genome;
  **all 60 coxL MAGs also encode a group 1h or 1l hydrogenase**, so the authors argue H₂ is the preferred donor.
  **339 MAGs (61%)** encode at least one aerotrophy gene (1h/1l, CODH or MMO). H₂ consumption is thermally
  resilient from −20 to +75 °C.
- **Use (§2.2 — the calibration we were missing, Alejandra's call whether to use it).** This is a
  genome-resolved, MAG-based, same-units comparator: Antarctic soil **58% hydrogenase / 11% CODH** against our
  **18% hydrogenase / 22% CODH**. It **inverts the naive reading**: our high-affinity hydrogenase fraction is
  *lower* than a cold desert's, while our CODH fraction is *twice* theirs. The honest headline is therefore not
  "unusually high" but **compositional**: the Rub' al-Khali guild is CO-weighted where cold deserts are
  H₂-weighted, which also sits well with our group-1l-over-1h observation (1,124 vs 445 proteins). Caveats
  before use: different binning/annotation pipelines, MAGs vs ANI95 species representatives, and their 554 is a
  much smaller and less dereplicated set — a same-pipeline re-count would be needed for a hard comparison.

### Coleine et al. 2024 — Antarctic endolithic bacteria. *Sci Total Environ* 917:170290 · doi:10.1016/j.scitotenv.2024.170290
`desert-mags-comparators/Antarctic_2023_endolithic_habitability_bioRxiv.pdf` · 12 pp · **read: results (2026-08-12)** · bib: `coleine2024endolithic` (**added**)
- **The PDF on disk is the 2023 bioRxiv preprint (2023.07.30.551190); the paper was published in 2024 in
  Science of the Total Environment with a changed title** ("potential adaptations", not "metabolic
  adaptations"). Numbers below are from the preprint and must be re-checked against the journal version before
  any of them enters the manuscript.
- Antarctic Rock Genomes Catalog: **2,636 HQ + 1,903 MQ bacterial MAGs → 2,278 species-level OTUs** at 95% ANI.
  **98.3% of species OTUs are distinct from GTDB**, i.e. 2,239 new candidate species.
- High-affinity [NiFe]-hydrogenases (1h, 1l, 1m, 2a) in **41.1% of dereplicated MAGs**, with a strong phylum
  gradient: Ca. Dormibacterota 88.9%, Eremiobacterota 80.2%, Actinobacteriota 59.1%, Gemmatimonadota 57.1%,
  Chloroflexota 53.0%, Acidobacteriota 43.9%, Verrucomicrobiota 25.8%, Planctomycetota 13.4%, Cyanobacteria
  7.5%, Bacteroidota 7.3%, Proteobacteria 6.1%, Armatimonadota 4.8%. Co-occurrence with form IE/ID RuBisCO in
  72.2% of Dormibacterota and 62.3% of Eremiobacterota MAGs.
- **Use (two sections).** (i) **§2.2**: a second genome-resolved comparator at 41.1% — with Hutchinson's 58%,
  our 18% is clearly at the low end of the cold-desert range, and the phylum breakdown shows why (their high
  fractions sit in Dormibacterota/Eremiobacterota, lineages that are minor in our catalogue). (ii) **§1.1**:
  the mirror-image contrast to ours, like Qaidam — 98.3% novel species from a reference-blank habitat against
  our 0 novel phyla or classes. Both are Antarctic/Mars-analogue rock or hyperarid systems, and both report
  novelty at the *species* rank where we report gaps at 3,755 species too. The distinction we should draw is
  rank-specific, not novelty-vs-conventionality in the abstract.

### Ehhalt & Rohrer 2009 — tropospheric H₂ budget. *Tellus B* 61:500–535 · doi:10.1111/j.1600-0889.2009.00416.x
`trace-gas-deserts/Ehhalt_2009_tropospheric_h2_cycle.pdf` · **read: abstract (2026-08-12)** · bib: `ehhalt2009h2cycle` (**added**)
- **The denominator the §2.2 \todo needs, from the abstract (p1) and Table 1 (p8):** soil uptake
  **60 +30/−20 Tg H₂ yr⁻¹** (asymmetric — this is how it is printed in both places; it is *not* 60 ± 30);
  loss to OH reaction 19 ± 5 Tg H₂ yr⁻¹; **total sinks 79 +30/−20 Tg H₂ yr⁻¹**; tropospheric burden
  **155 ± 10 Tg H₂**; tropospheric **lifetime 2.0 yr**. Sources: fossil fuel 11 ± 4, biomass burning 15 ± 6,
  ocean N-fixation 6 ± 3, land N-fixation 3 ± 2, photochemical from CH₄ 23 ± 8, photochemical from other VOC
  18 ± 7 Tg H₂ yr⁻¹.
- Soil uptake is the **dominant sink** and the largest single term in the budget; bottom-up turnovers (~75 Tg
  yr⁻¹) and top-down estimates (~105 Tg yr⁻¹) differ mainly because of it, and the authors favour the lower
  figure. Table 1 shows how unsettled it is across studies: 56 ± 41, 55, 58.3, 88 ± 11, 55 ± 8.3, 85 ± 5, and
  their own 60 +30/−20.
- **Use (§2.2):** cite as the **60 +30/−20 Tg H₂ yr⁻¹** soil sink against which our ~0.01 Tg yr⁻¹ Rub' al-Khali
  estimate is "well under one percent". Note that "60–80 Tg" as written in older repo files conflates the
  soil-uptake term with the total sink (79) — use the soil-uptake term. Also worth one clause: the sink's own
  uncertainty (+50%/−33% on the central value, and a 55–88 Tg yr⁻¹ spread across published budgets) is far
  larger than our biome's entire contribution, which is an argument for characterising unmeasured drylands
  rather than against it.

### Ocko & Hamburg 2022 — climate consequences of H₂ emissions. *Atmos Chem Phys* 22:9349–9368 · doi:10.5194/acp-22-9349-2022
`trace-gas-deserts/Ocko_2022_hydrogen_climate.pdf` · **read: abstract (2026-08-12)** · bib: `ocko2022hydrogen` (**added**)
- Hydrogen is an **indirect greenhouse gas** whose warming effect is "widely overlooked and underestimated";
  leakage across the hydrogen value chain can substantially offset the climate benefit of switching from
  fossil fuels.
- **Use (§2.2, one clause, cite-only):** why the soil H₂ sink is worth quantifying at all now. Keep it to the
  policy-relevance half-sentence; our data say nothing about leakage rates and the paper says nothing about
  deserts.

### Buzzard et al. 2022 — semiarid soil H₂ uptake vs moisture. *PeerJ* 10:e12966 · doi:10.7717/peerj.12966
`trace-gas-deserts/Buzzard_2022_semiarid_h2_uptake.pdf` · **read: abstract (2026-08-12)** · bib: `buzzard2022h2uptake` (**added**)
- Sonoran Desert / semiarid urban ecosystem; North American Monsoon (June–September) supplies roughly half the
  annual precipitation. Soil moisture dynamics strongly modulate atmospheric-H₂ uptake, with uptake peaking at
  intermediate moisture and suppressed when the soil is dry or saturated.
- **Use (§2.2 caveat, load-bearing for the flux limit):** the mechanistic reason our capacity number must not
  be read as a flux. In a hyperarid sand sea the moisture optimum is rarely met, so realised uptake will sit
  below the diffusion-limited estimate for most of the year. Pairs with Jordaan 2020 (hydration stimulates
  H₂-oxidiser activity in desert soils) — together they are the "capacity ≠ activity ≠ flux" citation pair.

### Ortiz et al. 2021 — multiple energy sources, Antarctic desert soils. *PNAS* 118:e2025322118 · doi:10.1073/pnas.2025322118
`trace-gas-deserts/Ortiz_2021_antarctic_energy.pdf` · **read: page 1 (2026-08-12)** · bib: `ortiz2021antarctic` (**added**)
- Trace-gas oxidation and chemolithotrophy, rather than photosynthesis, dominate the energy economy of
  Antarctic desert soils. Greening/Cowan/Makhalanyane consortium (16 authors, verified via Crossref).
- **Use (§2.2, cite-only alongside Ji 2017):** the OA companion to the paywalled founding paper, for the claim
  that trace gases are a primary rather than supplementary energy source in desert soil.

## §2.3 — UV, desiccation, salt

### Sharda, Badrinarayanan & Seshasayee 2020 — bacterial NHEJ distribution. *Genome Biol Evol* 12:2450–2466 · doi:10.1093/gbe/evaa223
`survival-genome/Sharda_2020_bacterial_NHEJ_evolution_GBE.pdf` · 17 pp · **retrieved + read this session via Europe PMC (2026-08-12)** · bib: `sharda2020nhej` (**added**)
- **The comparator the §2.3 \todo needs.** Comparative genomics across **~6,000 bacterial genomes**:
  "We found NHEJ in only **~1,300 (22%)** genomes studied here"; restated in the Discussion as "we confirm the
  nonubiquitous nature of NHEJ across the bacterial domain … with **22% bacteria coding for it**."
- Supporting detail: 75% of NHEJ-positive bacteria carry Ku and LigD within 10 kb of each other; 84% carry a
  single Ku copy; **138 (2.3%) encode Ku without LigD and 619 (10.3%) LigD without Ku** — so the operational
  definition matters, and a Ku-only or LigD-only call inflates prevalence. **96% of archaeal genomes** encode
  a LIG domain but no Ku.
- **The mechanism, and the caveat that must travel with the number.** NHEJ presence correlates with **large
  genome size and/or slow growth rate**, genome size being the dominant correlate; the pathway is sporadically
  distributed and spread by horizontal transfer. Their 22% is over mostly **cultured, complete genomes**, which
  are weighted toward small, fast-growing organisms — the opposite of the large-genome slow-growing soil
  Actinomycetota that dominate our catalogue.
- **Use (§2.3):** replaces the bare "a high prevalence for soil bacteria" with a real contrast — **74.5% of our
  5,438 species representatives versus 22% of ~6,000 sequenced bacterial genomes**, a 3.4-fold difference. But
  the sentence should not attribute the excess to desiccation selection alone, because Sharda supplies a
  competing explanation we cannot exclude: our catalogue is composition-biased toward exactly the genome-size
  and growth-rate class in which NHEJ is expected. The defensible claim is that the capacity is near-universal
  in this biome and far above the sequenced-genome baseline, with the lineage-composition confound stated. If
  Alejandra wants to strengthen it, the tractable test is internal: NHEJ presence versus genome size within
  our own catalogue, and whether the 74.5% survives conditioning on phylum.
- Pitcher et al. 2007 (`pitcher2007nhej`, already keyed) remains the correct **causal** citation — NHEJ protects
  stationary-phase mycobacteria against desiccation — but it supplies no prevalence number. Cite the two
  together: Sharda for the baseline, Pitcher for why it matters here.

### Moore et al. 2023 — desiccation metabolism of an Atacama cyanobacterium. *Sci Rep* 13:15767 · doi:10.1038/s41598-023-41879-8
`survival-genome/Atacama_2023_desiccation_pathways_cyanobacterium_SciRep.pdf` · 13 pp · **read: abstract + discussion (2026-08-12)** · bib: `moore2023desiccation` (**added**)
- Genome + genome-scale metabolic model (iGd895) for *Gloeocapsopsis dulcis* AAB1, an extremely xerotolerant
  Atacama hypolith that accumulates trehalose and sucrose at water activity as low as **a_w 0.4**.
- **The §2.3 caveat the \todo asks for, with its provenance (p10):** "Baran et al. proposed that cyanobacteria
  may respond to changes in salinity by **converting glycogen into compatible solutes, rather than synthesizing
  them de novo**, and this mechanism may also be relevant for *G. dulcis*." The glycogen-accumulation reaction
  appears **only in the water-constrained simulation** (p9), and trehalose is produced only under constrained
  water and only via trehalose phosphorylase (TreP, EC 2.4.1.64) rather than trehalose synthase (p5).
- Also: quinol oxidase appears only in the constrained simulation, implicating **metabolic water production**
  as a survival mechanism (p10).
- **Three caveats that must travel with this citation.** (i) It is **secondhand** — Moore attributes the
  glycogen→solute idea to Baran et al. and extends it; if the claim becomes load-bearing, Baran should be cited
  directly. (ii) It is a **flux-balance model prediction**, not measurement. (iii) The organism is a
  **photoautotrophic cyanobacterium under salinity/desiccation**, not a soil heterotroph.
- **Use (§2.3, one sentence):** the honest form is a limitation on our osmolyte capacity result — pathway
  presence does not establish de novo synthesis under stress, since stored carbon can be reallocated to
  solutes instead. That weakens nothing we claim (we claim capacity, not flux) and pre-empts the reviewer who
  asks whether trehalose genes are actually used.

### Oren 2008 — microbial life at high salt. *Saline Systems* 4:2 · doi:10.1186/1746-1448-4-2
**No PDF yet (OA, retrievable).** Abstract verified via Europe PMC REST, 2026-08-12 · bib: `oren2008halophiles` (**added**)
- The salt-in citation §2.3 needs, from the abstract: halophiles use **two strategies**. The first is
  **accumulation of molar concentrations of KCl**, which "requires adaptation of the intracellular enzymatic
  machinery… The proteome of such organisms is highly acidic, and most proteins denature when suspended in low
  salt. Such microorganisms generally cannot survive in low salt media." The second is to exclude salt and
  accumulate organic compatible solutes, used by most halophilic Bacteria and by halophilic methanogens.
  Halobacterium and relatives require **over 100–150 g/l salt** for growth and structural stability.
- **Use (§2.3):** exactly supports our sentence that 269 Halobacteriota species almost entirely lack the
  bacterial osmolyte pathways (ectoine 5/269, GG and MG in none) "consistent with their primary salt-in
  osmoregulation". Oren's acidic-proteome point is a bonus: it means the absence is a *coherent physiological
  syndrome*, not a MAG-completeness artefact — worth one clause, because incompleteness is the obvious
  alternative explanation for missing pathways.
- Optional companion: **Oren 1999**, *Microbiol Mol Biol Rev* 63:334–348, doi:10.1128/mmbr.63.2.334-348.1999
  (verified, not OA in Europe PMC) argues KCl accumulation is **energetically cheaper** than organic-solute
  synthesis. Only needed if §2.3 makes a cost argument.
- Mattimore & Battista 1996 (`mattimore1996radioresistance`) and Slade & Radman 2011 (`slade2011oxidative`) are
  **already in references.bib with correct metadata** — the §2.3 \todo for these two is a drafting task, not a
  bibliography gap.

## §2.4 — BGCs

### Andreani-Gerard, Cambiazo & González 2024 — Atacama soil BGCs. *mSphere* 9:e00192-24 · doi:10.1128/msphere.00192-24
`metabolism-bgc/AndreaniGerard_2024_Atacama_soil_BGCs_mSphere.pdf` · 15 pp · **read: abstract (2026-08-12)** · bib: `andreanigerard2024bgc` (**added**)
- Six bacterial communities along an Andean altitudinal gradient in the Atacama; **38 MAGs → 168 BGCs**, mostly
  NRPS, RiPPs and terpenes, concentrated in Acidobacteriota and Proteobacteria. Six MAG clusters by
  core-biosynthetic-gene composition, three dominated by a single phylum, two also site-specific.
- Accessory-gene associations: **RND multidrug efflux pumps with RiPPs**, **TonB iron transporter with
  terpenes**.
- **Their stated conclusion, which is ours too:** MAGs "did not show an obvious relationship between the
  composition of their core biosynthetic potential and taxonomy or geographic distribution."
- **Use (§2.4):** the hot-desert counterpart to Waschulin 2021's Antarctic long-read set, and a *short-read
  MAG-based* precedent — i.e. methodologically closer to us than Waschulin. Two concrete uses: (i) the
  taxonomy-independence of BGC composition is a second, independent line for the decoupling thesis, this time
  in secondary metabolism; (ii) their 168 BGCs from 38 MAGs sets the scale a reviewer will expect, so our
  antiSMASH census should be reported per-genome, not only as a total. Their resistance-transporter/BGC
  co-occurrence also connects §2.4's BGC census to the intrinsic-resistome half of the same subsection.

## §3.1 / contamination and specificity

### Carini et al. 2016 — relic DNA in soil. *Nat Microbiol* 2:16242 · doi:10.1038/nmicrobiol.2016.242
`contamination-armor/Carini_2016_relicDNA_soil_NatMicro.pdf` · 6 pp · **read: abstract (2026-08-12)** · bib: `carini2016relicdna` (**added**)
- Viability PCR (propidium monoazide) across a wide range of soils: **on average 40% of both prokaryotic and
  fungal DNA was extracellular or from cells that were no longer intact**. Relic DNA **inflates observed
  richness by up to 55%** and misestimates taxon relative abundances. It is not present in measurable amounts
  in all soils — more likely in soils with **low exchangeable base cations**, and its removal changes community
  structure most in **high-pH soils**.
- **Use (§3.1 and Methods, defensive):** the desk-reject risk our catalogue must answer, and the high-pH
  sensitivity makes it specifically relevant to alkaline desert sand. The honest position is that we cannot
  exclude a relic-DNA contribution to DNA-based abundances, but relic DNA does not manufacture the signals we
  actually claim: it cannot produce coherent 95%-ANI genome bins with consistent coding density, nor the
  purifying selection we measure in the dark families, nor the phylum-specific salinity coupling. Say this
  once, explicitly, rather than leaving it to a reviewer.

### Rabasco et al. 2025 — kitome of contemporary DNA extraction kits. bioRxiv 2025.09.23.677068
`contamination-armor/Kitome_2025_DNAextractionkits_bioRxiv.pdf` · 8 pp · **read: abstract (2026-08-12)** · bib: `rabasco2025kitome` (**added**)
- Updated characterisation of contamination in DNA extraction kits in current use, motivated by the fact that
  the pioneering surveys (Salter 2014, Karstens 2019) predate today's kits and chemistries; provides
  kit-selection guidance for low-biomass work. Preprint, not peer reviewed.
- **Use (cite-only, Methods):** cite alongside `salter2014contamination` and the 2025 low-biomass consensus
  statement to show the kit-contamination concern was addressed against *current* reagents rather than a
  decade-old survey. Flag as a preprint in the bib entry.

## Framing / comparators (cite-only)

### Schultz et al. 2026 — extreme Arabian environments and their microbiomes. *Extremophiles* 30:13 · doi:10.1007/s00792-026-01426-5
`desert-mags-comparators/Arabian_2026_extreme_microbiomes_astrobiology_Extremophiles.pdf` · 19 pp · **read: abstract (2026-08-12)** · bib: `schultz2026arabian` (**added**)
- Review from the KAUST/Rosado group (Schultz, García-Martínez, Altalhi, Kontis, dos Santos, Rosado — verified
  via Crossref). Argues Saudi Arabian deserts and salt flats are under-explored Mars analogues alongside the
  Atacama, McMurdo Dry Valleys and Río Tinto, and that their astrobiological potential is largely untapped.
- **Use (§1.1 or Discussion, one sentence):** the contemporary, KAUST-authored statement that Arabian extreme
  environments are under-sampled, which is the premise of our catalogue and better than asserting it ourselves.
  Handle the astrobiology angle lightly — our paper is a soil MAG catalogue, and a Mars-analogue frame invites
  reviewers to ask for habitability analyses we do not do.

### Zhou et al. 2025 — microbial life-history strategies and multifunctionality across aridity. *PNAS* 122:e2511071122 · doi:10.1073/pnas.2511071122
`dryland-ecology/Zhou_2025_lifehistory_aridity_multifunctionality_PNAS.pdf` · 8 pp · **read: abstract + significance (2026-08-12)** · bib: `zhou2025lifehistory` (**added**)
- 474 soil samples across a continental-scale aridity gradient. Microbial life-history strategies (high yield,
  resource acquisition, **stress tolerance**) covary with ecosystem multifunctionality, and there is an
  **aridity threshold beyond which plot-level multifunctionality declines sharply**.
- **Use (§2.3 or Discussion):** the trait-based framing that connects our stress-tolerance gene inventory to
  ecosystem consequence, and the current citation for aridity thresholds (fresher than Maestre 2015). One
  caution: their stress-tolerance trait proxies are not our KO panel, so use it for framing, not calibration.

### Li, He & Wang 2026 — small proteins at full ocean depth. *iScience* 29:114585 · doi:10.1016/j.isci.2025.114585
`smorfs-smallproteins/Li_2026_smallproteins_fulloceandepth_iScience.pdf` · 22 pp · **read: abstract + highlights (2026-08-12)** · bib: `li2026smallproteins` (**added**)
- Marine water-column small-protein catalogue: **4,307 high-confidence clusters** across depths; **87% are
  marine-specific**, **~70% have unknown domains**; includes prophage-associated and broadly conserved
  (≥5 phyla) clusters, with zone-specific transcription linking small proteins to deep-sea adaptation.
- **Use (§1.2, the closest structural analogue to our claim):** a single-biome small-protein catalogue that
  reports both a high habitat-specific fraction (87%) and a high unknown fraction (~70%) — precisely the shape
  of our dark-family specificity result (21.4% union-absent; 47.4% absent from GEM). Two things it buys us:
  (i) it shows single-biome specificity fractions of this magnitude are reported and accepted elsewhere, and
  (ii) their transcriptional evidence marks the gap in ours — they show expression, we show only capacity,
  which is the natural place to point at the companion annotator paper. Complements `duan2024gmsc` (global
  scale) with a habitat-scale comparator.

---

## Corrections and defects found while ingesting (2026-08-12)

1. **`bay2021trace` and `bay2021desert` are duplicate entries for the same paper, and both are cited in
   `main.tex`.** Crossref gives the seventh author as **Stefan K. Arndt**; `bay2021trace` says "Sandra K.
   Arndt" and is wrong, `bay2021desert` is correct. Needs a manuscript edit (consolidate to one key) — not
   done here.
2. **Delgado-Baquerizo networks paper: the previously promised key `delgadobaquerizo2019networks` carries the
   wrong year.** Verified: *Land Degradation & Development* **31(3):346–355**, doi:10.1002/ldr.3453, issue year
   **2020** (Crossref "issued" 2019 = online-first). Added as **`delgadobaquerizo2020networks`**. Nothing cited
   the old name, so no manuscript impact.
3. **Coleine endolithic: the PDF on disk is a preprint superseded by a 2024 journal article** with a changed
   title (doi:10.1016/j.scitotenv.2024.170290). Bib entry points at the journal version; the numbers in the
   note above come from the preprint and need re-checking before use.
4. **The five bib keys promised in the batch-1 footer had never been added.** All five are now in both bib
   files: `liu2025qaidam`, `contamination2025guidelines`, `delgadobaquerizo2020networks` (renamed, see 2),
   `zhang2022envselection`, `waschulin2021bgc`. Their metadata was verified against Crossref this session —
   note that `contamination2025guidelines` is a 30-author consensus statement led by **Fierer** and closing
   with **Greening**, so a first-author citation would be "Fierer et al. 2025", not an anonymous consensus.
5. **`ji2017atmospheric` is still PDF-less and paywalled** (*Nature* 2017, doi:10.1038/nature25014). It is
   cited in §2.2 for the founding claim. `ortiz2021antarctic` (now keyed, OA, same consortium) covers the same
   ground and can carry the claim if a reviewer challenges an unread citation.

_Still not resolvable from PDFs on disk: Makhalanyane 2015 (hot-desert light regimes, §2.3 photolyase framing)
and Roberts 2005 (compatible-solute hierarchy, §2.3 salt half) — both remain open rows in TO_DOWNLOAD.md._

---

## Citation audit — session 8 (2026-08-18). Sources verified for claim support, not yet distilled

Read to check whether cited papers support the claims attached to them in `paper/main.tex`
(§2.1, §3.1, Discussion, PI scaffold). Verdicts and fixes were returned to the main session as text.

### Richardson et al. 2023 — MGnify 2023. *NAR* 51:D753--D759 · doi:10.1093/nar/gkac1080
no PDF in repo · **read: abstract + p.D754 (2026-08-18)**, fetched via `europepmc.org/articles/PMC9825492?pdf=render` · bib: `richardson2023mgnify` (already present)
- Verbatim, p.D754: "Currently, 297 different biomes are represented in the database with over half of
  MGnify's analyses originating from merely nine of them: human-faecal, -oral, -digestive system, -skin,
  and unspecified human; marine; soil; mammalian digestive systems; and mixed biome samples."
- Scale: 475,390 analyses over 343,695 distinct samples in 4,601 studies; protein database "exceeds 2.4
  billion non-redundant sequences predicted from metagenomic assemblies".
- The paper never names **MGnify90**: the 717,738,164-sequence 90%-clustered subset used in §3.1 is a
  derived product, so Methods M.12 should state release and clustering rather than lean on this cite.
- **Use (§3.1):** the citable skew statement is *human-associated and marine dominate*, and **soil is one
  of the nine** -- so no fix may claim soil or desert absence from MGnify. See defect 6 below.

### Aplakidou et al. 2025 — "Quadrupling the protein family space with global metagenomics". bioRxiv
no PDF in repo · **read: abstract + Results (2026-08-18)**, fetched `biorxiv.org/content/10.64898/2025.12.21.695754v1.full` · bib: **owed** (`aplakidou2025quadrupling`)
- **DOI is `10.64898/2025.12.21.695754`, not `10.1101/...`** (bioRxiv's new prefix). The `\todo` in
  `main.tex` implies the 10.1101 form, which Crossref resolves as "Resource not found".
- 40.3 billion sequences over 40,446 metagenomes, 9,540 metatranscriptomes and 539 M proteins from 167,415
  reference genomes yield **608,258 previously uncharacterized families with >=100 members** (and 6.5 M with
  >=25), "none matching known Pfam domains or reference proteins".
- Verbatim: "Of these families, **385,791 (63.4%) are restricted to a single ecosystem**: primarily
  Environmental (57.5%), with smaller contributions from Host-Associated (4.9%) and Engineered (0.7%)."
  Single-ecosystem subtypes: Terrestrial 89,380 (14.6%), Freshwater 122,073 (20%), Marine 109,286 (17.9%).
- Two calibrations for §3.1: "ecosystem" = GOLD categories resolved into eight subtypes, not biome in a
  finer sense; and the 608,258 set was built under the criterion ">=100 members and at least 5% of samples
  originating from a single ecosystem type", which bakes a mild ecosystem skew into the denominator.
- **Use (§3.1):** carries the reframe -- biome restriction is the norm for a novel family, so our 28.7% is
  unremarkable *as a fraction* and the contribution is *which* biome the references are missing.

### Bay et al. 2021 — desert aridity gradient. *ISME J* 15:3339--3356 · doi:10.1038/s41396-021-01001-0
no PDF in repo · **read: abstract + methods/data availability (2026-08-18)**, EPMC PMC8528921 · bib: **owed** (`bay2021aridity`)
- **Distinct from `bay2021trace`** (*Nat Microbiol* 6:246--256). Different paper, different sites.
- Topsoils and biological soil crusts along a steep aridity gradient in the **Negev Desert (Israel)** across
  sub-humid, semi-arid, arid and hyper-arid zones. "All metagenomes sequenced for this project are available
  on the Sequence Read Archive under BioProject accession number **PRJNA657906**" -- matches the
  `desert_autotroph` leg (7,200,701 seqs) in `specificity_clean/README.md`.
- Cell-specific atmospheric H2 consumption rises **143-fold** along the gradient, correlating with increased
  abundance of high-affinity hydrogenases; chemosynthesis dominates in arid and hyper-arid soils,
  photosynthesis in biocrusts.
- **Use (§3.1):** the second published desert comparator, currently uncited. Also a legitimate §2.2 support
  for chemosynthesis dominating at the dry end -- but it is a *rate* measurement, so keep our capacity/flux
  fence. Note the Negev is the geographically nearest of the five comparators to the study area.

### Uritskiy et al. 2019 — Atacama halite after rainfall. *ISME J* 13:2737--2749 · doi:10.1038/s41396-019-0468-y
no PDF in repo · **read: abstract (2026-08-18)**, EPMC PMC6794293 · bib: `uritskiy2019atacama` (already present)
- "an extremophile community inhabiting **halite (salt rocks)** in the Atacama Desert"; 4-year longitudinal
  shotgun metagenomics of the response to a catastrophic rainfall.
- **This is endolithic salt rock, not soil**, which collides with §3.1's "the largest desert **soil**
  collections available". It is also only one of five comparator legs (6,193,359 of 147,531,178 seqs).
- Bib byline defect: Crossref gives **Getsin, Samantha** and **Glass, Brian**; `references.bib` has Samuel
  and Benjamin.

### D'Costa et al. 2011 — "Antibiotic resistance is ancient". *Nature* 477:457--461 · doi:10.1038/nature10388
no PDF, paywalled (no PMCID; PMID 21881561) · **read: abstract only (2026-08-18)** · bib: `dcosta2011ancient` (added by main session)
- Targeted metagenomics of authenticated ancient DNA from **30,000-year-old Beringian permafrost**: genes for
  beta-lactam, tetracycline and glycopeptide resistance; structure/function on the complete VanA element.
  "antibiotic resistance is a natural phenomenon that predates the modern selective pressure of clinical use."
- **Use (§2.4):** supports *pre-anthropogenic*, **not** the soil-competition mechanism. Permafrost sediment,
  not contemporary soil -- pair it with Van Goethem rather than letting it carry the sentence alone.

### Van Goethem et al. 2018 — pristine Antarctic soil resistome. *Microbiome* 6:40 · doi:10.1186/s40168-018-0424-5
no PDF in repo (OA: `europepmc.org/articles/PMC5824556?pdf=render`) · **read: abstract (2026-08-18)** · bib: `vangoethem2018antarctic` (added by main session)
- Verbatim background: "Soil bacteria naturally produce antibiotics as a competitive mechanism, with a
  concomitant evolution, and exchange by horizontal gene transfer, of a range of antibiotic resistance
  mechanisms." 17 pristine Mackay Glacier soils; **177 naturally occurring ARGs**, mostly efflux pumps, plus
  aminoglycoside, chloramphenicol and beta-lactam inactivation. ARG abundance correlates negatively with
  species richness (r = -0.49, P < 0.05); **no mobile genetic elements flanking ARGs** => "functional
  efficient historical genes that have since been vertically inherited over generations".
- **Use (§2.4):** the closer match of the two for L190 -- it is the paper that makes intrinsic resistance a
  *host-genome trait in pristine soil* rather than merely an ancient one.

### Lennon, Muscarella, Placella & Lehmkuhl 2018 — relic DNA. *mBio* 9(3):e00637-18 · doi:10.1128/mbio.00637-18
metadata only (Crossref, 2026-08-18) · bib: **owed** · not to be confused with `lennon2011seedbanks`
- Verified as the paper the `\todo` intends. The 2017 bioRxiv preprint `10.1101/131284` has a different title
  ("How, when, and where relic DNA **biases estimates of** microbial diversity") and a different author list
  (Muscarella SA, not Placella SA) -- cite the mBio version.

## Defects found in the citation audit (2026-08-18)

6. **`richardson2023mgnify` is mis-cited at `main.tex` L210 for "MGnify90 itself is dominated by aquatic and
   engineered biomes".** The paper documents the opposite: over half of MGnify's analyses come from nine
   biomes led by human-associated, plus marine, **soil**, mammalian digestive and mixed; "engineered" is not
   among them, and the passage counts analyses, not protein-database sequences. Provenance of the error: it
   descends from June's biome-of-**hits** breakdown (83.5% aquatic), which
   `specificity_clean/README.md` records as **not recomputable** from the current output -- and which
   measured where *our* families hit, never MGnify90's composition. Fix is a rewrite, not a re-cite, and it
   must stay inside the workbench's provenance-level rule ("no Arabian or Saharan site in any named desert
   comparator"), never "absent from all catalogues".
7. **`uritskiy2019atacama` alone cannot carry `main.tex` L210.** One cite for five comparators; the paper is
   halite endoliths rather than soil; and the negative "no Arabian or Saharan site appears in any of them"
   is a provenance fact no citation can carry -- three of the five legs are unpublished JGI studies. The
   superlative "the largest desert soil collections **available**" is also unsupportable while
   `DISCUSSION_NOTES.md` records Leung 2025 (PRJNA832417, 911 MAGs incl. **Kuwait**) outside the five.
8. **`matchado2024limits` is doing soil work it cannot do** (L168, L466, L438). Its experiments are human gut
   disease cohorts plus simulated data; its soil sentence is a second-hand summary of Djemiel et al. (its
   ref 45). `sun2020inference` is the paper with soil datasets and should be co-cited wherever the claim is
   about soil-microbiome inference.
9. **Box 2 (L1362): the contamination trio is attached to the wrong alternative.**
   `salter2014contamination,eisenhofer2019lowbiomass,davis2018decontam` are all the correct papers (Crossref-
   verified) but sit on "Reference bias overstated", while "Signal = contamination/kit-ome" one clause
   earlier carries no citation. Box 2's "436 matched samples" is also stale against §2.1's 448 -- scaffold
   text, retained by design, noted only.

---

## Figure-design survey (2026-08-23) — main-figure inventory of high-tier comparators

Purpose: identify the display genres that anchor comparable papers, to upgrade our display set (§09-figures).
Figure *content and type* only; captions read from the PDF, five anchor pages rendered and inspected.
Main figures only (Extended Data listed only where it carries a genre we lack).

### Rodríguez del Río et al. 2024 — FESNov, *Nature* 626:377 (`structure-darkmatter/`, on disk)
- Fig 1: workflow schematic (a) + three distribution panels — dN/dS (b), % identity (c), protein length vs
  known families (d). Pipeline and its sanity checks in one figure.
- Fig 2: presence/absence matrix — FESNov families (cols) x GTDB orders (rows), KEGG-pathway linked.
- Fig 3: **circular GTDB order-level tree with two annotation rings** (red bars = unknown families per
  genome; blue ring = uncultivated ratio), outlier clades labelled with counts; plus two small line/bar
  panels (b, c: habitat breadth and rank-confinement, dual y-axes with mobility lines). Big anchor + small
  supporting panels in one float.
- Fig 4: gene-neighbourhood schematics (arrow-gene diagrams) for four synapomorphic families.
- Fig 5: differential-abundance panel + classifier performance (disease application).

### Pavlopoulos et al. 2023 — functional dark matter, *Nature* 622:594 (PMC10584684; **not on disk**)
- Fig 1: clustering-overview schematic (a) + rarefaction curves (b) + cluster-composition bars (c,d).
- Fig 2: **UpSet plot** of novel families across 8 GOLD ecosystem types.
- Fig 3: UpSet of domain-level taxonomy + occurrence across MAGs.
- Fig 4: **filtering cascade drawn as a labelled flow diagram with counts at every step, integrated in the
  same panel with two structure galleries** — (b) 12 candidate novel-fold models coloured by pLDDT with a
  legend; (c) query/target structural superpositions annotated with family and SCOPe/PDB IDs.
- Extended Data adds: circos/chord plots per ecosystem, a **world map of sample locations**, and
  per-example composites (3D model | gene-neighbourhood network | ecosystem distribution) side by side.

### Durairaj et al. 2023 — new families and folds, *Nature* 622:646 (PMC10584680; **not on disk**)
- Fig 1: workflow schematic for classifying functionally dark proteins (UniProt -> AFDB).
- Fig 2: **large-scale sequence-similarity network**, 4.27M nodes / 10.3M edges reduced to 688,852
  communities (Cosmograph layout) — the "map of protein space" genre.
- Figs 3, 4: per-case deep dives — high-resolution CLANS network + structures + genomic context for one
  new family (Fig 3) and one new toxin-antitoxin superfamily (Fig 4).
- Fig 5: outlier-diagnosis panels (brightness, shape-mer diversity, length) separating fragments and
  repeats from genuine new folds — i.e. the false-positive anatomy shown, not just quoted.

### Yeo et al. 2025 — AFESM, bioRxiv (`structure-darkmatter/`, on disk)
- Fig 1: database composition (a) + clustering-pipeline schematic (b).
- Fig 2: within-cluster structural similarity distributions + taxonomic projection.
- Fig 3: biome analysis — LCB label breakdown, extreme-environment-specific clusters.
- Fig 4: novel-domain identification workflow with quality gates (segmentation -> CATH SSG5 -> metrics).
- Fig 5: **domain-combination analysis** — fraction of multi-domain proteins with novel co-occurring CATH
  topologies. The "novelty is architectural" display.
- Supp adds superposition galleries of environment-specific clusters vs their nearest non-specific match.

### Lau et al. 2024 — TED, *Science* 386:eadq4946 (`structure-darkmatter/`, on disk)
- Fig 1: overall workflow with counts.
- Fig 2: CATH-hierarchy classification — top-superfamily comparison + proportion-matched bars.
- Fig 3: **structure gallery** of high-symmetry domains and "extruded repeats" (examples, not statistics).
- Fig 4: novelty score vs sequence-cluster size scatter (a) + taxonomic distribution of novel clusters (b).
- Fig 5: enrichment log2FC + alignment-procedure schematic for internal structural repeats.

### Coelho et al. 2022 — GMGC, *Nature* 601:252 (PMC7613196; **not on disk**)
- Fig 1: **world map with sample dots sized by read count, with the clustering cascade drawn as arrows and
  counts around the map** (a) + habitat-sharing ribbon/alluvial plot (b) + gene-accumulation curves (c) +
  family-size distribution with a stacked proportion bar (d). Schematic, geography and three data views in
  one float — the model for a design figure that is not a box diagram.
- Fig 2: per-habitat density curves of conspecific gene counts.
- Fig 3: gene-prevalence histograms on log-log axes.
- Fig 4: operon preservation + fraction under positive selection vs detection count.

### Nayfach et al. 2021 — GEM, *Nat Biotechnol* 39:499 (PMC8041624; **not on disk**)
- Fig 1: **environmental and geographic distribution** of 52,515 MAGs (map + habitat composition).
- Fig 2: species-level clustering against >500,000 reference genomes; novelty vs references.
- Fig 3: **phylogenetic tree of 43,979 OTUs with green branches for new lineages and four concentric strip
  charts** (newly identified / uncultivated / environment / MAG count), plus (b) phylogenetic-diversity
  gain bars per phylum and environment.
- Fig 4: BGC frequency heat-strip by phylum and habitat + one annotated largest-BGC diagram.
- Fig 5: host-virus connectivity bars by phylum.

### Liu et al. 2025 — Qaidam Basin, 1,773 MAGs, *Sci Data* 12:1795 (`desert-mags-comparators/`, on disk)
- Fig 1: (a) **map of sampling sites**; (b) workflow schematic.
- Fig 2: completeness-vs-contamination scatter (a) + taxonomic-novelty bars by rank (b) + MAGs per phylum (c).
- Fig 3: same three-panel template for the viral MAGs.
- Fig 4: **iTOL circular phylogenomic tree of all 1,773 MAGs** — phylum colour strip + four concentric bar
  rings (completeness, contamination, genome size, GC), archaeal branch marked, phylum legend at left.
  This is the whole-catalogue tree figure at our scale, in a data-descriptor journal.

### Ortiz et al. 2021 — Antarctic energy sources, *PNAS* 118:e2025322118 (`trace-gas-deserts/`, on disk)
- Fig 1: community abundance/diversity/composition + phylogenetic divergence (boxplots + composition).
- Fig 2: **signature-gene heatmap** of metabolic potential (inorganic / organic / light) across samples.
- Fig 3: (a) maximum-likelihood **group-1 [NiFe] hydrogenase tree** with the group-1l clade highlighted,
  plus supporting expression/abundance panels.
- Fig 4: measured oxidation rates and energy yields (boxplots, gas chromatography).

### Hutchinson et al. 2025 — Antarctic trace gases, *ISME J* (`trace-gas-deserts/`, on disk)
- Fig 1: **map of Dronning Maud Land with sampling points** + in situ/ex situ oxidation panels.
- Fig 2: (A) **Antarctica-wide map of metagenome sampling locations**; (C) trace-gas gene abundance across
  metagenomes and phyla.
- Fig 3: oxidation rates across a temperature range.
- Figs 4-5: **predicted-rate maps of ice-free Antarctica** (baseline and seasonal/emissions scenarios) —
  gene survey extrapolated onto geography.

### Ji et al. 2017 — *Nature* 552:400 — figure inventory NOT retrievable
No PMC record (isOA: N), nature.com `/figures` returns a 303 to `idp.nature.com` from this node. Already on
TO_DOWNLOAD for its text; the figure list stays unknown rather than guessed. Ortiz 2021 and Hutchinson 2025
(both on disk, same consortium/biome) carry the trace-gas display genre instead.

### Bay et al. 2021 — *Nat Microbiol* 6:246 — figure inventory NOT retrievable
No PMCID; paywalled. Numbers already distilled above from the abstract; figures unseen.

### Filing trap found (2026-08-23)
`trace-gas-deserts/Bay_2020_chemosynthesis_poles.pdf` is **not a Bay paper**. Page 1 identifies it as
**Ray, Zhang, Terauds, Ji, Kong & Ferrari 2020, "Soil Microbiomes With the Genetic Capacity for Atmospheric
Chemosynthesis Are Widespread Across the Poles...", *Front Microbiol* 11:1936, doi:10.3389/fmicb.2020.01936**
(5 figures: site map, qPCR abundance, gene distribution, ordination/ANOSIM, random-forest predictors).
Same trap class as `Constant_2016_streptomyces_h2.pdf` (= Liot & Constant 2016). Neither Bay 2021 paper is
on disk. Rename or re-file before anything cites it.

### Recurring genres, ranked by how often they anchor these papers
1. **Annotated circular phylogeny with concentric data rings** — FESNov 3, GEM 3, Qaidam 4.
2. **Workflow/cascade schematic fused with its own data panels** — Pavlopoulos 1 and 4, Durairaj 1,
   AFESM 1 and 4, TED 1, GMGC 1a, Qaidam 1b. Never a bare box diagram: counts and distributions are inside it.
3. **Geographic map of samples, sometimes with a predicted-quantity surface** — GMGC 1a, GEM 1,
   Qaidam 1a, Hutchinson 1/2/4/5, Pavlopoulos ED 7.
4. **Protein-structure gallery / superposition** — Pavlopoulos 4b-c, TED 3, Durairaj 3-4, AFESM supp.
5. **Set-overlap displays (UpSet, circos, alluvial, network)** — Pavlopoulos 2-3, GMGC 1b, Durairaj 2.
6. **Marker-gene tree with clade calls** — Ortiz 3A, Hutchinson (supp), FESNov 2 matrix-on-tree.
7. **Distributions/scatters as supporting panels, rarely as the anchor** — this is where our current set
   sits almost entirely.

### What reads as "Nature quality" beyond chart type (from the five rendered anchors)
- **One dominant panel plus small satellites.** FESNov 3 and GMGC 1 give ~60% of the area to one rich
  anchor; the remaining panels are deliberately small and simple. Our floats give equal weight to 2-3
  panels of the same type.
- **Numbers live inside the graphic.** Pavlopoulos 4a prints the count at every filtering step; GMGC 1a
  prints 2,007,736,046 -> 966,108,540 -> 302,655,267 along the arrows. The reader never leaves the figure.
- **Multiple encodings stacked on one geometry.** Qaidam 4 and GEM 3 hang 4-5 variables off one tree
  (colour strip, bar rings, branch colour); FESNov 3a uses two rings plus labelled outliers.
- **Named exemplars alongside the statistics.** TED 3, Durairaj 3-4, Pavlopoulos 4c show specific
  structures with IDs. Statistics establish the claim; exemplars make it concrete.
- **Restrained palette with semantic colour.** Categorical hues for taxa only; a single sequential ramp
  for one continuous variable (pLDDT, aridity); grey for controls/background. No decorative colour.
- **Legends inside the plotting area**, small type, and dual y-axes used freely (FESNov 3b-c).

### Asset check for the genre upgrade (verified 2026-08-23, workbench paths)
- Catalogue ring-tree: **no de novo tree of the 5,438 exists**, but an R232 GTDB-Tk 2.7.2 user MSA does —
  `M/rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/gtdbtk_2.7.2_r232/run/align/gtdbtk.{bac120,ar53}.user_msa.fasta.gz`
  (4,920 + 621 = 5,541 sequences; **5,251 of the 5,438 reps present, 96.6%**). A ring-tree is a subset +
  FastTree away, not a new GTDB-Tk run. The 187 absent reps need accounting before the tree is quoted.
- Structure gallery: `M/06_esmfold_novelfold/pdb/` (68,543 PDB files); `module load pymol/3.1.6` or
  `chimerax/1.7` available on IBEX (no py3Dmol/nglview in the base python).
- Map: `M/metadata/site_coords.tsv` (lat/lon, 81 rows incl. non-canonical sites; 60 canonical).
- Marker trees: `M/06_atmospheric/classify/{hydB,coxL}.nwk` (already the source of the §2.2 option figure).
