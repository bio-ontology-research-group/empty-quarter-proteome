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

## Regional and national soil genome catalogues (L1, ingested 2026-09-22)

A shelf that did not exist before this pass: genome-resolved soil surveys with a defined
geography and environmental metadata, as opposed to the pooled global family catalogues
(`pavlopoulos2023darkmatter`, `coelho2022gmgc`, `duan2024gmsc`, `aplakidou2025quadrupling`) and the
hot/cold-desert comparators already on the shelves. PDFs in
`resources/literature/regional-soil-catalogues/`. Every number below was read off the PDF or the
publisher record this session; page pointers are **printed journal pages** unless marked "PDF p.".

### Singleton et al. 2026 — Microflora Danica. *Nature* 649:971–981 · doi:10.1038/s41586-025-09794-2 (`singleton2026danica`)
On disk (`Singleton_2026_Microflora_Danica_atlas_Nature.pdf`). Online 3 Dec 2025, printed issue
22 Jan 2026; **Crossref indexes it as 2025**, the printed running head says 2026 — the bib entry
carries a note.
- **Design.** 10,683 shotgun metagenomes, mean 4.5 ± 2.1 Gb, 48.2 Tb total; sample types soil
  8,361, sediment 1,828, water 222, other 181 (p.971–972, Fig. 1a). Every sample carries GPS
  coordinates and a curated five-level habitat ontology (MFDO1–3; 28 MFDO1 categories) mapped to
  EMPO, Natura 2000 and EUNIS. Nationwide, one country, one campaign — **no temporal replication**:
  nothing in the paper resamples a site.
- **Genomes.** 19,253 MAGs of at least medium quality → **5,518 species at 95% ANI**, of which
  **4,604 (83.4%) novel vs GTDB R220** (p.977). CheckM2 1.0.2, GTDB-Tk 2.4.0, dRep `-sa 0.95 -nc 0.4`,
  SingleM-guided targeted recovery, Bakta 1.8.1.
- **Novelty is 16S-based, not proteome-based.** 605,861 ASVs → 141,252 species-level (98.7%) OTUs;
  **82.5% new species vs SILVA 138.1 but only 1.9% new families** (p.972, Fig. 2a). Eukaryotes: 77%
  new species vs PR2, 32% below 93% identity. The "novelty collapses above species" shape is exactly
  ours, reached by a different marker.
- **No protein catalogue, no gene catalogue, no unannotated/dark fraction anywhere.** UniRef100
  appears once in Methods, only to classify contigs during a targeted reassembly of two nitrifier
  type genomes. Function in the paper = one targeted guild (nitrifiers; GraftM amoA/nxrA models,
  SingleM, sylph, DRAM/KEGG r109).
- **Ecology.** Habitat separation of the prokaryotic community: **PERMANOVA R² = 0.27, P = 0.001;
  ANOSIM R = 0.69, n = 9,643** (p.974). Headline ecological result is an alpha-vs-gamma contrast:
  fields 12,797 common species vs grassland formations 26,721, i.e. human disturbance lowers gamma
  while alpha is unchanged ("hidden homogeneity"). pH invoked for bog/fen dispersion; random-forest
  habitat classification scored by PR-AUC.
- **Archaea.** Only Nitrososphaerota/Nitrososphaeraceae AOA are treated (new genera TA-21, TH5893,
  TH5896, TH1177; one undescribed species at median 4.6% of field samples, 320 MAGs). Halobacterota
  appears once as a tree label. **Thermoplasmatota is not mentioned in the main text or Methods.**
  This is main-text absence only — whether the 19,253 MAGs contain Thermoplasmatota is unchecked and
  would need Supplementary Note 7 / the released genome table.
- **Verdict.** Supports and sharpens §"The organisms are ordinary at every rank above order"; forces
  a qualification in Introduction ¶3 and Discussion ¶1. Danica is a large, non-pooled, fully
  metadata-linked catalogue, so "pooled catalogues strip context" cannot be stated unqualified —
  what Danica lacks is a proteome layer, not context.

### Ma et al. 2023 — SMAG. *Nat Commun* 14:7318 · doi:10.1038/s41467-023-43000-z (`ma2023smag`)
On disk (`Ma_2023_SMAG_soil_genomic_catalogue_NatCommun.pdf`). Author Correction (14:8079,
doi:10.1038/s41467-023-44072-7) is **typographical only** ("reservo)ir"; "SAMG"→"SMAG"); not filed.
- **Design.** 3,304 soil metagenomes, 363 in-house + 2,941 public, global, mainly grassland,
  cultivated land and forest (p.2, Fig. 1a–b). Single-sample assembly. **No environmental metadata
  analysis of any kind** — no pH, no climate, no drivers, no ordination. It is a resource-and-mining
  paper.
- **Genomes.** 40,039 MAGs at MIMAG medium-or-better; 3,641 (9.1%) high quality with rRNAs and ≥18
  tRNAs, a further 5,184 (13%) ≥90%/<5% without the rRNA complement (p.2). 95% ANI → **21,077 SGBs**;
  vs ~500k reference genomes at 95% ANI/30% AF, **16,530 (78.4%) are "unknown" uSGBs**; 70.8% of SGBs
  are singletons (p.3). Taxonomy GTDB **r202** (old), 88 bacterial + 11 archaeal phyla.
- **Higher-rank novelty, and this is the number that will be quoted against us.** By RED in GTDB
  r202: 6,392 unannotated genus-level, 1,166 family-level, 258 order-level, **31 class-level** and
  **2 candidate phylum-level** genome bins (p.3). Read carefully: these are *unannotated* bins from
  automatic RED, with no chimera control, no de novo tree adjudication and no release calibration —
  precisely the step our §1.1 performs and SMAG does not.
- **Dark fraction.** Reported only as a COG category: "function unknown" over-represented in the
  catalogue and in uSGBs in particular (Supp. Fig. 3b) — no denominator, no percentage, no protein
  clusters, no reference-absence call. **Pangenome layer:** 107 pangenomes over 2,200 SGBs with >10
  HQ MAGs; core genes better annotated than accessory genes against all five databases
  (Wilcoxon P < 0.001), and uSGB core genes less well annotated than kSGB core genes (eggNOG
  P = 0.054, KEGG P = 0.0005, GO P = 0.041) (p.5). **Selection layer:** 582,519,530 intraspecies SNVs
  from 2,448 SGBs; most phyla pN/pS < 1 (p.5–6).
- **Also released:** 43,169 BGCs and 8,545 CRISPR-Cas genes (abstract).
- **Archaea.** 11 archaeal phyla; only Thermoproteota (median genome 1.60 Mb) and Asgardarchaeota are
  discussed. **Thermoplasmatota is not named in the main text**, but the 11 archaeal phyla are never
  listed, so it is probably among them; the GTDB r202 strings in Supplementary Data 2 / the released
  catalogue (Zenodo record 8223844) settle it. Do not treat this as catalogue-level absence.
- **Verdict.** Pre-empts nothing but reframes three sections. (i) §1.1: the largest dedicated soil
  MAG catalogue reports 78.4% unknown species and 31 candidate classes; our "0 novel phyla or
  classes" needs the adjudication contrast stated explicitly or it reads as a weaker result.
  (ii) §2.1 / "the gap is in the reference set": SMAG is a **soil-specific** reference of 40,039
  genomes, and the dark-family search used MGnify90, GEM, SPIRE and five desert catalogues — a
  reviewer will ask why the largest soil-only catalogue is not among them. (iii) §1.2: SMAG's
  pN/pS < 1 and the core-vs-accessory annotation gap are the closest published precedents for our
  purifying-selection and annotation-gradient statements.

### Diamond et al. 2019 — Mediterranean grassland. *Nat Microbiol* 4:1356–1367 · doi:10.1038/s41564-019-0449-y (`diamond2019grassland`)
On disk (`Diamond_2019_grassland_soil_rainfall_depth_NatMicrobiol.pdf`). Banfield lab, Angelo Coastal
Range Reserve, California.
- **Design.** 60 metagenomes + 20 metaproteomes, three depths (10–20, 20–30, 30–40 cm), six sites,
  three of them under a **14-year spring-rainfall extension experiment** (p.1357). 1.2 Tb raw,
  67 Gbp assembled.
- **Genomes.** **793 near-complete genomes from 18 phyla**, accounting for about one third of
  organisms detected. rpS3 census: 10,158 sequences → 3,325 non-redundant species groups across 26
  phylum-level lineages (p.1357). No dereplication-to-species-representative framing, no novelty
  statistics, no protein-cluster or dark-matter layer at all.
- **Ecology, and this is the methodological precedent for our §"Encoded function shifts with soil
  position".** 391 genomes increased and 179 decreased with depth — **72% of assembled genomes are
  depth-stratified** (p.1361). Function then tested on the *differentially abundant genome sets*:
  C1 and CAZy diversity up near the surface, inorganic N turnover up at depth; rainfall extension
  flattens that partitioning. Statistics: DESeq2 with covariates (`plot_replicate + treatment + date
  + depth`), **MRPP in vegan (10,000 permutations)** for metadata influence, and custom 10,000-fold
  permutation tests for phylum/function enrichment against the 793-genome background (Reporting
  Summary, PDF p.15). Sampling **date is a modelled covariate**, so there is limited temporal structure,
  but no multi-year resampling design.
- **Archaea.** Euryarchaeota, Thaumarchaeota, Bathyarchaeota (2019 names, pre-GTDB). All archaeal
  lineages enriched at depth; all have very low CAZy diversity, handled by re-running the CAZy
  analysis without archaea. **No Thermoplasmatota** — but the taxonomy predates GTDB, so the name
  could not have appeared; Euryarchaeota as used here contains what GTDB later split out.
- **Verdict.** Supports §2's depth/position argument and gives it a citable precedent for the
  "differential abundance first, function enrichment second, permutation background" design. Also
  the cleanest precedent for our "lineage turnover, not within-lineage gene gain" reading, because
  Diamond's archaea-removed control is the same move as our within-phylum control.
- **There is no Diamond 2022 *Nat Microbiol* soil paper.** The 2022 Diamond-coauthored soil paper is
  Crits-Christoph, Diamond, Al-Shayeb, Valentin-Alvarado & Banfield, *ISME Commun* 2:70,
  doi:10.1038/s43705-022-00140-5 (46 MAGs, Acidobacteria BGC enrichment). The population-genetics
  companion on the same meadow is Crits-Christoph, Olm, Diamond, Bouma-Gregson & Banfield,
  *ISME J* 14:1834–1846 (2020), doi:10.1038/s41396-020-0655-x — 19 abundant species, nucleotide
  identity significantly associated with local geography for over half of them. Neither is filed;
  the 2020 paper is the better §1.2/§3 population-genetics comparator if that argument grows.

### Song et al. 2026 — continental-scale NEON soils + organic-matter chemistry. *Nat Commun* 17:5290 · doi:10.1038/s41467-026-71453-5 (`song2026continental`)
On disk (`Song_2026_NEON_continental_SOM_chemistry_NatCommun.pdf`). Read at abstract-plus-results
level, not deeply.
- 47 standardised US soil cores selected by respiration rate from 106 NEON soils; 0.76 Tbp of
  metagenome; **828 MAGs dereplicated to 358 genomes**, 48 high quality; paired FTICR-MS giving
  66,727 SOM molecules. Archaeal side is Nitrososphaerales only.
- The point of interest is the **design**, not the counts: this is the published precedent for
  pairing a genome-resolved soil survey with measured chemistry, which is what our
  §"Geochemistry structures known and unknown gene content beyond distance" does with XRF. Their
  genome yield is 47 cores to 358 dereplicated genomes.
- **Verdict.** Cite-only, in §3 and in the Introduction's "genomes plus measured chemistry" claim.
  Does not pre-empt: no dark-matter layer, no reference-absence call, chemistry is organic-matter
  composition rather than elemental.

### Li et al. 2023 — culturomics-based metagenomics of desert soils. *npj Biofilms Microbiomes* 9:67 · doi:10.1038/s41522-023-00439-8 (`li2023desertcbm`)
On disk (`Li_2023_desert_culturomics_metagenomics_npjBiofilms.pdf`). Small; cite-only.
- Two original Xinjiang desert soils (bulk + rhizosphere) plus 360 culture-enriched samples;
  300.14 Gbp shotgun; 580 bins → **148 MAGs (33 high, 115 medium quality)**, 146 bacterial and 2
  archaeal. 1,941 of 3,779 ASVs (51.4%) called potentially novel taxa, 1,095 of them culturable.
- Uses "microbial dark matter" in the taxonomic sense (unisolated taxa), not the sequence sense.
  **Term collision to keep straight**, like `harteveld2024darkfolds`.
- **Verdict.** Supports §1.2's cultured-genome control by precedent — culturing recovers taxa that
  direct sequencing misses in desert soil — and is the only desert paper found that pairs isolates
  with metagenomes. Too small to be a catalogue comparator.

### Leung, Bay, … Greening et al. 2025 — global desert survey (preprint, NOT on disk) · doi:10.1101/2025.11.29.691269 (`leung2025aridity`)
bioRxiv returned HTTP 429 to this node on four attempts; row added to `TO_DOWNLOAD.md`.
Metadata verified via Europe PMC (PPR1128095) and the bioRxiv API, abstract read in full.
- Standardised field survey of cold, hot and polar deserts on seven continents, with **measured**
  respiration, trace-gas oxidation and carbon-fixation activities alongside metagenomes. Desert
  communities dominated by streamlined Actinobacteriota and Chloroflexota; enriched in stress
  tolerance genes, mobile elements and antiviral systems; depleted in carbohydrate and protein
  degradation; enriched in chemosynthetic carbon fixation, atmospheric trace-gas and light energy
  harvesting. **Aridity is the primary driver** in their driver analysis.
- **Verdict — the most consequential item in this pass for §2.2.** It is a direct, activity-measured
  successor to `ji2017atmospheric` / `bay2021trace` covering seven continents, from the same group,
  and it will be published before or alongside us. Our §2.2 is capacity-only and states that
  explicitly; that framing survives, but the Discussion's trace-gas paragraph should cite it rather
  than be overtaken by it. Whether it reports MAGs, and whether any Arabian site is included, cannot
  be established without the PDF.

### Negative results of this search (one search round, not exhaustive: read as "not found")
- **No new dryland MAG catalogue since 2023** beyond what is already on the shelves: nothing for
  Atacama, Namib, Negev, Sahara, Gobi or Tengger surfaced past `liu2025qaidam` and
  `andreanigerard2024bgc`. The live dryland item is the Leung preprint above.
- **No Australian soil MAG catalogue exists.** BASE (Bissett et al., *GigaScience* 5:21, 2016) is a
  16S/ITS amplicon database with soil chemistry, not genomes; the Australian Microbiome publication
  list carries no continental MAG set.
- **No Chinese national soil MAG catalogue distinct from SMAG** was found; the Tibetan Plateau
  catalogue (*Nat Commun* 15, 2024, doi:10.1038/s41467-024-45895-8) is aquatic, not soil.
- **No NEON MAG catalogue as such** — Song 2026 is the largest NEON genome-resolved study and is 358
  dereplicated genomes.
- No European soil MAG catalogue beyond Danica; LUCAS has no genome-resolved counterpart.
- Cross-ecosystem MAG repositories found but not worth PDFs: `gcMeta 2025` (*NAR* 2026,
  doi:10.1093/nar/gkaf1115) and `MAGdb` (*Genome Biol* 26, 2025, doi:10.1186/s13059-025-03711-6).

### What this changes for main.tex
1. **Introduction ¶3 and Discussion ¶1 overstate the case against prior work.** Both say pooled
   catalogues have "no repeated sampling, no measured chemistry and no cultured neighbour". Danica
   (`singleton2026danica`) has metadata, geography and a habitat ontology; Diamond
   (`diamond2019grassland`) has depth, a 14-year rainfall manipulation and modelled sampling date;
   Song (`song2026continental`) has measured organic-matter chemistry. The defensible line is
   narrower and stronger: regional catalogues restore ecology but stop at the genome — **none of
   them carries a proteome layer with a matched control**, and none reports a reference-absence
   fraction at all. Recommend rewriting the three-gap sentence as a two-part claim: the global
   catalogues have the proteome and no context; the regional catalogues have the context and no
   proteome.
2. **Introduction ¶1 needs a "what a regional genome catalogue looks like" sentence.** As written it
   jumps from drylands to "the Rub' al-Khali had no genome-resolved microbial survey", which invites
   the reader to supply the wrong comparison. One sentence naming Danica (10,683 samples → 5,518
   species) and SMAG (3,304 metagenomes → 21,077 SGBs) sets the scale the reader should hold, and
   makes 5,438 species from 551 samples of one sand sea read as dense rather than small.
3. **§1.1's "0 novel phyla or classes" must be stated against SMAG's 31 candidate classes and 2
   candidate phyla** (`ma2023smag`), otherwise it reads as a null result rather than as the product
   of an adjudication step SMAG did not perform. The contrast is method, not biology: automatic RED
   in GTDB r202 versus GUNC chimera screening plus release-calibrated de novo trees in R232.
4. **§"The gap is in the reference set" has a named omission, and it is fixable.** SMAG is 40,039
   soil genomes and is not in the 1.258-billion-sequence search, even though the argument is that
   biome-matched references beat depth and SMAG is the largest biome-matched reference that exists.
   The catalogue is released (Zenodo record 8223844, surfaced in this search but **not verified by
   download**), so this is one DIAMOND pass rather than a caveat to write around. This is the single
   most likely reviewer objection arising from this pass.
5. **The §2 depth/position analysis has a precedent to cite** (`diamond2019grassland`): differential
   abundance first, function enrichment on the differential sets second, permutation background,
   and an archaea-removed control that mirrors our within-phylum control.
6. **Thermoplasmatota is absent from the main text of all three, which is suggestive but not yet
   evidence.** SMAG discusses only Thermoproteota and Asgardarchaeota among its 11 archaeal phyla,
   Danica treats only Nitrososphaerota, and Diamond predates the GTDB name. **None of this has been
   checked against the released genome tables**, and Thermoplasmatota is common enough in soil that
   it is probably present in both catalogues. The claim §1.4 could make — that these classes are
   poorly sampled rather than merely unnamed — needs the GTDB strings from SMAG (Zenodo 8223844 /
   Supplementary Data 2) and Danica (Supplementary Note 7). Hand-off to the Thermoplasmatota brief;
   do not write it up from main-text silence.

## Thermoplasmatota brief (L2, ingested 2026-09-22)

Read in full this session: `maeke2025thermoplasmatota` (already on disk), plus four papers downloaded
into `resources/literature/thermoplasmatota/` — `sheridan2022lutacidiplasmatales`,
`zhou2019thermoprofundales`, `diamond2022cummo`, `finstad2017atacama`. Databases queried live on
**2026-09-22**: GTDB **R232** (`gtdb-api.ecogenomic.org`), NCBI Datasets + BioSample, MGnify genomes
API, SMAG Supplementary Data 3 (`ma2023smag`, Springer MOESM5). Every number below was pulled this
session unless marked UNVERIFIED.

### Q1 — what the two classes are, and four historical names that do NOT belong to them

**Correction to the premise of the request.** Most of the historical names in the brief map to a
*different* class. Checked directly against GTDB R232:

| historical name | GTDB R232 placement (verified) |
|---|---|
| TMEG (Terrestrial Miscellaneous Euryarchaeota Group) = Ca. Lutacidiplasmatales | `c__Thermoplasmata; o__Lutacidiplasmatales`, 186 species reps |
| MBG-D = Thermoprofundales (Zhou 2019) = Ca. Proteinoplasmatales | `c__Thermoplasmata; o__SG8-5`, 22 reps. GTDB carries no `Thermoprofundales`/`Thermoprofundaceae`; `o__SG8-5` is anchored by *Candidatus Proteinoplasmatales archaeon SG8-5*. The MBG-D = SG8-5 equation is a name match, **UNVERIFIED** by phylogeny |
| "Group E2" of 16S clone libraries (Iino 2013, Ca. *Methanogranum caenicola*, strain Kjm51a, digester sludge) | became **Methanomassiliicoccales**, `c__Thermoplasmata`, 784 genomes. **This is not GTDB `c__E2`**. The GTDB order count is verified; the Group E2 to Methanomassiliicoccales identification is from the Iino 2013 abstract, **paper not read**, no bib entry made |
| Ca. Sysuiplasmatales | `c__Thermoplasmata; o__Sysuiplasmatales`, 10 reps |
| Ca. Angelarchaeales (Diamond 2022, CuMMO-encoding) | formerly RBG-16-68-12; `c__Thermoplasmata; o__RBG-16-68-12`, 111 genomes |
| Rice Cluster III | within Thermoplasmata; exact GTDB order **UNVERIFIED** |

The entire soil-Thermoplasmata literature a referee will reach for (TMEG, MBG-D, Sysuiplasmatales,
Angelarchaeales) is **class Thermoplasmata**. It is neighbour evidence for R4, never a synonym.
Do not write "formerly TMEG" or "formerly MBG-D" for these two classes.

**What they positively are.**
- `c__SW-10-69-26` is named for **GCA_003009755.1, "Thermoplasmatales archaeon SW_10_69_26"**, a
  Banfield-lab MAG from **BioProject PRJNA351262 = 26 halite salt crusts along a fog gradient at
  Salar Llamara, hyperarid Atacama Desert, Chile** (BioSample SAMN06020844: `env_broad_scale =
  desert`, `env_local_scale = salt flat (salar)`, `env_medium = salt crust`, 21.21 S 69.84 W,
  elev 836 m, surface). That is `finstad2017atacama`. **The type genome of this class is already a
  hyperarid-desert genome.** Finstad's own text assigns its single Thermoplasmatales MAG to "Marine
  Group II"; that the sentence refers to this genome is inference (same BioProject, same sample
  naming), so do not quote it.
- `c__E2` is a GTDB alphanumeric placeholder of undocumented provenance (GTDB publishes no note).
  Its largest order is `o__DHVEG-1` (118 of 197 reps), Deep-sea Hydrothermal Vent Euryarchaeotic
  Group 1 of Takai & Horikoshi 1999; subgroup DHVE2 is *Aciduliprofundum boonei*, a thermoacidophile.
  Other orders: `o__UBA202` 39, `o__JACPAO01` 38, `o__GCA-002838935` 1, `o__UBA9212` 1.
- **Neither class has a cultured or type-material genome**: 0 of 21 and 0 of 197 species reps carry
  `isNcbiTypeMaterial` in GTDB R232. R4's "no cultured genome reaches these classes" is verified.
- **Tree position (GTDB R232 `ar53.tree`, 10,122 archaeal reps, parsed this session).** Within
  Thermoplasmatota the branching order is B47-G6, UBA186, **E2**, **SW-10-69-26**, JASEFT01,
  DAOVMN01, Poseidoniia, JAKEGO01, **Penumbrarchaeia**, DTKX01, Thermoplasmata. `c__E2` is sister to
  everything after it (support 35/39), `c__SW-10-69-26` sister to everything after *that* (61): a
  grade of two early-branching lineages, not a clade. `c__Penumbrarchaeia` sits on the far side of
  Poseidoniia, sister to `c__JAKEGO01`. **Backbone support is 35–61, so all of this is weak.**
- Environments of the GTDB reps (NCBI BioSample, mined this session). **SW-10-69-26 (21):** 16
  soil/terrestrial, 4 freshwater/engineered/saline, 1 coastal; USA 9, China 3, Chile 2, South Africa 1,
  Russia 1, Canada 1. **E2 (197):** 46 hydrothermal/hot spring (27 of them Guaymas/Pescadero, Mexico),
  41 freshwater/engineered/saline, 33 "terrestrial" unspecified, 24 soil, 23 coastal/tidal, 18 marine.
  So SW-10-69-26 is already a terrestrial class in the public record; E2 is a mixed one with a large
  hydrothermal component.
- Largest single donor of reference genomes to both classes is **PRJEB82487, "Bin Chicken: targeted
  coassembly of metagenomes"** (Aroney/Woodcroft, registered 2024-12-12): 63 of 197 E2 reps and 5 of
  21 SW-10-69-26 reps. A novelty-mining exercise, not an ecological survey.

**Reported metabolism, for the neighbours (all class Thermoplasmata).**
- `zhou2019thermoprofundales` (ISME J 13:885–901): MBG-D/Thermoprofundales, 146 genomes, five
  subgroups A–E with distinct habitats; peptide transport and assimilation, acetate and ethanol by
  fermentation, secreted extracellular peptidases (M09B collagenase highly expressed in mangrove
  metatranscriptomes), cosmopolitan in global sediment. This is the closest published analogue of
  R4's "aerobic scavenger of peptides", but it is **anaerobic and sedimentary**.
- `sheridan2022lutacidiplasmatales` (Nat Commun 13:4110): 35 TMEG genomes from **9 surface (0–15 cm)
  and 6 subsurface (30–60 cm) acidic Scottish soils**, 22 surface + 13 subsurface genomes; average
  completeness 78% (45–98%), contamination 3%. Complex IV in 19 of 22 higher-completeness genomes =
  "a common aerobic metabolism"; obligate heterotrophs, no complete carbon-fixation pathway; sulfite
  oxidase rather than reductase; some facultatively anaerobic (acetate/ethanol). **Lutacidiplasmatales
  are predicted non-motile — they lack archaellum genes present in Thermoplasmatales and
  Poseidoniales.** Aerobic respiration, acid tolerance and glycolysis are **convergent**, acquired
  independently by divergent Thermoplasmatota lineages. Their coxA tree includes
  *Thermoplasmatota archaeon SW-10-69-26* as a single-genome branch: **our class already has a
  published cytochrome c oxidase**.
- `diamond2022cummo` (ISME J 16:1348–1362): Ca. Angelarchaeales (= RBG-16-68-12) in soils and
  sediments encode copper membrane monooxygenases; 20 of the genome set carry pmo/amo clusters.
  Aerobic, soil/sediment Thermoplasmata that may oxidise ammonia. Relevant as a **check on R4's
  "no ammonia oxidation"**: the claim is about our two classes, but a referee knows archaeal amo
  exists outside Nitrososphaerales in this phylum.

### Q2 — soil and desert archaea surveys: do they see Thermoplasmatota?

- The received story is exactly "Nitrososphaerota dominate soil archaea". `bates2011soilarchaea`
  (ISME J 5:908–917, 146 globally distributed soils, barcoded pyrosequencing): soil archaea are
  **dominated by just two phylotypes from a constrained clade within Crenarchaeota (now
  Nitrososphaeria), together >70% of all archaeal sequences**; archaea average **2%** of the community
  (0 to >10%); soil C:N is the only consistent correlate.
- Hyperarid soils repeat it. `hwang2021nitrosodeserticola` (Microbiome 9:234, PDF now in
  `dryland-ecology/`, abstract and conclusions read, body not): 11 metagenomes from three Atacama
  hyperarid-core boulder fields; **archaea are significantly more abundant in soil below boulders
  than beside them at all sites across up to 205 km**, and the taxon responsible is a new
  ammonia-oxidising Thaumarchaeota genus, Ca. *Nitrosodeserticola* (8 mid-to-high-quality genomes;
  amo, 3-hydroxypropionate/4-hydroxybutyrate carbon fixation, acetate use, oxidative-stress and
  DNA-repair genes). The closest published "sheltered-microsite archaeal enrichment in a hyperarid
  desert" result, and it is Nitrososphaerota, not Thermoplasmatota. Useful as the contrast R4's
  depth result is measured against.
- **No dryland, hot-desert, Arabian, Saharan, Atacama, Namib or Australian-arid soil survey found
  this session reports Thermoplasmatota as the dominant or a subsurface-enriched archaeal group.**
  Searches returned Thermoplasmatota in soil only as a minor phylum, or in hydrothermal, geothermal
  and deep-sediment settings. That is a negative result from ~6 targeted searches, not an exhaustive
  review; it is enough to say the pattern is unreported, not enough to say it has been looked for.
- **Amplicon blind spot, for the companion-survey sentence.** `eloefadrosh2016gaps` (Nat Microbiol
  1:15032): primer fidelity checked against SSU genes recovered from **>6,000 assembled metagenomes**;
  **~10% of environmental sequences are missed by classical PCR SSU surveys, mostly CPR bacteria and
  uncharacterised Archaea**. This is the clean citation for "amplicons cannot reach the archaea"
  (main.tex line 56). The companion survey is V3–V4; in-silico primer evaluations report very poor
  archaeal coverage for 341F/805R, but the one retrieved this session scores against an *oral*
  archaeal database, so it is not usable here. Do not cite a V3–V4 archaeal-coverage percentage
  without a soil-relevant source.

### Q3 — precedents for auxotrophy, cross-feeding, CRISPR and archaella

- **Peptide/protein scavenging by uncultured archaea is well precedented** (`zhou2019thermoprofundales`;
  `maeke2025thermoplasmatota` p.8: amino-acid degradation in all 35 MAGs, extracellular peptidases
  C11A/M14B/S08A in sediment families, oligopeptide transporters and aminopeptidases in all families,
  deamination via aspB/alaA, 2-oxoacids to acyl-CoA via por/ior/kor/vor, ATP via ADP-forming acetyl-CoA
  synthetase). R4's protease/BCAA-transporter profile is **expected** for this phylum.
- **What is not precedented in what was read: an explicit branched-chain amino-acid auxotrophy call
  for a Thermoplasmata lineage.** Maeke reports *ilvE* (BCAA aminotransferase) in single MAGs but
  makes no auxotrophy claim; Sheridan calls Lutacidiplasmatales obligate heterotrophs but does not
  itemise amino-acid auxotrophy. The nearest framed precedents are DPANN gene loss (loss of amino-acid
  and cofactor biosynthesis is common in symbionts relative to free-living siblings) and Bathyarchaeia
  as a cobalamin *donor* rather than auxotroph. **Both abstract-level only, neither read.** R4's
  "keeps the transaminase, loses the carbon skeleton" is a sharper statement than anything found in
  this literature, which is an opportunity and an exposure: state it as encoded-pathway absence, not
  as demonstrated auxotrophy.
- **Archaellum is a real contrast, not a shared trait.** Sheridan: Lutacidiplasmatales are non-motile
  and lack archaellum genes that Thermoplasmatales and Poseidoniales have. R4's "all 237 carry an
  archaellum" therefore separates these classes from the best-studied soil Thermoplasmata order.
- **No CRISPR-rich Thermoplasmata precedent was found.** Maeke does not report CRISPR-Cas; Sheridan
  does not; Diamond does not. R4's 79% CRISPR-Cas has no comparator in this literature, so it needs
  a within-catalogue baseline (what fraction of the other 5,201 reps carry a system?) rather than a
  citation.

### Q4 — Maeke 2025's rarity argument, stated precisely (for contrast with our matched null)

PDF on disk, page numbers are the printed article pages.
1. **The question is posed on p.10**, section "Rare MAGs in the phylum Thermoplasmatota hold higher
   numbers of genes encoding functionally unknown genes": *"Considering the low relative abundance of
   Ca. Penumbrarchaeia MAGs in the screened metagenomic runs, this raises the question of whether the
   unexpectedly high number of genes encoding hypothetical and unknown proteins in this novel class
   may be related to its rarity."*
2. **The test (p.10, Fig. 6a–b).** Rarity index = median relative abundance weighted by fraction of
   occurrence across screened runs; MAGs below the median are "rare", above "common". *"The percentages
   of genes encoding unknown proteins in the rare group were higher than in the common group."*
   73% of their Penumbrarchaeia MAGs are rare. Fig. 6c gives unknown % per Thermoplasmatota order.
3. **The mechanism they name is database composition, not genome quality (p.13).** *"It is a common
   feature of uncultivated and rare taxa to contain high numbers of unknown protein families with up
   to 70% of genes that cannot be annotated by databases such as KEGG and NR, as these databases use
   complete genomes and functionally characterized genes as reference."* The only genome-quality
   caveat in the paper (p.12) is about **missing metabolic genes** in water-column families
   ("missing genes in some families could also result from the incompleteness of the analyzed MAGs"),
   not about the unknown fraction.
   **→ main.tex line 423 currently attributes "rarity and genome quality" to this paper. "Genome
   quality" is not their argument. Reword to rarity and reference-database composition, or the
   citation is a misattribution a referee who knows the paper will catch.**
4. **Their definitions differ from ours.** Maeke: *hypothetical* = at least one hypothetical hit in
   KEGG or NR; *unknown* = **no database match at all** (NR, not UniRef50). Ours: no UniRef50
   homolog. NR is the larger and more MAG-saturated target, so their "unknown" is systematically
   smaller than our "dark". Their numbers: unknown 1.6–63% per MAG, **average 26%**; AGNOSTOS genomic
   unknowns 14–29%, environmental unknowns 0.24–17%.
5. **Their rarity axis is aquatic only, and they say so (p.10).** MAGs in the "not detected" group
   *"originated from soil, biodigesters, or human- and animal-associated habitats … low or missing
   relative abundances for these nonaquatic Thermoplasmatota do not necessarily classify these taxa
   as overall rare but as rare in the screened environments."* Their screen was 8,287 aquatic runs;
   a soil lineage is unmeasurable in it by construction.
6. **Scale.** 35 MAGs, completeness 80.3–97.9%, contamination 0.05–4.35%, completeness-adjusted
   genome sizes 1.14–3.94 Mb; detected in 128 of 8,573 runs (1.49%); relative abundance 0–0.32%,
   mean 0.031%. Our classes are detected in 481 of 549 samples at a median 0.39% of reads.
7. **The number that matters most — Maeke's own Table S16 (405 public Thermoplasmatota genomes;
   downloaded and recomputed this session).** Median % genes with no database match, by class:

   | class | n | median | mean | max |
   |---|---:|---:|---:|---:|
   | UBA186 | 2 | 33.2 | 33.2 | 35.7 |
   | EX4484-6 (Penumbrarchaeia) | 22 | 20.5 | 26.7 | 61.5 |
   | SW-10-69-26 | 2 | 5.4 | 5.4 | 7.2 |
   | Thermoplasmata | 266 | 5.2 | 9.1 | 59.7 |
   | E2 | 36 | 4.0 | 5.9 | 41.5 |
   | Poseidoniia | 77 | 1.6 | 2.4 | 19.0 |

   By orders: DHVEG-1 4.1 (n=21), UBA202 3.3 (n=12), UBA9212 6.8 (n=3), SW-10-69-26 7.2 (n=1),
   JACQPN01 3.6 (n=1). Habitats of their E2 genomes: 17 "extreme", 4 groundwater, 4 marine sediment,
   **2 soil**; their two SW-10-69-26 genomes are "extreme" and groundwater.
8. **That table looks like a counter-argument and is not, because deposition status confounds it.**
   Maeke's "unknown" means no hit in NR, and NR contains the deposited proteins of the very genomes
   being scored. Splitting S16 by whether the genome identifier is a GenBank accession or a
   study-internal bin name (recomputed this session):

   | S16 genomes | n | median % unknown | mean |
   |---|---:|---:|---:|
   | deposited (GCA/GCF, therefore in NR) | 383 | 3.9 | 7.4 |
   | newly binned by Maeke, not in NR | 22 | 22.5 | 29.6 |

   Within Penumbrarchaeia alone: 4 deposited genomes 2.8% median, 18 new ones 25.5%. Within
   Thermoplasmata: 262 deposited 5.1%, 4 new 14.5%. **All 36 E2 and both SW-10-69-26 genomes in S16
   are deposited.** Their 4% is a measurement of being in GenBank, not of E2 darkness, and the same
   confound sits underneath their rarity-to-unknown correlation, since rare lineages are the ones
   least likely to have deposited genomes. The comparison that is actually like-for-like is their
   **new** MAGs (mean 26–30%) against our **new** MAGs (median 31%), and those agree.
   **→ R4's matched null is what removes this confound**: every comparator is an undeposited genome
   from the same soils, the same binning and the same annotation run, so deposition status is held
   constant by construction. One sentence in the Discussion stating the definitional difference
   (no NR match versus no UniRef50 homolog) and this confound converts the apparent contradiction
   into support. NR being the larger target is a second, weaker reason their numbers run low;
   the self-inclusion mechanism is the directional one.
9. **"Sister class" is wrong.** From the R232 tree above, Penumbrarchaeia is not adjacent to either
   class; Poseidoniia and two other classes lie between. Use "another Thermoplasmatota class" or
   "a different class in the same phylum" at main.tex line 423.

### Q5 — is a desert-dominant set of 237 from these two classes unreported? Counts per catalogue

All rows queried 2026-09-22. Units differ and are stated; do not merge them.

| catalogue | version / basis | c__SW-10-69-26 | c__E2 | source |
|---|---|---:|---:|---|
| **This catalogue** | 5,438 ANI95 species reps | **185** | **51** | R4 |
| GTDB | R232, species representatives | 21 | 197 | `gtdb-api.ecogenomic.org/search/gtdb`, `searchField=gtdb_tax`, `gtdbSpeciesRepOnly=true` |
| GTDB | R232, all genomes | 21 | 248 | `gtdb-api.ecogenomic.org/taxon/p__Thermoplasmatota` (`total` field) |
| GTDB | R232, type-material genomes | 0 | 0 | same query, `isNcbiTypeMaterial` |
| MGnify genomes | all catalogues, species reps | 10 (8 soil-v1-0, 2 marine-sediment-v1-0) | 73 (58 marine-sediment, 11 soil-v1-0, 4 marine-v2-0) | `ebi.ac.uk/metagenomics/api/v1/genomes?taxon_lineage=` |
| SMAG (`ma2023smag`) | 40,039 soil MAGs (bins, not species reps) | 14 | 35 | Supplementary Data 3, Springer MOESM5, recounted from the GTDB string column. **Floor**: 1,351 SMAG bins are Thermoplasmatota and 372 of them carry no class assignment, so some could belong to either class |
| Maeke's public screen (`maeke2025thermoplasmatota`) | 405 Thermoplasmatota genomes from 11,479 assemblies | 2 | 36 | Table S16 |
| GEM | — | in GenBank, therefore inside the GTDB counts | | not queried separately |
| SPIRE | — | **not queryable**: no public API or taxon route from this node | | `spire.embl.de` returns the SPA shell |
| Microflora Danica (`singleton2026danica`) | — | **no genome-level GTDB table** in the published supplementary (Supplementary Data 1–6 are core genera, AOA/nitrifier annotation, sampling and sample metadata) | | Springer MOESM3–8 downloaded and inspected |
| Sandpiper (SRA-wide singlem survey) | — | **not queryable**: JS app behind reCAPTCHA, no JSON route found | | `sandpiper.qut.edu.au` |

Context for those counts: GTDB R232 Thermoplasmatota is 4,961 genomes in 12 classes; Thermoplasmata
2,643 and Poseidoniia 1,953 hold 93% of them. **SMAG, 3,304 soil metagenomes and 40,039 MAGs, yields
14 + 35 = 49 bins from these two classes; this catalogue yields 236 species representatives from 551
samples of one sand sea.** That comparison is the defensible form of "unreported at this scale".
The **Microflora Danica and Sandpiper rows are open items** — both would settle the soil-versus-desert
question better than anything above, and neither is reachable from this node.

### 2024–2026 Thermoplasmatota work a referee would expect cited
- `maeke2025thermoplasmatota` — already cited.
- `sheridan2022lutacidiplasmatales`, `zhou2019thermoprofundales`, `diamond2022cummo` — the three
  standard soil/sediment Thermoplasmatota genomics papers; **all three now on disk**.
- Ca. Sysuiplasmatales, AEM 2021, doi:10.1128/aem.01065-21 — acid-mine-drainage Thermoplasmata order.
  Cite-only if at all; **not read, not on disk**.
- "Genomic Evidence for the Recycling of Complex Organic Carbon by Novel Thermoplasmatota Clades in
  Deep-Sea Sediments", mSystems 2022 — cite-only. **Not read**, no bib entry made.
- "Recovery of nearly 3,000 archaeal genomes from 152 terrestrial geothermal spring metagenomes",
  Sci Data 2025, doi:10.1038/s41597-025-04493-z (Tengchong; Thermoplasmatota among the most abundant
  archaeal phyla). Plausible source of the three Tengchong hot-spring E2 reps. Cite-only; **not read**.
  No bib entries were made for the last three; make them only if the text actually cites them.

### What this changes for R4 and the Discussion

1. **Add the strongest new fact to R4.** The type genome of `c__SW-10-69-26` is itself a
   hyperarid-desert genome, from Atacama halite salt crust (`finstad2017atacama`). One clause turns
   "a class nobody has sampled" into "a class first seen in another hyperarid desert and never since
   sampled at scale". Sixteen of its 21 GTDB reps are terrestrial.
2. **Fix the Maeke attribution (main.tex line 423).** "rarity and genome quality" should be rarity
   and reference-database composition. Genome quality is not Maeke's argument for darkness.
3. **Drop "sister" (main.tex line 423).** Penumbrarchaeia is not sister to either class in the R232
   archaeal tree; two classes and Poseidoniia lie between. "Another Thermoplasmatota class" is safe.
4. **Use Table S16 as support, with one sentence of method.** It looks like a contradiction (their
   public E2 genomes score 4.0% unknown against our 31% dark) and is not: their metric is no NR hit,
   NR contains the deposited genomes being scored, and across all 405 of their genomes deposited ones
   sit at 3.9% and newly binned ones at 22.5%. All of their E2 and SW-10-69-26 genomes are deposited
   and none of ours are. Their new MAGs (26–30%) and ours (31%) agree. The Discussion sentence that
   already says the matched null removes their explanation should name what it removes: reference
   composition and deposition status, held constant because every comparator is an undeposited genome
   from the same soils and the same annotation run.
5. **Two claims now have published precedent and should cite it.** The aerobic reading
   (cytochrome c oxidase; `sheridan2022lutacidiplasmatales` already places SW-10-69-26 in its coxA
   tree, and shows aerobic respiration is convergent in this phylum) and the peptide-scavenger
   reading (`zhou2019thermoprofundales`). Neither is novel for Thermoplasmatota; what is new is the
   combination with an archaellum, aerobic soil habit and subsurface enrichment.
6. **Two claims are genuinely without a comparator and should be phrased carefully.** The archaellum
   in all 237 is a contrast with Lutacidiplasmatales, which are non-motile; 79% CRISPR-Cas has no
   published Thermoplasmata baseline, so it needs a within-catalogue denominator rather than a
   citation. BCAA auxotrophy should stay an encoded-pathway absence, not a demonstrated auxotrophy.
7. **Concede the amo question in one clause if a reviewer would.** `diamond2022cummo` shows
   soil/sediment Thermoplasmata (Ca. Angelarchaeales = RBG-16-68-12) that encode CuMMOs. R4's
   "no ammonia oxidation" is a within-class statement and survives, but the phylum-level claim does not.
8. **Never write "formerly TMEG" or "formerly MBG-D".** Those are `c__Thermoplasmata` orders. Also
   do not connect `c__E2` to the 16S "Group E2" of Iino 2013, which became Methanomassiliicoccales.
9. **The Q5 table is citable as a count comparison**, with SMAG as the load-bearing row: 49 bins from
   40,039 soil MAGs against 236 species representatives here.
