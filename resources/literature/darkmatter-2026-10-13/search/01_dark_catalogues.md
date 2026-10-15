# Search Synthesis & Literature Expansion: 01_dark_catalogues

**Direction:** Unknown microbial protein families, AGNOSTOS, FESNov, GMGC, ecosystem specificity, novelty of longitudinal biome sampling versus pooled catalogues.  
**Target Manuscript:** P4 Dark-Proteome Manuscript (`manuscript/main.tex`, `notes/DARKFAM_PRIMER.md`, `notes/dossiers/02_dark_definition_families.md`, `notes/dossiers/04_reference_gap.md`).  
**Date:** 2026-10-13  
**Output Files:**
- Manifest: [`manifest.tsv`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/manifest.tsv)
- Papers Folder: [`resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/)
- Synthesis: [`01_dark_catalogues.md`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/search/01_dark_catalogues.md)

---

## 1. Executive Summary & Core Conceptual Tension

The central thesis of the P4 manuscript is that the microbial "dark proteome" (proteins lacking reference homologs in databases such as UniRef50) is not primarily an exotic hallmark restricted to candidate phyla, nor an artifact of assembly fragmentation, but rather represents the vast, uncataloged **strain-variable accessory genome** of dominant, ordinary microbial lineages. 

In the Rub' al-Khali sand sea survey (551 soil samples, 60 sites, 5 expeditions over 31 months, 5,438 species-level MAGs), 106,042 dark families (holding 536,437 member genes) were identified under strict admission criteria ($\ge 3$ diverged members, $\ge 3$ samples, and $\ge 2$ sites). Four independent lines of evidence validate that these families represent bona fide coding biology:
1. **Recurrence:** The median family is assembled in 12 samples across 9 sites.
2. **Longitudinal Persistence:** Read-mapping recruitment detects 79.4% of families in all five expeditions spanning 31 months (28.9% independently assembled in all five trips).
3. **Purifying Selection:** Median pairwise $dN/dS = 0.425$ (with 96.7% $< 1$), closely anchored by an identically passed, length- and size-matched UniRef50-matched control ($0.388$), while a +1 nucleotide frameshift control yields $dN/dS = 1.087$.
4. **Lineage Fidelity & Strain-Level Decoupling:** 97.9% of attributable families sit within a single GTDB phylum (dominated by ordinary Actinomycetota at 32% and Thermoplasmatota at 19%). Crucially, across 23 species cultured from the same basin, isolates and MAGs share 80% of their annotated proteins but only 16% of their dark proteins.

### The Literature Gap: Pooled Planetary Catalogues vs. Longitudinal Biome Sampling
The dominant paradigm in microbial dark-matter discovery is built on massive **pooled cross-sectional catalogues** (e.g., AGNOSTOS [Vanni et al. 2022], GMGC [Coelho et al. 2022], Pavlopoulos et al. 2023, FESNov [Rodríguez del Río et al. 2024], SPIRE [Schmidt et al. 2024], and Aplakidou et al. 2026). While these catalogues have performed the invaluable service of documenting that tens of thousands to millions of novel protein families exist globally, their pooled design introduces five fundamental structural limitations that P4 directly addresses and challenges:

1. **Stripping Spatio-Temporal and Environmental Context:** Pooled meta-catalogues aggregate samples across thousands of unrelated projects. Consequently, they cannot test whether an unannotated gene cluster persists in a resident population over time or represents transient moribund/relic DNA, contaminant spillover, or single-sample assembly chimeric noise.
2. **Conflating Abundance/Rarity with Biological Reality:** In pooled catalogues (like GMGC and Pavlopoulos et al.), rarity in global cross-sectional datasets is conflated with insignificance or artifact risk (hence Pavlopoulos et al. and NMPFamsDB require $\ge 100$ global members to define a family, systematically discarding localized or niche-specific protein families). P4 demonstrates that families with modest member counts are highly persistent (79.4% across 5 expeditions over 31 months) when measured in a dedicated longitudinal transect.
3. **Misattributing Darkness Exclusively to Uncultivated "Exotic" Phyla:** Studies like FESNov (Rodríguez del Río et al. 2024) frame unknown gene families (UGFs) as synapomorphies of deeply branching, uncultivated taxa (e.g., CPR/Patescibacteria). P4 challenges this narrative by demonstrating that in hyperarid soil, the dark proteome is overwhelmingly carried by dominant, well-known phyla (Actinomycetota, Thermoplasmatota), where darkness is a property of the open, strain-variable accessory genome that isolate genomes never capture.
4. **Coarse Ecological Categorization:** Multi-biome catalogues (like Aplakidou et al. 2026) report that 63.4% of novel families are restricted to a single ecosystem, but their ecological definitions rely on high-level GOLD categories (e.g., "Terrestrial: Soil"). This aggregates hyperarid deserts, temperate forests, and agricultural lands into one bucket, completely missing the edaphic and geochemical drivers of localized functional dark matter.
5. **Absence of Calibrated Controls for Novelty Calling:** Cross-sectional catalogues report "no-hit against database" as novel biological discovery without passing a matched set of known genes through the identical clustering, alignment, and threshold steps to measure false-negative reference dropouts and catalogue-size E-value sensitivity.

---

## 2. Search Methodology & Execution Log

Searches were conducted across PubMed, PubMed Central (PMC), Europe PMC, bioRxiv, Google Scholar, and publisher platforms (Nature, Oxford Academic/NAR, ASM).

### Search Queries Executed
1. `"AGNOSTOS" "unifying the known and unknown" Vanni eLife`
2. `"FESNov" "Functional and evolutionary significance of unknown genes" OR "AGNOSTOS" 2024 OR 2025 OR 2026`
3. `"Functional and evolutionary significance of unknown genes from uncultivated taxa" OR "FESNov" 2024 OR 2025 OR 2026 microbial protein families`
4. `"NMPFamsDB" OR "Novel Metagenome Protein Families Database"`
5. `"NMPFamsDB" "Baltoumas" OR "Pavlopoulos" Nucleic Acids Research OR NAR`
6. `"GMGC" OR "Global Microbial Gene Catalog" 2024 OR 2025 OR 2026 "Coelho" OR "Bork"`
7. `"dark proteome" OR "unknown protein families" OR "environmental unknowns" metagenome 2024 OR 2025 OR 2026`
8. `"AGNOSTOS" "environmental unknowns" "genomic unknowns"`
9. `"unknown protein families" metagenomics 2024 OR 2025 OR 2026`
10. `"Quadrupling the protein family space with global metagenomics" "Aplakidou" Nucleic Acids Research`
11. `"Quadrupling the protein family space with global metagenomics" "doi.org"`
12. `"nf-core/proteinfamilies" OR "Baltoumas" "proteinfamilies"`
13. `"WasteFams" OR "metagRoot" OR "ArthroVerse" "protein families"`
14. `"EnvoFams" OR "metagRoot" "Pavlopoulos" doi`
15. `"WasteFams" "Pavlopoulos" OR "Chasapi" OR "Lamari" doi OR biorxiv`
16. `longitudinal soil metagenomics "temporal stability" OR "persistence" "gene" OR "accessory"`
17. `"soil" "pangenome" "accessory genome" "strain" "metagenome" 2024 OR 2025 OR 2026`
18. `Choudoir "Comparative pangenomics" mSphere doi`
19. `"metagenomic catalog" OR "gene catalog" "cross-sectional" OR "pooled" limitations OR bias "longitudinal" OR "temporal"`
20. `metagenomics "spurious ORFs" OR "false positive" "protein families" OR "dark matter"`

---

## 3. Verified Core Literature (12 Detailed Papers)

Below are the 12 verified core papers evaluating unknown microbial protein families, catalogue methodologies, ecosystem specificity, and longitudinal sampling. All 12 papers have been fully downloaded and verified in [`resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/).

### 1. Vanni et al. (2022) — The AGNOSTOS Framework
- **Title:** Unifying the known and unknown microbial coding sequence space
- **Authors:** Chiara Vanni, Matthew S. Schechter, Silvia G. Acinas, Albert Barberán, Pier Luigi Buttigieg, Emilio O. Casamayor, Tom O. Delmont, Carlos M. Duarte, A. Murat Eren, Robert D. Finn, Renzo Kottmann, Alex Mitchell, Pablo Sánchez, Kimmo Siren, Martin Steinegger, Frank Oliver Glöckner, Antonio Fernández-Guerra
- **Year / Journal:** 2022 | *eLife*, 11:e67667
- **DOI:** [`10.7554/eLife.67667`](https://doi.org/10.7554/eLife.67667)
- **Publisher URL:** [https://elifesciences.org/articles/67667](https://elifesciences.org/articles/67667)
- **Full-Text Local Path:** [`Vanni_2022_AGNOSTOS_eLife_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Vanni_2022_AGNOSTOS_eLife_preprint.pdf) (5.2 MB) & [`Vanni_2022_AGNOSTOS_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Vanni_2022_AGNOSTOS_fulltext.xml) (253 KB)
- **Already Cited:** Yes (`vanni2022agnostos` in `manuscript/references.bib:13`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Established AGNOSTOS to systematically partition coding sequence space into Knowns with Pfam, Knowns without Pfam, Genomic Unknowns (GUs, found in genomes but lacking annotations), and Environmental Unknowns (EUs, found only in metagenomes). Documented that ~30% of microbial genes are unknown and that unknown gene clusters show high lineage and habitat restriction.
  - *Limitation / Challenge for P4:* More than **80% of environmental unknown gene clusters in AGNOSTOS could not be associated with any MAG**, leaving them as unbinned, fragmented contigs with unverified organismal origin. P4 overcomes this by requiring that all 4.29M candidate dark representatives sit on contigs belonging to HQ+MQ bins (genome-resolved), demonstrating that dark matter can be reliably tied to organismal lineages and traced across space and time.

### 2. Rodríguez del Río et al. (2024) — The FESNov Catalogue
- **Title:** Functional and evolutionary significance of unknown genes from uncultivated taxa
- **Authors:** Álvaro Rodríguez del Río, Joaquín Giner-Lamia, Carlos P. Cantalapiedra, Jorge Botas, Ziqi Deng, Ana Hernández-Plaza, Martí Munar-Palmer, Saray Santamaría-Hernando, José J. Rodríguez-Herva, Hans-Joachim Ruscheweyh, Lucas Paoli, Thomas S. B. Schmidt, Shinichi Sunagawa, Peer Bork, Emilia López-Solanilla, Luis Pedro Coelho, Jaime Huerta-Cepas
- **Year / Journal:** 2024 | *Nature*, 626(7998):377–384
- **DOI:** [`10.1038/s41586-023-06955-z`](https://doi.org/10.1038/s41586-023-06955-z)
- **Publisher URL:** [https://www.nature.com/articles/s41586-023-06955-z](https://www.nature.com/articles/s41586-023-06955-z)
- **Full-Text Local Path:** [`RodriguezdelRio_2024_FESNov_Nature.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/RodriguezdelRio_2024_FESNov_Nature.pdf) (4.9 MB) & [`RodriguezdelRio_2024_FESNov_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/RodriguezdelRio_2024_FESNov_fulltext.xml) (116 KB)
- **Already Cited:** Yes (`rodriguezdelrio2024fesnov` in `manuscript/references.bib:77`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Compiled 404,085 FESNov gene families exclusive to uncultivated bacteria and archaea from ~150,000 environmental genomes; demonstrated that they evolve under purifying selection ($dN/dS < 1$), triple known prokaryotic family diversity, and provide lineage-specific synapomorphic markers for uncultivated clades.
  - *Limitation / Challenge for P4:* FESNov conceptualizes unknown gene families as markers of deeply uncultivated, exotic taxonomic branches (Candidate Phyla Radiation, etc.). P4 challenges this by demonstrating that in hyperarid terrestrial soils, dark families are overwhelmingly carried by ordinary, dominant named phyla (Actinomycetota 32%, Thermoplasmatota 19%), and constitute the **strain-variable accessory genome** rather than lineage-defining core markers (proven by the 80% annotated vs 16% dark sharing between cultured strains and MAGs of the identical species).

### 3. Coelho et al. (2022) — The Global Microbial Gene Catalog (GMGC)
- **Title:** Towards the biogeography of prokaryotic genes
- **Authors:** Luis Pedro Coelho, Renato Alves, Álvaro Rodríguez del Río, Pernille Neve Myers, Carlos P. Cantalapiedra, Joaquin Giner-Lamia, Thomas S. B. Schmidt, Daniel R. Mende, Askarbek Orakov, Ivica Letunic, Falk Hildebrand, Thea Van Rossum, Sofia K. Forslund, Supriya Khedkar, Oleksandr M. Maistrenko, Shaojun Pan, Longhao Jia, Pamela Ferretti, Shinichi Sunagawa, Xing-Ming Zhao, Henrik Bjørn Nielsen, Jaime Huerta-Cepas, Peer Bork
- **Year / Journal:** 2022 | *Nature*, 601(7892):252–256
- **DOI:** [`10.1038/s41586-021-04233-4`](https://doi.org/10.1038/s41586-021-04233-4)
- **Publisher URL:** [https://www.nature.com/articles/s41586-021-04233-4](https://www.nature.com/articles/s41586-021-04233-4)
- **Full-Text Local Path:** [`Coelho_2022_GMGC_Nature_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Coelho_2022_GMGC_Nature_preprint.pdf) (667 KB) & [`Coelho_2022_GMGC_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Coelho_2022_GMGC_fulltext.xml) (82 KB)
- **Already Cited:** Yes (`coelho2022gmgc` in `manuscript/references.bib:721`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Clustered 2 billion ORFs from 13,174 metagenomes into 303 million unigenes across 14 habitats; established that the vast majority of genes are rare and habitat-specific, with only a tiny fraction showing cosmopolitan distribution.
  - *Limitation / Challenge for P4:* Cross-sectional pooling cannot distinguish whether a "rare" gene observed in 1–2 samples is genuine, resident biology or transient relic DNA, moribund air-spillover, or an assembly artifact. P4 resolves this limitation by substituting cross-sectional rarity with longitudinal persistence, demonstrating that 79.4% of Rub' al-Khali dark families recur across all five field expeditions spanning 31 months.

### 4. Pavlopoulos et al. (2023) — Global Metagenomic Functional Dark Matter
- **Title:** Unraveling the functional dark matter through global metagenomics
- **Authors:** Georgios A. Pavlopoulos, Fotis A. Baltoumas, Sirui Liu, Oguz Selvitopi, Antonio Pedro Camargo, Stephen Nayfach, Ariful Azad, Simon Roux, Lee Call, Natalia N. Ivanova, I-Min A. Chen, David Paez-Espino, Evangelos Karatzas, Ioannis Iliopoulos, Konstantinos Konstantinidis, James M. Tiedje, Jennifer Pett-Ridge, David Baker, Axel Visel, Christos A. Ouzounis, Sergey Ovchinnikov, Aydin Buluc, Nikos C. Kyrpides
- **Year / Journal:** 2023 | *Nature*, 622(7983):594–602
- **DOI:** [`10.1038/s41586-023-06583-7`](https://doi.org/10.1038/s41586-023-06583-7)
- **Publisher URL:** [https://www.nature.com/articles/s41586-023-06583-7](https://www.nature.com/articles/s41586-023-06583-7)
- **Full-Text Local Path:** [`Pavlopoulos_2023_darkmatter_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Pavlopoulos_2023_darkmatter_fulltext.xml) (207 KB)
- **Already Cited:** Yes (`pavlopoulos2023darkmatter` in `manuscript/references.bib:7`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Clustered 3 billion metagenomic proteins from 26,931 metagenomes into 1.17 billion non-redundant sequences; identified 106,198 novel protein families (each with $\ge 100$ members) with no Pfam or reference genome matches; predicted 3D folds and genomic neighborhoods.
  - *Limitation / Challenge for P4:* Enforcing a strict cutoff of $\ge 100$ sequences across global pooled metagenomes filters out localized, ecosystem-specific, and low-biomass protein families. P4 shows that an admission rule of $\ge 3$ diverged variants across $\ge 3$ samples and $\ge 2$ sites in a single hyperarid transect uncovers 106,042 robust, purifying-selected dark families that are fundamentally missing from Pavlopoulos et al.'s global catalogue (52.68% absent from MGnify90; 28.74% absent from all global databases combined).

### 5. Aplakidou et al. (2026) — "Quadrupling" the Protein Family Space (NAR 2026)
- **Title:** “Quadrupling” the protein family space with global metagenomics
- **Authors:** Eleni Aplakidou, Fotis A. Baltoumas, Maria N. Chasapi, Efthimia Lamari, Iro N. Chasapi, Ilias Georgakopoulos-Soares, Evangelos Karatzas, Ioannis Iliopoulos, Aydin Buluç, Robert D. Finn, Antonio Pedro Camargo, Nikos C. Kyrpides, Georgios A. Pavlopoulos
- **Year / Journal:** 2026 | *Nucleic Acids Research*, 54(18):gkag938 (October 2026)
- **DOI:** [`10.1093/nar/gkag938`](https://doi.org/10.1093/nar/gkag938) *(bioRxiv preprint DOI: 10.64898/2025.12.21.695754)*
- **Publisher URL:** [https://academic.oup.com/nar/article/doi/10.1093/nar/gkag938](https://academic.oup.com/nar/article/doi/10.1093/nar/gkag938)
- **Full-Text Local Path:** [`Aplakidou_2026_quadrupling_NAR_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Aplakidou_2026_quadrupling_NAR_preprint.pdf) (5.5 MB)
- **Already Cited:** Cited ONLY as a bioRxiv preprint (`aplakidou2025quadrupling` in `references.bib:1097`). The formal October 2026 *Nucleic Acids Research* journal article is new!
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Analyzed 40.3 billion ORFs across 40,446 metagenomes and 9,540 metatranscriptomes, recovering 608,258 novel protein families ($\ge 100$ members) and 6.5 million smaller families ($\ge 25$ members). Reported that **63.4% of novel families are restricted to a single ecosystem** (under GOLD 3/8-level classification).
  - *Limitation / Challenge for P4:* Aplakidou et al. confirms that single-ecosystem restriction is the norm across microbial sequence space (directly supporting P4's finding of 28.74% complete absence from global catalogues). However, because its ecosystem labels are broad GOLD classes (e.g. "Terrestrial: Soil"), it aggregates disparate soil types, ignoring hyperaridity and micro-niche chemistry. Furthermore, it operates on cross-sectional data with no longitudinal recurrence validation.

### 6. Baltoumas et al. (2024) — NMPFamsDB
- **Title:** NMPFamsDB: a database of novel protein families from microbial metagenomes and metatranscriptomes
- **Authors:** Fotis A. Baltoumas, Evangelos Karatzas, Maria N. Chasapi, Ioannis Iliopoulos, Nikos C. Kyrpides, Georgios A. Pavlopoulos
- **Year / Journal:** 2024 | *Nucleic Acids Research*, 52(D1):D502–D512
- **DOI:** [`10.1093/nar/gkad800`](https://doi.org/10.1093/nar/gkad800) *(PMCID: PMC10767849)*
- **Publisher URL:** [https://academic.oup.com/nar/article/52/D1/D502/7330925](https://academic.oup.com/nar/article/52/D1/D502/7330925)
- **Full-Text Local Path:** [`Baltoumas_2024_NMPFamsDB_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Baltoumas_2024_NMPFamsDB_fulltext.xml) (111 KB)
- **Already Cited:** No (New reference).
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Established the dedicated database resource housing the 106,198 novel protein families from Pavlopoulos et al., providing searchable MSAs, HMM profiles, AlphaFold2 structural coordinates, and ecosystem distributions.
  - *Limitation / Challenge for P4:* NMPFamsDB is the global benchmark for uncharacterized metagenomic protein families. However, its database entries were derived from pooled contigs without requiring binning into MAGs, and its search interface does not account for localized desert endemism or strain-level pangenomic turnover.

### 7. Karatzas et al. (2026) — nf-core/proteinfamilies
- **Title:** nf-core/proteinfamilies: a scalable pipeline for the generation of protein families
- **Authors:** Evangelos Karatzas, Martin Beracochea, Fotis A. Baltoumas, Eleni Aplakidou, Lorna Richardson, James A. Fellows Yates, Daniel Lundin, nf-core community, Aydin Buluç, Nikos C. Kyrpides, Ilias Georgakopoulos-Soares, Georgios A. Pavlopoulos, Robert D. Finn
- **Year / Journal:** 2026 | *GigaScience*, 15:giag009
- **DOI:** [`10.1093/gigascience/giag009`](https://doi.org/10.1093/gigascience/giag009)
- **Publisher URL:** [https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giag009](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giag009)
- **Full-Text Local Path:** [`Karatzas_2026_proteinfamilies_GigaScience_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Karatzas_2026_proteinfamilies_GigaScience_preprint.pdf) (2.8 MB)
- **Already Cited:** No (New reference).
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Standardized de novo protein family generation into a Nextflow pipeline combining MMseqs2 clustering, FAMSA/MAFFT multiple alignment, and HMMER profile generation with linear scaling.
  - *Limitation / Challenge for P4:* Validates that the choice of clustering identity and coverage mode (e.g. MMseqs2 `--min-seq-id 0.5 -c 0.8`) is the critical determinant of family boundary definitions. Directly relevant to P4's robustness test (§1.4), which reclustered the 4.29M dark representatives at 30% identity to demonstrate that the dark universe does not artificially collapse or fragment at different clustering thresholds (moving smoothly from 108,826 to 130,843 families).

### 8. Chasapi et al. (2026) — metagRoot
- **Title:** metagRoot: a comprehensive database of protein families associated with plant root microbiomes
- **Authors:** Maria N. Chasapi, Fotis A. Baltoumas, Alexandros Galaras, Efthimia Lamari, Iro N. Chasapi, Eleni Aplakidou, Evangelos Karatzas, Ilias Georgakopoulos-Soares, Nikos C. Kyrpides, Georgios A. Pavlopoulos
- **Year / Journal:** 2026 | *Nucleic Acids Research*, 54(D1):D1733–D1742
- **DOI:** [`10.1093/nar/gkaf862`](https://doi.org/10.1093/nar/gkaf862) *(PMCID: PMC12808537)*
- **Publisher URL:** [https://academic.oup.com/nar/article/54/D1/D1733/8282361](https://academic.oup.com/nar/article/54/D1/D1733/8282361)
- **Full-Text Local Path:** [`Chasapi_2026_metagRoot_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Chasapi_2026_metagRoot_fulltext.xml) (111 KB)
- **Already Cited:** No (New reference).
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Demonstrated that isolating a specific plant-associated sub-habitat (rhizosphere/root) yields 71,091 enriched protein families ($\ge 100$ members) with distinct structural folds and functional signatures that are overlooked in bulk planetary catalogues.
  - *Limitation / Challenge for P4:* Supports P4's micro-habitat sampling design (surface sand vs. subsurface vs. root-associated rhizosphere). However, metagRoot pools cross-sectional root metagenomes globally without longitudinal monitoring, missing the temporal stability of root-associated dark families under extreme hyperarid conditions.

### 9. Galaras et al. (2026) — WasteFams
- **Title:** WasteFams: A database of protein families from global wastewater microbiomes
- **Authors:** Alexandros Galaras, Iro N. Chasapi, Eleni Aplakidou, Maria N. Chasapi, Efthimia Lamari, Sophia Diplari, Ilias Georgakopoulos-Soares, Evangelos Karatzas, Fotis A. Baltoumas, Nikos C. Kyrpides, Georgios A. Pavlopoulos
- **Year / Journal:** 2026 | *bioRxiv*, preprint
- **DOI:** [`10.64898/2026.05.08.723720`](https://doi.org/10.64898/2026.05.08.723720)
- **Publisher URL:** [https://www.biorxiv.org/content/10.64898/2026.05.08.723720v1](https://www.biorxiv.org/content/10.64898/2026.05.08.723720v1)
- **Full-Text Local Path:** [`Galaras_2026_WasteFams_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Galaras_2026_WasteFams_preprint.pdf) (12.7 MB)
- **Already Cited:** No (New reference).
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Highlights that extreme, highly specialized biomes harbor unique sets of unknown protein families frequently associated with mobile genetic elements, resistance mechanisms, and environmental stress defense.
  - *Limitation / Challenge for P4:* Wastewater is an open, high-flux artificial system where continuous human/industrial influx confounds resident adapted genes with transient transients. In contrast, the Rub' al-Khali is an undisturbed hyperarid natural end-member where dark families persist over evolutionary timeframes under severe moisture and nutrient scarcity.

### 10. Schmidt et al. (2024) — SPIRE
- **Title:** SPIRE: a searchable, planetary-scale microbiome resource
- **Authors:** Thomas S. B. Schmidt, Anthony Fullam, Pamela Ferretti, Askarbek Orakov, Oleksandr M. Maistrenko, Hans-Joachim Ruscheweyh, Ivica Letunic, Yiqian Duan, Thea Van Rossum, Shinichi Sunagawa, Daniel R. Mende, Robert D. Finn, Michael Kuhn, Luis Pedro Coelho, Peer Bork
- **Year / Journal:** 2024 | *Nucleic Acids Research*, 52(D1):D777–D783
- **DOI:** [`10.1093/nar/gkad943`](https://doi.org/10.1093/nar/gkad943) *(PMCID: PMC10767986)*
- **Publisher URL:** [https://academic.oup.com/nar/article/52/D1/D777/7342080](https://academic.oup.com/nar/article/52/D1/D777/7342080)
- **Full-Text Local Path:** [`Schmidt_2024_SPIRE_NAR_preprint.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Schmidt_2024_SPIRE_NAR_preprint.pdf) (3.7 MB) & [`Schmidt_2024_SPIRE_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Schmidt_2024_SPIRE_fulltext.xml) (101 KB)
- **Already Cited:** Yes (`schmidt2024spire` in `manuscript/references.bib:474`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Reconstructed 102,374 MAGs from ~500,000 metagenomes to create a planetary-scale microbiome database with 269M representative proteins.
  - *Limitation / Challenge for P4:* SPIRE serves as one of the primary genome-resolved comparators in P4 (§3.1). Despite its massive size, 39.9% of Rub' al-Khali dark families have no hit in SPIRE. In P4's rarefaction analysis, SPIRE leaves **+21.30 percentage points less unexplained** than expected from a generic database of equivalent size, proving that genome-resolved reference catalogues still exhibit severe habitat-specific dropouts when hyperarid desert interiors are omitted.

### 11. Choudoir et al. (2025) — Soil Pangenomics and Accessory Genome Evolution
- **Title:** Pangenomes suggest ecological-evolutionary responses to experimental soil warming
- **Authors:** Mallory J. Choudoir, Achala Narayanan, Damayanti Rodriguez-Ramos, Rachel Simoes, Alon Efroni, Abigail Sondrini, Kristen M. DeAngelis
- **Year / Journal:** 2025 | *mSphere*, 10(2):e00059-25
- **DOI:** [`10.1128/msphere.00059-25`](https://doi.org/10.1128/msphere.00059-25) *(PMCID: PMC12039271)*
- **Publisher URL:** [https://journals.asm.org/doi/10.1128/msphere.00059-25](https://journals.asm.org/doi/10.1128/msphere.00059-25)
- **Full-Text Local Path:** [`Choudoir_2025_soilwarming_fulltext.xml`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Choudoir_2025_soilwarming_fulltext.xml) (208 KB)
- **Already Cited:** No (New reference).
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Demonstrated via comparative pangenomics of dominant soil bacterial isolates (including Actinobacteriota) that ecological and climate adaptations are concentrated in the **accessory genome**, where variable gene clusters undergo rapid turnover while core housekeeping genes remain static.
  - *Limitation / Challenge for P4:* Choudoir et al. investigated cultured isolates in a single temperate forest warming plot. P4 scales this conceptual framework to 5,438 uncultivated MAGs across 60 hyperarid sites, providing the empirical proof that the unannotated "dark proteome" *is* the strain-variable accessory genome of dominant soil taxa (which explains why a few isolate reference genomes fail to capture 84% of their dark genes).

### 12. Duan et al. (2024) — Global Small Protein Dark Matter (GMSC)
- **Title:** A catalog of small proteins from the global microbiome
- **Authors:** Yiqian Duan, Célio Dias Santos-Júnior, Thomas S. B. Schmidt, Anthony Fullam, Breno L. S. de Almeida, Chengkai Zhu, Michael Kuhn, Xing-Ming Zhao, Peer Bork, Luis Pedro Coelho
- **Year / Journal:** 2024 | *Nature Communications*, 15:7563
- **DOI:** [`10.1038/s41467-024-51894-6`](https://doi.org/10.1038/s41467-024-51894-6)
- **Publisher URL:** [https://www.nature.com/articles/s41467-024-51894-6](https://www.nature.com/articles/s41467-024-51894-6)
- **Full-Text Local Path:** [`Duan_2024_GMSC_NatCommun.pdf`](file:///home/leechuck/Documents/papers/empty-quarter-mags/resources/literature/darkmatter-2026-10-13/papers/01_dark_catalogues/Duan_2024_GMSC_NatCommun.pdf) (2.5 MB)
- **Already Cited:** Yes (`duan2024gmsc` in `manuscript/references.bib:480`)
- **Precise Claim / Limitation Relevant to P4:**
  - *Claim:* Identified 96 million small open reading frames (smORFs, $\le 100$ amino acids) from global metagenomes (GMSC), finding that small proteins are dramatically enriched in uncharacterized sequence space and display extreme habitat specificity.
  - *Limitation / Challenge for P4:* In P4, 76% of dark representatives are shorter than 100 residues (compared to 30% of annotated representatives). However, Duan et al. pooled cross-sectional samples, leaving open the skeptical critique that small environmental ORFs are spurious gene-caller artifacts. P4 refutes this by proving that short dark families persist across multi-year expeditions (89.7% of families $<60$ aa detected in all 5 expeditions) and display strong purifying selection ($dN/dS < 1$).

---

## 4. Unresolved Leads & Ongoing Catalogues

To preserve literature integrity, the following leads were discovered during search sweeps but remain unresolved or pending formal validation:
1. **GMGCv2 / Pan-Microbial Unigene Expansions:** The Big Data Biology Lab has discussed ongoing expansions of GMGC (integrating longer read assemblies and deeper marine/soil datasets). However, as of late 2026, no standalone publication or versioned data release replacing GMGCv1.0 has been published. Status: *Unresolved lead; retain GMGCv1 (Coelho et al. 2022) as primary.*
2. **ArthroVerse (Chasapi et al. 2026, bioRxiv: 10.64898/2026.07.11.737898):** A companion database in the EnvoFams ecosystem mapping protein family diversity in arthropod-associated microbiomes. While technically interesting, its biological scope (insect and arachnid microbiomes) has minimal direct bearing on hyperarid soil edaphic ecology. Status: *Deprioritized lead.*
3. **Metaproteomics Dark Matter Benchmarks (CAMPI / uMetaP):** Several recent preprints from the Metaproteomics Initiative benchmark ultra-sensitive mass spectrometry and *de novo* peptide sequencing to detect translated dark proteins in low-biomass soils. These provide exciting future validation for expressing dark ORFs, but current datasets lack paired longitudinal metagenomic bins. Status: *Methodological lead.*

---

## 5. Challenging Novelty Claims & Skeptical Analysis

Reviewers frequently challenge claims of massive "functional dark matter" in metagenomics. Below is a structured critique evaluating how prior literature addresses these challenges and where P4's evidence stands:

### Challenge 1: Spurious Open Reading Frames (ORFs) vs. True Genes
* **The Skeptic's Attack:** Standard gene callers (Prodigal, Pyrodigal, MetaGeneMark) over-predict short ORFs, particularly on fragmented metagenomic contigs. Many "novel" families are simply recurrent non-coding DNA, intergenic fragments, or shadow ORFs.
* **Literature Status:** AGNOSTOS (Vanni et al. 2022) filters against AntiFam and shadow ORFs, but still retains hundreds of thousands of unbinned clusters. Pavlopoulos et al. (2023) relies on $\ge 100$ members to filter noise, which conflates family size with coding validity.
* **P4's Resolution:** P4 introduces a 4-leg defense:
  1. Admitted families require $\ge 3$ diverged variants independently assembled across $\ge 3$ samples and $\ge 2$ sites (assembly noise does not independently evolve).
  2. Pairwise $dN/dS$ over 97,200 dark families yields a median of $0.425$ ($96.7\% < 1$).
  3. A shifted +1 nucleotide frame control yields median $dN/dS = 1.087$ ($92.4\%$ of families have higher $dN/dS$ in the shifted frame than the annotated frame).
  4. The median $0.425$ is calibrated against an identically clustered, length- and size-matched UniRef50-matched control ($0.388$).

### Challenge 2: Rarity in Pooled Catalogues vs. Longitudinal Persistence
* **The Skeptic's Attack:** In GMGC (Coelho et al. 2022) and Pavlopoulos et al. (2023), most novel genes are found in only one or two samples. These could be moribund relic DNA (Carini et al. 2016) or transient airborne contaminants rather than living, reproducing organisms.
* **Literature Status:** Cross-sectional catalogues cannot test persistence because each site is sampled exactly once.
* **P4's Resolution:** Longitudinal sampling across 5 expeditions over 31 months proves that 79.4% of dark families are detected by read mapping in all 5 expeditions (and 28.9% independently assembled in all 5 trips). Dead or transient contamination does not consistently re-assemble and recruit reads over multi-year intervals.

### Challenge 3: Uncultivated Phyla Markers vs. Strain-Variable Accessory Genomes
* **The Skeptic's Attack:** FESNov (Rodríguez del Río et al. 2024) asserts that unknown genes represent clade-defining markers (synapomorphies) of uncultivated, deep phyla (CPR, etc.). If dark genes belong to ordinary bacteria, why weren't they discovered in sequenced isolates?
* **Literature Status:** FESNov focused heavily on cross-sectional MAGs from uncultivated clades, creating an observational bias that darkness equals deep taxonomic novelty.
* **P4's Resolution:** P4 demonstrates that 97.9% of attributable dark families belong to ordinary, dominant soil phyla (Actinomycetota 32%, Thermoplasmatota 19%). By comparing MAGs with 23 strains cultured from the exact same basin, P4 reveals that isolates and MAGs share 80% of annotated genes but only 16% of dark genes. Dark matter is the open, strain-variable accessory genome that isolate sequencing misses.

---

## 6. Recommendations for P4 Manuscript & BibTeX Updates

1. **Update `aplakidou2025quadrupling`:** Update the citation from the bioRxiv 2025 preprint to the formal October 2026 *Nucleic Acids Research* publication (`10.1093/nar/gkag938`).
2. **Incorporate `baltoumas2024nmpfamsdb` and `karatzas2026proteinfamilies`:** Cite NMPFamsDB (`10.1093/nar/gkad800`) as the reference repository for global metagenomic dark families, and Karatzas et al. (2026) when discussing family clustering stability and Nextflow pipeline reproducibility.
3. **Contrast Single-Ecosystem Restriction:** Contrast Aplakidou et al.'s 63.4% single-ecosystem rate with P4's finding that 28.74% of dark families are absent from all global databases combined. Note that Aplakidou's GOLD categories are too coarse to capture hyperarid soil edaphic specialization.
4. **Leverage Choudoir et al. (2025):** Cite Choudoir et al. (`10.1128/msphere.00059-25`) to reinforce that soil microbial environmental adaptation is concentrated in the accessory genome, bridging isolate pangenomics with P4's metagenome-resolved dark proteome.
