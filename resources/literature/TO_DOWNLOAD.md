# To download (relevant papers this IBEX node could NOT fetch — PMC/paywall-blocked)

The compute node can pull bioRxiv + nature.com + author PDFs, but **PMC and hard-paywalled journals
(Science, Cell, older Nature Microbiology, some OUP) are blocked**. Below are the relevant papers not yet
retrieved. Download each and drop it into `resources/literature/` (any filename); the agent identifies,
renames, and files them.

Identifiers are best-known; confirm the DOI on the publisher page before it reaches a reference list.

| Paper | ID | → folder | why it matters |
|---|---|---|---|
| Makhalanyane et al. 2015, hot-desert edaphic microbial ecology, *FEMS Microbiol Rev* | doi:10.1093/femsre/fuu011 | `survival-genome/` | hypolith/subsurface light regimes — decides §2.3 photolyase framing (sun exposure vs surface layer) |
| Roberts 2005, organic compatible solutes of halotolerant/halophilic microorganisms, *Saline Systems* | doi:10.1186/1746-1448-1-5 (OA) | `survival-genome/` | solute cost / water-activity ordering for §2.3's salt half (trehalose vs ectoine vs GG/MG hierarchy) |
| Ji et al. 2017, atmospheric trace gases support primary production, *Nature* 552:400 | doi:10.1038/nature25014 | `trace-gas-deserts/` | cited in §2.2 as the founding precedent but **never read** — the one load-bearing citation in that subsection with no PDF on disk. Ortiz 2021 (now keyed, OA, same consortium) can carry the claim if it stays unread |
| Coleine et al. 2024, Antarctic endolithic bacteria, *Sci Total Environ* 917:170290 | doi:10.1016/j.scitotenv.2024.170290 | `desert-mags-comparators/` | **only the 2023 preprint is on disk.** The 41.1%-of-MAGs high-affinity hydrogenase figure and the 2,278-OTU/98.3%-novel counts we would cite in §2.2/§1.1 come from the preprint and must be re-checked against the journal version (title changed between the two) |
| Baran et al. 2017, extensive turnover of compatible solutes in cyanobacteria (D₂O SIP), *ACS Chem Biol* 12:674–681 | (reference 70 of Moore 2023; DOI not yet verified) | `survival-genome/` | the **primary** source for the glycogen→compatible-solute claim that §2.3 wants as a caveat. Moore 2023 only cites it secondhand and extends it via a metabolic model; if the caveat becomes load-bearing, cite Baran directly. Low priority while the sentence stays a one-line limitation |
| D'Costa et al. 2011, "Antibiotic resistance is ancient", *Nature* 477:457--461 | doi:10.1038/nature10388 (Crossref-verified; **no PMCID**, PMID 21881561) | `contamination-armor/` | cited at §2.4 L190 for the intrinsic/pre-anthropogenic resistome. Abstract verified from this node (30,000-yr Beringian permafrost; beta-lactam/tetracycline/glycopeptide ARGs; complete VanA), but Nature is hard-paywalled so no page-level quote is available. Low priority: `vangoethem2018antarctic` is OA (`europepmc.org/articles/PMC5824556?pdf=render`) and is the closer match to the soil-competition wording; D'Costa only needs its PDF if the ancient/pre-anthropogenic half becomes load-bearing |

_Rows added 2026-08-12 (§2.3 consult) and 2026-08-12 batch-2 ingest. Makhalanyane and Roberts DOIs remain
UNVERIFIED against publisher pages; the Ji, Coleine and Baran identifiers were verified this session
(Crossref / the citing PDF's own reference list)._

**Resolved without the PDF (2026-08-12):**
- **Bay et al. 2021** (*Nat Microbiol* 6:246–256, doi:10.1038/s41564-020-00811-w) — abstract retrieved and
  verified via the Europe PMC REST API; the numbers §2.2 needed are recorded in READING_NOTES.md, including
  the scope trap (">70%" is *inorganic energy sources*, not trace gases; the trace-gas claim is "19 phyla").
  The PDF is still worth having eventually for the thermodynamic-modelling detail, but it no longer blocks
  §2.2. Row retired from the pending table.
- **Sharda et al. 2020**, bacterial NHEJ across ~6,000 genomes (*Genome Biol Evol* 12:2450–2466) — **retrieved
  this session** via the Europe PMC render route (`europepmc.org/articles/PMC7719229?pdf=render`), which works
  from this node even though `ncbi.nlm.nih.gov/pmc` does not. Filed at
  `survival-genome/Sharda_2020_bacterial_NHEJ_evolution_GBE.pdf`. This supplies the 22% NHEJ baseline the
  §2.3 \todo asked for. **Worth remembering: the EPMC render URL is a working route to OA PMC PDFs.**

## Downloaded 2026-08-12 (previous pending list fully cleared — filed and renamed)

| final path | paper |
|---|---|
| `structure-darkmatter/RodriguezdelRio_2024_FESNov_unknown_genes_Nature.pdf` | Rodríguez del Río et al. 2024, FESNov, *Nature* 626:377 |
| `structure-darkmatter/Lau_2024_TED_encyclopedia_domains_Science.pdf` | Lau et al. 2024, TED, *Science* (adq4946) |
| `structure-darkmatter/Waman_2025_CATHv44_expansion_NAR.pdf` | Waman et al. 2025, CATH v4.4, *NAR* 53:D348 |
| `structure-darkmatter/Harteveld_2024_darkmatter_folds_deeplearning_CellSystems.pdf` | Harteveld et al. 2024, Genesis VAE dark folds, *Cell Systems* |
| `contamination-armor/Carini_2016_relicDNA_soil_NatMicro.pdf` | Carini et al. 2016, relic DNA, *Nat Microbiol* |
| `reference-bias/Sun_2020_PICRUSt_accuracy_sampletypes_Microbiome.pdf` | Sun, Jones & Fodor 2020, *Microbiome* |
| `reference-bias/Matchado_2024_16S_prediction_limits_MicrobGenom.pdf` | Matchado et al. 2024, *Microbial Genomics* 10:001203 |
| `reference-bias/Wright_2025_PICRUSt2SC_reference_update_Bioinformatics.pdf` | Wright & Langille 2025, PICRUSt2-SC, *Bioinformatics* (btaf269) |
| `pangenome-decoupling/Whelan_2021_accessory_selection_pangenome_MBE.pdf` | Whelan, Hall & McInerney 2021, *Mol Biol Evol* (msab139) |
| `pangenome-decoupling/Guerra_2026_pangenome_statistical_model_BioinfAdv.pdf` | Guerra 2026, *Bioinformatics Advances* 6(1) vbag069 (list previously said 2025 — actual 2026) |
| `metabolism-bgc/AndreaniGerard_2024_Atacama_soil_BGCs_mSphere.pdf` | Andreani-Gerard et al. 2024, *mSphere* (filed under metabolism-bgc, not desert-mags-comparators: §2.4 BGC comparator next to Waschulin 2021) |
| `desert-mags-comparators/Arabian_2026_extreme_microbiomes_astrobiology_Extremophiles.pdf` | Extreme Arabian environments & their microbiomes review, *Extremophiles* 30:13, 2026 (KAUST/Rosado group) |
| `dryland-ecology/Zhou_2025_lifehistory_aridity_multifunctionality_PNAS.pdf` | Zhou et al. 2025, *PNAS* 122(41) e2511071122 |
| `smorfs-smallproteins/Li_2026_smallproteins_fulloceandepth_iScience.pdf` | Li, He & Wang, *iScience* 29:114585 (2026, not 2025 as previously listed) |

_Earlier arrivals (2026-08-09 routes): AFESM 2025, Qaidam 2025, GMSC 2024, kitome 2025, low-biomass
guidelines 2025, Antarctic endolithic 2023, Delgado-Baquerizo aridity 2019, Atacama desiccation 2023._

**Ingest status (2026-08-12).** All 14 papers in the table above, plus the three 2026-08-09 arrivals that had
not yet been distilled (kitome, Antarctic endolithic, Atacama desiccation), are now read and recorded in
`READING_NOTES.md` under "Batch 2", with bib entries in `paper/references.bib` and `manuscript/references.bib`
(the two files are byte-identical). Five June trace-gas PDFs that had never been distilled were also
verified and keyed: Ehhalt & Rohrer 2009, Ocko & Hamburg 2022, Buzzard 2022, Ortiz 2021, Hutchinson 2025.
Nothing in `resources/literature/` is now un-ingested.

_Row added 2026-08-18 (citation audit, session 8). Also confirmed working from this node: Crossref REST for metadata, the Europe PMC render route for OA PDFs, and **bioRxiv `.full` HTML** -- the last rate-limits with HTTP 429, so use `curl --retry 4 --retry-delay 40 --retry-all-errors`. Note bioRxiv's **new DOI prefix `10.64898`** for late-2025 preprints; the old `10.1101` form does not resolve._

## Added 2026-08-18 (session 8, citations sweep) — claim-bearing cites that could NOT be verified at T1

Full audit: `resources/literature/CITATION_AUDIT.md` §5. These are the papers whose *claim support* could not be
checked because no PDF is on disk and the fulltext is not reachable from this node. Metadata for all of them is
Crossref-verified; it is the content that is unread.

| paper | identifier | why it is load-bearing |
|---|---|---|
| Edgell 2006, *Arabian Deserts*, Springer | 10.1007/1-4020-3970-0 — ch. 7 "Arabian sand seas" pp. 125-166 = 10.1007/1-4020-3970-0_7 | the "world's largest contiguous sand sea" claim (Intro L79). Springer fulltext auth-walled |
| Huang et al. 2016, *Nat Clim Change* 6:166-171 | 10.1038/nclimate2837 | "drylands are the largest terrestrial biome" (audit D-7). nature.com blocked; no abstract via OpenAlex/S2 either |
| Ji et al. 2017, *Nature* 552:400-403 | 10.1038/nature25014 | load-bearing for S2.2's energy argument (audit D-5). Abstract only |
| Greening et al. 2022, *Nat Rev Microbiol* 20:513-528 | 10.1038/s41579-022-00724-x | carries the "absent from global budgets" negative (audit D-3). Paywalled, no PMCID |
| King & Weber 2007, *Nat Rev Microbiol* 5:107-118 | 10.1038/nrmicro1595 | the CO-oxidiser breadth claim (audit D-4). Paywalled review |
| Constant et al. 2010, *Environ Microbiol* 12:821-829 | 10.1111/j.1462-2920.2009.02130.x | "classic mesic-soil group 1h" (audit D-6). NOTE: the on-disk `Constant_2016_streptomyces_h2.pdf` is **Liot & Constant 2016**, a different paper |
| D'Costa et al. 2011, *Nature* 477:457-461 | 10.1038/nature10388 | pre-anthropogenic half of S2.4's resistome claim. Not OA; abstract sufficient for its cite-only role |
| Parks et al. 2017, *Nat Microbiol* 2:1533-1542 | 10.1038/s41564-017-0012-7 | optional 2nd anchor for the L134 deep-lineage expectation. Not added -- `liu2025qaidam` (on disk) carries it |

**Worth having but not blocking** (both on disk already, just unkeyed -- see audit D-3):
Reji et al., "Global implications of a low soil moisture threshold for microbial hydrogen uptake", 10.1038/s41467-025-67208-3
(`trace-gas-deserts/2025_soil_moisture_threshold_H2_uptake_NatComm.pdf`); Ray et al. 2023, *MMBR* 87:e0004823,
10.1128/mmbr.00048-23 (`trace-gas-deserts/Ray_2023_atmospheric_chemosynthesis_review.pdf`).

**Route note (confirmed again this session):** the Europe PMC render URL
`https://europepmc.org/articles/PMC<ID>?pdf=render` still works for OA PMC PDFs. **OpenAlex** (`api.openalex.org`)
resolves records Crossref does not -- it supplied Holm 1979, Mantel 1967 and the BWA-MEM preprint. arXiv's export
API and Semantic Scholar were NOT reachable (empty response / HTTP 429).

## Added 2026-08-18 — citation audit (session 8)

Every one of these is **cited in `paper/main.tex`** and could only be verified against an abstract or
metadata from this node. Full audit: `resources/literature/CITATION_AUDIT.md` §5. Priority order.

| paper | identifier | file into | why it is wanted |
|---|---|---|---|
| Greening & Grinter 2022, *Nat Rev Microbiol* 20:513–528 | doi:10.1038/s41579-022-00724-x (bib `greening2022review`) | `trace-gas-deserts/` | carries §2.2's "hot hyperarid drylands are **effectively absent** from global trace-gas budgets" — a strong negative claim the abstract cannot support (audit W-2). Paywalled, no PMCID |
| Ji et al. 2017, *Nature* 552:400–403 | doi:10.1038/nature25014 (bib `ji2017atmospheric`) | `trace-gas-deserts/` | cited **five** times; carries W-4 (cited for "maintenance rather than growth" while its title is "…support **primary production**"). Standing request, now load-bearing for the Introduction too |
| **Constant et al. 2010, *Environ Microbiol* 12:821–829** | doi:10.1111/j.1462-2920.2009.02123.x (bib `constant2010streptomycetes`) | `trace-gas-deserts/` | ⚠ **filing trap**: the on-disk `Constant_2016_streptomyces_h2.pdf` is **Liot & Constant 2016, *MicrobiologyOpen*** — a different paper. The cited 2010 paper is not on disk at all |
| D'Costa et al. 2011, *Nature* 477:457–461 | doi:10.1038/nature10388 (bib `dcosta2011ancient`) | `survival-genome/` or a new `resistome/` | newly cited this session for §2.4's intrinsic-resistome claim; abstract-verified only. Paywalled, no PMCID |
| King & Weber 2007, *Nat Rev Microbiol* 5:107–118 | doi:10.1038/nrmicro1595 (bib `king2007coxidizers`) | `trace-gas-deserts/` | carries W-3, the CO-vs-H₂ "more broadly distributed" comparison |
| Huang et al. 2016, *Nat Clim Change* 6:166–171 | doi:10.1038/nclimate2837 (bib `huang2016dryland`) | `dryland-ecology/` | opens the Introduction. **Abstract not retrievable at all** from this node — verified on metadata only (T3), the weakest tier of any claim-bearing cite |
| Edgell 2006, *Arabian Deserts*, Springer, **ch. 7 "Arabian sand seas" pp. 125–166** | book doi:10.1007/1-4020-3970-0; chapter doi:10.1007/1-4020-3970-0_7 | `desert-mags-comparators/` | the "world's largest contiguous sand sea" wording could not be confirmed (Springer redirects to auth). Needs a printed/library copy; cite the chapter, not the book |

**Routes confirmed working from this node (2026-08-18):** Crossref REST + content negotiation
(`curl -LH "Accept: application/x-bibtex" https://doi.org/<DOI>`), **OpenAlex** (`api.openalex.org` — the
fallback that resolved Holm 1979, Mantel 1967 and the BWA-MEM arXiv record, none of which Crossref has),
bioRxiv fulltext, and the Europe PMC render route. **Blocked:** `ncbi.nlm.nih.gov/pmc`, `export.arxiv.org`
(returns empty), nature.com and link.springer.com (auth redirects). Semantic Scholar rate-limits (429).

## Added 2026-09-22 — regional soil genome catalogues (L1 pass)

| paper | identifier | file into | why it is wanted |
|---|---|---|---|
| Leung, Bay, … Delgado-Baquerizo, Greening 2025, "Aridity drives global convergence of desert microbiomes and biogeochemical activities", bioRxiv (32 authors) | doi:10.1101/2025.11.29.691269 (bib `leung2025aridity`) | `regional-soil-catalogues/` or `trace-gas-deserts/` | the activity-measured successor to `ji2017atmospheric` / `bay2021trace`: standardised desert survey on seven continents with measured respiration, trace-gas oxidation and carbon fixation, aridity as the primary driver. Bears directly on §2.2 (our trace-gas census is capacity-only) and on the Discussion's desert-comparison paragraph. **Cannot establish from the abstract whether it reports MAGs or includes any Arabian site** — that is the blocking question |

**Route note (2026-09-22):** the Europe PMC render URL `europepmc.org/articles/PMC<id>?pdf=render` now
returns a **Cloudflare challenge (HTTP 403, `cf-mitigated: challenge`)** from this node, with or without a
browser user agent — the August route is dead. Two replacements work: **`nature.com` open-access PDFs**
(`https://www.nature.com/articles/<article-id>.pdf`) download cleanly with a browser user agent, and the
**Europe PMC REST full text** (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC<id>/fullTextXML`)
returns complete JATS XML without any challenge. **bioRxiv now rate-limits this node (HTTP 429)** on both
`*.full.pdf` and `*.full-text`.

## Added 2026-09-22 — Thermoplasmatota brief (L2 pass)

| paper | identifier | file into | why it is wanted |
|---|---|---|---|
| Eloe-Fadrosh, Ivanova, Woyke & Kyrpides 2016, "Metagenomics uncovers gaps in amplicon-based detection of microbial diversity", *Nat Microbiol* 1:15032 | doi:10.1038/nmicrobiol.2015.32 (bib `eloefadrosh2016gaps`) | `reference-bias/` | the citation for main.tex line 56, "amplicons cannot reach … the archaea": ~10% of environmental SSU sequences are missed by classical PCR surveys, mostly CPR and uncharacterised Archaea, measured across >6,000 assembled metagenomes. **Paywalled, no PMCID, nature.com PDF route returns HTML.** Verified from Crossref plus the Europe PMC abstract only; the "mostly Archaea" split cannot be quantified further without the paper |
| Bates, Berg-Lyons, Caporaso, Walters, Knight & Fierer 2011, "Examining the global distribution of dominant archaeal populations in soil", *ISME J* 5:908–917 | doi:10.1038/ismej.2010.171 (bib `bates2011soilarchaea`) | `dryland-ecology/` | the baseline R4 is measured against: two Nitrososphaeria phylotypes are >70% of soil archaeal sequences across 146 global soils, archaea average 2% of the community. Needed if R4 or the Discussion states that Thermoplasmatota being the cosmopolitan archaea of this sand sea is atypical for soil. **Paywalled, no PMCID**; numbers above come from the abstract, the per-biome breakdown (are any drylands in the 146?) is the blocking question |
| Maeke et al. 2025 supplementary Fig. S16/S17 (AGNOSTOS unknown fractions per MAG) | part of doi:10.1186/s40168-025-02140-8 MOESM1 (`.docx`, 10 MB, downloaded to `/tmp` this session but **not filed**) | `thermoplasmatota/` | Table S16 was mined and is recorded in READING_NOTES; the AGNOSTOS panels S17 are the only figure in that paper directly comparable to our dark-family axis. Retrieve with the same Springer MediaObjects URL if the comparison is ever put in the text |
