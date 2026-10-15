# Citation audit — empty-quarter-mags

Session 8 ("citations-sweep"), 2026-08-18. Scope: every citation in `paper/main.tex` — verify the paper
is **real** and that it **says what we cite it for** — plus resolution of every `\todo{cite}`.

**Hard gate applied:** no bibkey entered `main.tex` unless verified against a PDF on disk or a fetched
open-access source this session. Unverifiable → the `\todo` stays and the paper is listed in
`TO_DOWNLOAD.md`.

## How to read this file

| tier | what was actually checked |
|---|---|
| **T1** | PDF on disk, relevant passage read |
| **T2** | OA fulltext or abstract fetched this session |
| **T3** | Crossref/OpenAlex metadata match only — **sufficient for tool/method cites**, where the claim is "we used X" and there is no scientific assertion to misread |
| **T4** | unverifiable from this node → `\todo` retained + `TO_DOWNLOAD.md` row |

| verdict | meaning |
|---|---|
| SUPPORTS | the paper demonstrably makes/supports the claim |
| WEAK | topically right, does not establish the claim *as worded* |
| MISCITE | wrong paper for the claim |
| STALE | cites a tool/step that was not actually used |

**The papers you still need to read yourself** = the claim-bearing rows below T1. They are collected in
§5 so you do not have to scan the whole table.

---

## 1. Existence layer — 80/80 verified, 0 errors

Every distinct bibkey cited in `paper/main.tex` was resolved against the Crossref API by DOI. **All 80
DOIs resolve, and all 80 titles match the bib entry.** No fabricated, hallucinated or wrong-paper
references exist in the bibliography.

Metadata discrepancies flagged by the automated comparison were each inspected and are **all benign**:

- **14 year mismatches** — every one is Crossref's *online-first* date versus the *issue* year
  (NAR database issues, Nature online-first). Confirmed against `published-print` and `journal-issue`
  for all 14: the bib carries the correct citable issue year. No change needed.
- **4 first-author mismatches** — `rodriguezdelrio2024fesnov` and `sondergaard2016hyddb` are LaTeX
  accent escapes (`{\'\i}`, `{\o}`) versus Unicode; `bowers2017mimag` and `nei1986dnds` have an empty
  author family name in Crossref's own record. Not errors.
- **1 journal mismatch** — `edgell2006arabian` is a book, so Crossref has no container title.

| bibkey | year | venue | DOI | uses | lines |
|---|---|---|---|---|---|
| `anderson2001permanova` | 2001 | Austral Ecology | 10.1111/j.1442-9993.2001.01070.pp.x | 1 | 972 |
| `aroney2025coverm` | 2025 | Bioinformatics | 10.1093/bioinformatics/btaf147 | 5 | 392,678,793,864,+ |
| `barriohernandez2023clustering` | 2023 | Nature | 10.1038/s41586-023-06510-w | 1 | 77 |
| `bay2021trace` | 2021 | Nature Microbiology | 10.1038/s41564-020-00811-w | 2 | 174,252 |
| `berman2000pdb` | 2000 | Nucleic Acids Research | 10.1093/nar/28.1.235 | 3 | 328,722,1046 |
| `blin2023antismash` | 2023 | Nucleic Acids Research | 10.1093/nar/gkad344 | 4 | 188,274,918,1052 |
| `borcard1992varpart` | 1992 | Ecology | 10.2307/1940179 | 1 | 974 |
| `bowers2017mimag` | 2017 | Nature Biotechnology | 10.1038/nbt.3893 | 2 | 575,1026 |
| `buchfink2021diamond` | 2021 | Nature Methods | 10.1038/s41592-021-01101-x | 2 | 654,927 |
| `cantalapiedra2021eggnog` | 2021 | Molecular Biology and Evolution | 10.1093/molbev/msab293 | 2 | 653,1042 |
| `carini2016relicdna` | 2016 | Nature Microbiology | 10.1038/nmicrobiol.2016.242 | 1 | 158 |
| `chaumeil2022gtdbtk` | 2022 | Bioinformatics | 10.1093/bioinformatics/btac672 | 5 | 81,134,222,596,+ |
| `chklovski2023checkm2` | 2023 | Nature Methods | 10.1038/s41592-023-01940-w | 2 | 574,1026 |
| `coelho2022gmgc` | 2021 | Nature | 10.1038/s41586-021-04233-4 | 2 | 77,92 |
| `constant2010streptomycetes` | 2010 | Environmental Microbiology | 10.1111/j.1462-2920.2009.02130.x | 1 | 172 |
| `cordero2019atmospheric` | 2019 | The ISME Journal | 10.1038/s41396-019-0479-8 | 2 | 174,858 |
| `davis2018decontam` | 2018 | Microbiome | 10.1186/s40168-018-0605-2 | 1 | 1362 |
| `douglas2020picrust2` | 2020 | Nature Biotechnology | 10.1038/s41587-020-0548-6 | 5 | 164,438,835,1053,+ |
| `durairaj2023folds` | 2023 | Nature | 10.1038/s41586-023-06622-3 | 6 | 75,83,89,112,+ |
| `eddy2011hmmer` | 2011 | PLoS Computational Biology | 10.1371/journal.pcbi.1002195 | 1 | 1051 |
| `edgell2006arabian` | 2006 |  | 10.1007/1-4020-3970-0 | 1 | 79 |
| `eisenhofer2019lowbiomass` | 2019 | Trends in Microbiology | 10.1016/j.tim.2018.11.003 | 1 | 1362 |
| `feldgarden2021amrfinderplus` | 2021 | Scientific Reports | 10.1038/s41598-021-91456-0 | 3 | 190,918,1053 |
| `fierer2017soil` | 2017 | Nature Reviews Microbiology | 10.1038/nrmicro.2017.87 | 3 | 79,96,463 |
| `greening2015persistence` | 2015 | Proceedings of the National Academ | 10.1073/pnas.1508385112 | 2 | 174,176 |
| `greening2016genomic` | 2015 | The ISME Journal | 10.1038/ismej.2015.153 | 3 | 172,252,856 |
| `greening2022review` | 2022 | Nature Reviews Microbiology | 10.1038/s41579-022-00724-x | 1 | 176 |
| `grinter2023structural` | 2023 | Nature | 10.1038/s41586-023-05781-7 | 1 | 172 |
| `heinzinger2024prostt5` | 2024 | NAR Genomics and Bioinformatics | 10.1093/nargab/lqae150 | 3 | 323,715,1045 |
| `huang2016dryland` | 2015 | Nature Climate Change | 10.1038/nclimate2837 | 2 | 79,94 |
| `hyatt2010prodigal` | 2010 | BMC Bioinformatics | 10.1186/1471-2105-11-119 | 1 | 623 |
| `jain2018fastani` | 2018 | Nature Communications | 10.1038/s41467-018-07641-9 | 2 | 585,1026 |
| `ji2017atmospheric` | 2017 | Nature | 10.1038/nature25014 | 5 | 79,97,176,252,+ |
| `kautsar2021bigfam` | 2020 | Nucleic Acids Research | 10.1093/nar/gkaa812 | 2 | 274,1052 |
| `khan2020aridsoil` | 2020 | SOIL | 10.5194/soil-6-513-2020 | 1 | 79 |
| `kim2026binning` | 2026 | Nature Communications | 10.1038/s41467-026-71521-w | 1 | 242 |
| `king2007coxidizers` | 2007 | Nature Reviews Microbiology | 10.1038/nrmicro1595 | 1 | 174 |
| `kryazhimskiy2008dnds` | 2008 | PLoS Genetics | 10.1371/journal.pgen.1000304 | 2 | 382,1050 |
| `larralde2022pyrodigal` | 2022 | Journal of Open Source Software | 10.21105/joss.04296 | 2 | 361,623 |
| `lin2023esmfold` | 2023 | Science | 10.1126/science.ade2574 | 3 | 323,731,1046 |
| `liu2022vfdb` | 2021 | Nucleic Acids Research | 10.1093/nar/gkab1107 | 3 | 192,920,1053 |
| `maestre2015aridity` | 2015 | Proceedings of the National Academ | 10.1073/pnas.1516684112 | 2 | 79,96 |
| `makhalanyane2015microbial` | 2015 | FEMS Microbiology Reviews | 10.1093/femsre/fuu011 | 2 | 79,96 |
| `matchado2024limits` | 2024 | Microbial Genomics | 10.1099/mgen.0.001203 | 5 | 83,115,168,438,+ |
| `mcinerney2017openpangenome` | 2017 | Nature Microbiology | 10.1038/nmicrobiol.2017.40 | 1 | 81 |
| `murrell2013fubar` | 2013 | Molecular Biology and Evolution | 10.1093/molbev/mst030 | 2 | 407,1050 |
| `nayfach2021gem` | 2020 | Nature Biotechnology | 10.1038/s41587-020-0718-6 | 6 | 77,206,456,931,+ |
| `nei1986dnds` | 1986 | Molecular Biology and Evolution | 10.1093/oxfordjournals.molbev.a040410 | 4 | 156,382,783,1050 |
| `olm2017drep` | 2017 | The ISME Journal | 10.1038/ismej.2017.126 | 2 | 584,1026 |
| `orakov2021gunc` | 2021 | Genome Biology | 10.1186/s13059-021-02393-0 | 1 | 1026 |
| `pan2023semibin2` | 2023 | Bioinformatics | 10.1093/bioinformatics/btad209 | 2 | 557,1026 |
| `parks2022gtdb` | 2021 | Nucleic Acids Research | 10.1093/nar/gkab776 | 5 | 81,134,222,597,+ |
| `pavlopoulos2023darkmatter` | 2023 | Nature | 10.1038/s41586-023-06583-7 | 2 | 75,89 |
| `ray2022atmospheric` | 2022 | The ISME Journal | 10.1038/s41396-022-01298-5 | 2 | 79,97 |
| `richardson2023mgnify` | 2022 | Nucleic Acids Research | 10.1093/nar/gkac1080 | 8 | 77,92,206,210,+ |
| `rodriguezdelrio2024fesnov` | 2023 | Nature | 10.1038/s41586-023-06955-z | 5 | 83,89,110,156,+ |
| `salter2014contamination` | 2014 | BMC Biology | 10.1186/s12915-014-0087-z | 2 | 158,1362 |
| `sberro2019smorfs` | 2019 | Cell | 10.1016/j.cell.2019.07.016 | 1 | 317 |
| `schmidt2024spire` | 2023 | Nucleic Acids Research | 10.1093/nar/gkad943 | 3 | 206,932,1339 |
| `schwengers2021bakta` | 2021 | Microbial Genomics | 10.1099/mgen.0.000685 | 2 | 652,1042 |
| `sieber2018dastool` | 2018 | Nature Microbiology | 10.1038/s41564-018-0171-1 | 2 | 571,1026 |
| `sondergaard2016hyddb` | 2016 | Scientific Reports | 10.1038/srep34212 | 5 | 172,252,856,1051,+ |
| `steen2019uncultured` | 2019 | The ISME Journal | 10.1038/s41396-019-0484-y | 2 | 75,89 |
| `steinegger2017mmseqs2` | 2017 | Nature Biotechnology | 10.1038/nbt.3988 | 2 | 629,1042 |
| `steinegger2018linclust` | 2018 | Nature Communications | 10.1038/s41467-018-04964-5 | 3 | 296,361,629 |
| `sun2020inference` | 2020 | Microbiome | 10.1186/s40168-020-00815-y | 2 | 83,115 |
| `suzek2015uniref` | 2014 | Bioinformatics | 10.1093/bioinformatics/btu739 | 4 | 298,377,632,1044 |
| `tettelin2005pangenome` | 2005 | Proceedings of the National Academ | 10.1073/pnas.0506758102 | 1 | 81 |
| `uritskiy2019atacama` | 2019 | The ISME Journal | 10.1038/s41396-019-0468-y | 2 | 210,456 |
| `vankempen2024foldseek` | 2023 | Nature Biotechnology | 10.1038/s41587-023-01773-0 | 3 | 323,718,1045 |
| `vanni2022agnostos` | 2022 | eLife | 10.7554/eLife.67667 | 2 | 75,89 |
| `varadi2022alphafolddb` | 2021 | Nucleic Acids Research | 10.1093/nar/gkab1061 | 1 | 323 |
| `varadi2024alphafolddb` | 2023 | Nucleic Acids Research | 10.1093/nar/gkad1011 | 4 | 328,351,721,1046 |
| `wang2024comebin` | 2024 | Nature Communications | 10.1038/s41467-023-44290-z | 1 | 242 |
| `wilson1927` | 1927 | Journal of the American Statistica | 10.1080/01621459.1927.10502953 | 1 | 342 |
| `xu2010tmscore` | 2010 | Bioinformatics | 10.1093/bioinformatics/btq066 | 1 | 331 |
| `yeo2025afesm` | 2025 | bioRxiv | 10.1101/2025.04.23.650224 | 3 | 83,112,356 |
| `zelezniak2015smetana` | 2015 | Proceedings of the National Academ | 10.1073/pnas.1421834112 | 2 | 444,1054 |
| `zhang2005tmalign` | 2005 | Nucleic Acids Research | 10.1093/nar/gki524 | 1 | 331 |
| `zimmermann2021gapseq` | 2021 | Genome Biology | 10.1186/s13059-021-02295-1 | 2 | 444,1054 |
---

## 2. Tool / method layer (T3) — correct, with one version miscite

Every Methods citation was checked for **key ↔ tool correspondence**: does the cited paper actually
describe the tool named in that sentence? All 46 tool-cite instances resolve to the right software
paper. **One is wrong** (the antiSMASH version — the fix applied is recorded in §4b), and one initially
looked wrong but is correct:

| # | line | issue | verdict | action |
|---|---|---|---|---|
| T-1 | `paper/main.tex:918` (M.11) and `:188` (§2.4) | The text states **antiSMASH v8.0.4** but cites `blin2023antismash`, which is the **antiSMASH 7.0** paper (NAR 2023). The v8 paper exists and was verified: Blin et al. 2025, *NAR* 53:W32–W38, doi 10.1093/nar/gkaf334. | MISCITE (version) | add `blin2025antismash8` and cite it where this paper's v8.0.4 run is described |

### Checked and found NOT to be a defect

`orakov2021gunc` at `main.tex:1026` first read as a stale cite, because session 7a established on disk
that **GUNC never ran**. It is not a defect. That citation sits inside the `% PI pointer-level Methods`
block, which carries an explicit `\status{HISTORICAL; … superseded by M.1--M.14. Where they disagree ---
GUNC, the 573-sample basis, HyPhy-FUBAR, gapseq/SMETANA and the PICRUSt2 gating --- the numbered
subsections govern.}` marker at line 1021. The block is retained by design (session 7b) as the frozen
drafting source, GUNC is *already named* there as a known disagreement, and M.2 states plainly that no
dedicated chimerism screen was run. The citation correctly identifies the GUNC paper for the word "GUNC"
in a passage the manuscript already flags as superseded. **Do not remove it** — doing so would damage a
deliberately preserved historical record.

`blin2023antismash` also remains correct at `:1052` (same historical block) and `:274` (PI scaffold),
where it refers to antiSMASH generally rather than to the v8.0.4 run.

---

## 3. `\todo{cite}` resolution

28 `\todo{...}` blocks exist in `main.tex`. **17 are citation todos** (this session's scope); 11 are not
(numeric placeholders for §1.3, session-9 figure notes, and two items blocked on Alejandra) and are
listed in §6 so they are visibly *not* forgotten.

### 3a. New bib entries — 21 added, every one fetched and verified this session

Generated from authoritative Crossref records (content negotiation / REST), **not written from memory**;
author lists truncated to the repo's existing 10-plus-`others` convention.

| bibkey | paper | venue | DOI | verified via |
|---|---|---|---|---|
| `chen2018fastp` | fastp: an ultra-fast all-in-one FASTQ preprocessor | Bioinformatics 34:i884 | 10.1093/bioinformatics/bty560 | Crossref |
| `li2015megahit` | MEGAHIT | Bioinformatics 31:1674 | 10.1093/bioinformatics/btv033 | Crossref |
| `gurevich2013quast` | QUAST | Bioinformatics 29:1072 | 10.1093/bioinformatics/btt086 | Crossref |
| `li2018minimap2` | Minimap2 | Bioinformatics 34:3094 | 10.1093/bioinformatics/bty191 | Crossref |
| `li2013bwamem` | Aligning sequence reads … with BWA-MEM | arXiv:1303.3997 | 10.48550/arXiv.1303.3997 | OpenAlex |
| `molder2021snakemake` | Sustainable data analysis with Snakemake | F1000Research 10:33 | 10.12688/f1000research.29032.2 | Crossref |
| `harris2020numpy` | Array programming with NumPy | Nature 585:357 | 10.1038/s41586-020-2649-2 | Crossref |
| `virtanen2020scipy` | SciPy 1.0 | Nat Methods 17:261 | 10.1038/s41592-019-0686-2 | Crossref |
| `cock2009biopython` | Biopython | Bioinformatics 25:1422 | 10.1093/bioinformatics/btp163 | Crossref |
| `benjamini1995fdr` | Controlling the False Discovery Rate | JRSS B 57:289 | 10.1111/j.2517-6161.1995.tb02031.x | Crossref |
| `holm1979multiple` | A Simple Sequentially Rejective Multiple Test Procedure | Scand J Stat 6:65 | 10.2307/4615733 | OpenAlex (not in Crossref) |
| `mantel1967` | The Detection of Disease Clustering … | Cancer Res 27:209 | — (PMID 6018555) | OpenAlex; **no DOI exists** |
| `zhang2004tmscore` | Scoring function for … template quality | Proteins 57:702 | 10.1002/prot.20264 | Crossref |
| `mcardle2001dbrda` | Fitting Multivariate Models to Community Data | Ecology 82:290 | 10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2 | Crossref |
| `dcosta2011ancient` | Antibiotic resistance is ancient | Nature 477:457 | 10.1038/nature10388 | Crossref |
| `vangoethem2018antarctic` | 'Historical' ARGs in pristine Antarctic soils | Microbiome 6:40 | 10.1186/s40168-018-0424-5 | Crossref |
| `lennon2018relicdna` | How, When, and Where Relic DNA Affects Microbial Diversity | mBio 9:e00637-18 | 10.1128/mbio.00637-18 | Crossref |
| `bay2021aridity` | Chemosynthetic and photosynthetic bacteria … desert aridity gradient | ISME J 15:3339 | 10.1038/s41396-021-01001-0 | Crossref |
| `blin2025antismash8` | antiSMASH 8.0 | NAR 53:W32 | 10.1093/nar/gkaf334 | Crossref |
| `jordaan2020gobi` | Hydrogen-oxidizing bacteria are abundant in desert soils | mSystems 5:e01131-20 | 10.1128/msystems.01131-20 | Crossref (PDF on disk) |
| `aplakidou2025quadrupling` | Quadrupling the protein family space with global metagenomics | bioRxiv (posted 2025-12-23) | 10.64898/2025.12.21.695754 | Crossref + bioRxiv fulltext |

**Note on `mantel1967`:** OpenAlex returns a DOI for this record that belongs to a *Journal of
Biogeography* article — an upstream metadata error. The entry therefore carries the PMID and no DOI.
Do not "fix" it by pasting that DOI in.

**Two names from session 7b's owed list are not needed:** `seqkit` and `samtools` are never mentioned in
`main.tex`, so no bibkey is owed for either. BWA-MEM *is* mentioned (twice) and is covered above.

### 3b. `\todo{cite}` markers resolved — 17 of 17

| line | what it asked for | resolved with | tier |
|---|---|---|---|
| 134 | expectation of novel deep lineages in extreme/undersampled environments | `liu2025qaidam` (Qaidam Mars-analog desert, 94.5% of 1,773 MAGs novel taxa) | T1 |
| 134 | MIMAG, Bowers 2017 | `bowers2017mimag` (already in bib) | T3 |
| 146 | ProstT5, FoldSeek, AFDB, PDB | `heinzinger2024prostt5`, `vankempen2024foldseek`, `varadi2024alphafolddb`, `berman2000pdb` | T3 |
| 148 | ESMFold | `lin2023esmfold` | T3 |
| 158 | Lennon 2018 relic DNA | `lennon2018relicdna` (new) | T2 |
| 176 | Ehhalt & Rohrer 2009 global soil H₂ sink denominator | `ehhalt2009h2cycle` (already in bib) | T3 |
| 180 | desiccation–radiation cross-resistance | `mattimore1996radioresistance`, `slade2011oxidative` | T2 |
| 180 | bacterial NHEJ, Pitcher 2007 | `pitcher2007nhej` + `sharda2020nhej` for the prevalence clause | T2 / T1 |
| 184 | Oren, salt-in strategy | `oren2008halophiles` (**not** `oren1999bioenergetics` — that is a cost argument) | T2 |
| 190 | intrinsic environmental resistome | `vangoethem2018antarctic`, `dcosta2011ancient` (new, Van Goethem first — it is the load-bearing pristine-soil result) | T2 |
| 208 | Aplakidou et al. 2025 | `aplakidou2025quadrupling` (new) | T2 |
| 210 | bay2021aridity, distinct from bay2021trace | `bay2021aridity` (new) | T2 |
| 531 | fastp | `chen2018fastp` (new) | T3 |
| 554 | MEGAHIT, QUAST, minimap2, BWA-MEM, Snakemake | 5 new keys, placed at first mention of each tool | T3 |
| 749 | TM-score normalisation, Zhang & Skolnick 2004 | `zhang2004tmscore` (new) | T3 |
| 811 | McArdle & Anderson 2001 distance-based RDA | `mcardle2001dbrda` (new — **distinct** from the existing `anderson2001permanova`, which is Anderson 2001 *Austral Ecology*) | T3 |
| 975 | NumPy, SciPy, Biopython, Benjamini–Hochberg, Holm, Mantel | 6 new keys, placed at first mention | T3 |

⚠ **The `\todo` at line 208 named the wrong DOI.** It said bioRxiv `10.1101/2025.12.21.695754`; that 404s. bioRxiv
has changed its prefix and the live DOI is **`10.64898/2025.12.21.695754`** (checked: 200 vs 404). The bib entry
carries the correct one and a comment saying so.

---

## 4. Claim layer — does each citation say what we cite it for?

Verified by three literature agents against PDFs on disk, Europe PMC and publisher OA fulltext. **Nothing here
is a fabricated reference** — every paper is real (§1). These are *attachment* problems: the right paper cited
for a claim it does not actually make.

### 4a. Confirmed correct (SUPPORTS) — no action

Introduction: `steen2019uncultured`, `vanni2022agnostos` (eLife trap re-confirmed — bib is right), `durairaj2023folds`
(for "no matched structure"), `coelho2022gmgc`, `richardson2023mgnify`, `nayfach2021gem`, `barriohernandez2023clustering`,
`maestre2015aridity`, `ji2017atmospheric` + `ray2022atmospheric` (cold-desert chemosynthesis), `khan2020aridsoil`,
`chaumeil2022gtdbtk`, `parks2022gtdb`, `mcinerney2017openpangenome`, `rodriguezdelrio2024fesnov`, `yeo2025afesm`,
`sun2020inference`. §1.4: `nei1986dnds`, `rodriguezdelrio2024fesnov`, `carini2016relicdna`, `salter2014contamination`.
§2.2: `greening2015persistence`, `bay2021trace` (the ">70% = inorganic energy sources" trap is **not** triggered —
L174 makes no percentage claim). §2.3 todos: `mattimore1996radioresistance`, `slade2011oxidative`, `pitcher2007nhej`,
`sharda2020nhej`, `oren2008halophiles`. §3.1: `richardson2023mgnify`/`nayfach2021gem`/`schmidt2024spire` as comparator
identities. Box 2: `salter2014contamination`, `eisenhofer2019lowbiomass`, `davis2018decontam` are the correct canonical trio.

### 4b. Fixed this session (citations added to close an uncited claim; **no prose was rewritten**)

| line | the gap | fix applied |
|---|---|---|
| 172 + 854 | **HydDB group "1l" appears in neither cited paper.** `sondergaard2016hyddb` and `greening2016genomic` document group 1h; the published [NiFe] group 1 scheme runs 1a–1k and a regex for `1l` returns 0 hits in both PDFs. | added `ortiz2021antarctic` (names and characterises group 1l / HylL) at both the Results claim and the M.9 definition, and `jordaan2020gobi` (Gobi desert: "almost all community members … encode the group 1l") |
| 172 | the "…rather than the heterotrophic **carbon-monoxide-dehydrogenase** groups" half was uncited — neither hydrogenase paper mentions CODH | added `cordero2019atmospheric`, which draws exactly that form-I/form-II line and is already used for the same purpose in M.9 |
| 918 + 188 | antiSMASH **v8.0.4** cited to the antiSMASH **7.0** paper | `blin2025antismash8` |

### 4c. ⚠ DECISIONS FOR ALEJANDRA — claims whose citation does not support the wording

**9 of the 15 rows below are still open.** Each needs a wording call, and every affected
section is already shipped, so the call is yours. Ordered by severity. The evidence is quoted so you do not
have to re-open the papers.

⚠ **Exception 2 — D-9, D-10, D-11 and D-13 are MOOT as of the 2026-08-20 Introduction rewrite** (verified
2026-08-23, gap-review). All four quote Introduction sentences that no longer exist in `paper/main.tex`:
none of "genuinely new folds are", "no experimentally characterised function", "not fixed by phylogeny",
"actinobacterial diversity" or "close to saturation" survives, and `tettelin2005pangenome` (the D-11
anchor) is now cited 0 times anywhere in the file. They are marked ⬜ MOOT below. Re-open only if the
rewritten Introduction is reverted. **Live rows: D-1…D-6, D-12, D-14, D-15, plus D-16 in §9.**

⚠ **Exception — D-7 and D-8 were APPLIED by a parallel session on 2026-08-18, not by the session that wrote
this audit.** Verified on disk at `paper/main.tex:79`. D-7 was fully applied by that session; **D-8 was completed on 2026-08-18** (`fierer2017soil` and
`makhalanyane2015microbial` dropped from L79). Both are marked ✅ below and are excluded from the 13.

| # | line | § | problem | options |
|---|---|---|---|---|
| **D-1** | 210 | 3.1 | **The claim is factually inverted.** "MGnify90 itself is dominated by aquatic and engineered biomes \citep{richardson2023mgnify}". Richardson 2023 p.D754 actually says: "297 different biomes … with over half of MGnify's analyses originating from merely nine of them: human-faecal, -oral, -digestive system, -skin, and unspecified human; marine; soil; mammalian digestive systems; and mixed biome samples." So the dominant biomes are **human/host-associated**; **soil is itself in the top nine**; **"engineered" is never named**. The workbench notes that June's "83.5% aquatic" figure is *not recomputable* and described where *our families hit*, not what MGnify contains — this looks like that dead number re-attributed to Richardson. | (a) restate to what the paper says — "over half of MGnify's analyses come from nine biomes, led by human- and other host-associated samples, marine and unspecified soil" (still makes your point: no hot-hyperarid category); or (b) **delete the clause** — the desert-comparator sentence already carries §3.1's argument |
| **D-2** | 210 | 3.1 | "the five desert comparators are the largest desert **soil** collections available … and **no Arabian or Saharan site appears in any of them**" cites only `uritskiy2019atacama`, which is (i) one of five, (ii) *halite salt-rock endoliths*, not soil, and (iii) cannot support a negative. `bay2021aridity` is now added, but three legs (arid, antarctic_colddesert, mojave) are unpublished JGI studies with no paper at all. | drop "soil" → "desert metagenome collections"; attribute the set as "five desert metagenome studies assembled in SPIRE \citep{schmidt2024spire}" with per-study accessions in M.12; rest the negative on the M.12 accession list, not on literature. Consider "no site on the Arabian Peninsula and none in the Sahara" — the Negev is the nearest of the five |
| **D-3** | 176 | 2.2 | "hot hyperarid drylands … are **effectively absent from those budgets**" cites `greening2022review`, whose abstract says nothing about biome coverage of global trace-gas budgets. A strong negative on an abstract-only source. | soften to "poorly constrained in those budgets", and/or cite sources that actually address it — both on disk, neither keyed: Reji et al., "Global implications of a low soil moisture threshold for microbial hydrogen uptake" (*Nat Commun*, 10.1038/s41467-025-67208-3) — "the soil microbial H₂ sink — the largest and most uncertain term in the global H₂ budget"; and Ray et al. 2023 *MMBR* 87:e0004823 (10.1128/mmbr.00048-23) — observations "restricted to cold deserts", calls for work "extending to hot desert environments" |
| **D-4** | 174 | 2.2 | "consistent with aerobic CO oxidation being **the more broadly distributed trait**" — neither `cordero2019atmospheric` nor `king2007coxidizers` makes a CO-vs-H₂ comparison (Cordero says "more widespread **than previously reported**"). Two desert MAG surveys on disk run the *opposite* way: Hutchinson 2025 Antarctic = 322/554 MAGs with 1h/1l vs 60/554 with coxL. | make it local: "…in this catalogue CO capacity is slightly the more prevalent of the two, consistent with aerobic CO oxidation being taxonomically widespread in soils (685 species across seven dominant soil phyla)" |
| **D-5** | 176 | 2.2 | `ji2017atmospheric` cited for "maintenance and persistence **rather than growth**", but Ji's title is "Atmospheric trace gases **support primary production**" — the cite is attached to the opposite of its headline. Ji *does* support the sufficiency half ("rates sufficient to sustain their theoretical maintenance energy"). | either reword to the sufficiency claim, or let `greening2015persistence` + `bay2021trace` carry "maintenance rather than growth" (Bay's thermodynamic modelling states it without the primary-production framing) |
| **D-6** | 172 | 2.2 | "the classic **mesic-soil** group 1h" — no source on disk calls group 1h "mesic". Constant 2010 isolated its streptomycetes from "temperate forest and agricultural soils". | "the classic group 1h described from temperate forest and agricultural soils" |
| ✅ **D-7** *(APPLIED)* | 79 | Intro | **Done — `main.tex:79` now reads `\citep{huang2016dryland,maestre2015aridity}`.** Original finding: "Drylands are **the largest terrestrial biome** and are expanding under climate change \citep{huang2016dryland}". Huang supports *expansion* (title: "Accelerated dryland expansion under climate change"); no evidence found that it states "largest terrestrial biome". | keep Huang for the expansion clause; anchor the size claim on `maestre2015aridity`, already in bib and cited, which states drylands "occupy ∼41% of Earth's surface" (no new key needed). Alternative on disk: Coleine et al. 2024 *ISME J* wrae056 (10.1093/ismejo/wrae056) — "Drylands account for 45% of the Earth's land area" |
| ✅ **D-8** *(APPLIED)* | 79 | Intro | **Done 2026-08-18 — cut to the two keys that carry the claim: `\citep{maestre2015aridity,delgadobaquerizo2013decoupling}`.** `fierer2017soil` and `makhalanyane2015microbial` were dropped here (both WEAK for this sentence); `fierer2017soil` remains cited at L96 and L463, so nothing is orphaned. Original finding: "aridity … lowers microbial diversity and reorganises soil communities **and nutrient cycling**" cites `maestre2015aridity` (SUPPORTS), `makhalanyane2015microbial` (a *review*, WEAK) and `fierer2017soil` (whose abstract has **nothing** on aridity, drylands or nutrient cycling). The nutrient-cycling half is unsupported. | `delgadobaquerizo2013decoupling` — "Decoupling of soil nutrient cycles as a function of aridity in global drylands" — is **already in references.bib** and is the exact anchor. Suggested: `\citep{maestre2015aridity,delgadobaquerizo2013decoupling}`; move or drop `fierer2017soil` here |
| ⬜ **D-9** *(MOOT — Intro rewritten 2026-08-20)* | 83 | Intro | "global structural surveys indicate that genuinely new folds are **rare** … close to saturation" cites `yeo2025afesm` (SUPPORTS precisely) **and `durairaj2023folds`, which is a fold-*discovery* paper** making no saturation or rarity claim (grep for `saturat`/`rarely`/`few novel` → no hits). | drop `durairaj2023folds` from this sentence; `yeo2025afesm` carries it alone. Do **not** substitute `lau2024ted` — TED reports 7,427 putative novel folds, arguing the opposite |
| ⬜ **D-10** *(MOOT — Intro rewritten 2026-08-20)* | 75 | Intro | "most of the proteins they encode have **no experimentally characterised function**" cites `pavlopoulos2023darkmatter`, which measures *reference-invisibility* ("1.17 billion protein sequences … with no similarity to any sequences from 102,491 reference genomes or the Pfam database"), not experimental characterisation. | reword to "…fall outside all reference families and have no assigned function". Keeping "experimentally characterised" needs a Swiss-Prot/reviewed-fraction statistic, which is not in the bib |
| ⬜ **D-11** *(MOOT — Intro rewritten 2026-08-20)* | 81 | Intro | "accessory content is **not fixed by phylogeny**" cites `tettelin2005pangenome`, which establishes *openness* from six *S. agalactiae* genomes only. | `whelan2021accessory` is on disk and in the bib (86.7% of common accessory genes in coincident relationships) — the right anchor for the not-fixed-by-phylogeny half |
| **D-12** | 168 / 466 | 2.1 / Disc | `matchado2024limits` carries "cannot reach the majority" and "affects soil-microbiome inference broadly", but its data are **human gut (CRC/obesity/T2D) + simulation — no soil**. Its only soil statement is second-hand. | pair it with `sun2020inference`, already in bib and on disk, the one source with **soil** data (for both soil datasets, unpermuted correlations were not significantly different from permuted). Suggested `\citep{sun2020inference,matchado2024limits}` at both lines |
| ⬜ **D-13** *(MOOT — Intro rewritten 2026-08-20)* | 79 | Intro | "Genome-resolved desert microbiology has concentrated on … notably **actinobacterial diversity** and atmospheric chemosynthesis" — the actinobacterial clause is carried by neither `ji2017atmospheric` nor `ray2022atmospheric` (Ji's dominant taxa are WPS-2/AD3). | in-bib, Crossref-verified options: `bull2016atacama`, `idris2017atacama`, `goodfellow2018raretaxa` — but all three are culture/amplicon work, so citing them means dropping "genome-resolved" from that clause. The genome-resolved option is `liu2025qaidam` (565 of 1,773 MAGs Actinomycetota) |
| **D-14** | 1362 | Box 2 | The contamination trio is on the **"Reference bias overstated"** alternative, while the **"Signal = contamination/kit-ome"** alternative one clause earlier carries no citation. | move `\citep{salter2014contamination,eisenhofer2019lowbiomass,davis2018decontam}` to the kit-ome clause it actually supports |
| **D-15** | 210 | 3.1 | Scope note on the Aplakidou concession (not an error): its "ecosystem" is GOLD's 3 top-level / 8 sub-type classification, coarser than "biome" as §3.1 uses it elsewhere; and the 608,258 families were *defined* as those with ≥100 members and ≥5% of samples from one ecosystem, so a mild restriction bias is built into the denominator. Neither undermines the concession — both make it safer — but a reviewer who opens the preprint will see them. | optional one-clause hedge in M.12 |

### 4d. Bib metadata repaired

- **`grinter2023structural`** — byline was **12 of 19 authors**; Crossref (on the entry's own DOI) confirms 19. Added Barlow, Belousoff, Shafaat, Schittenhelm, Vincent, Khalid, Berggren; cleared the stale `note = {verify ordering against PDF}`.
- **`uritskiy2019atacama`** — given names wrong: **Getsin, Samuel → Samantha**; **Glass, Benjamin → Brian**. (Crossref-confirmed. Suggests a reconstructed byline; the rest of the entry is correct.)

### 4e. Reported defect that did NOT reproduce

An agent reported that `schulzemakuch2018transitory` carries Lebre 2017's journal/volume/pages/DOI. **This is not
true of the current file.** The entry holds *PNAS* 115:2670–2675, doi `10.1073/pnas.1714341115`, which Crossref
resolves to "Transitory microbial habitat in the hyperarid Atacama Desert". No change made. Recorded here so the
claim is not re-litigated.

---

## 5. Papers you still need to read yourself

These are the **claim-bearing** citations that could not be verified at T1 (PDF in hand). Everything else in the
manuscript is either a tool cite (T3 is sufficient) or was read this session. Rows added to `TO_DOWNLOAD.md`.

| paper | why it matters | status from this node |
|---|---|---|
| **Edgell 2006, *Arabian Deserts* (Springer)** | carries "the world's largest contiguous sand sea" — the opening geographic claim of the paper | **T4 UNVERIFIABLE.** Book confirmed (10.1007/1-4020-3970-0); Springer fulltext redirects to `idp.springer.com` auth. Cite the chapter so a reader can check: ch. 7 "Arabian sand seas", pp. 125–166 (10.1007/1-4020-3970-0_7). The wording needs confirming off-node |
| **Huang et al. 2016, *Nat Clim Change* 6:166** | the "largest terrestrial biome" claim (D-7) | T3 only — nature.com blocked, OpenAlex abstract null, S2 abstract elided by publisher. Metadata confirmed, content not |
| **Ji et al. 2017, *Nature* 552:400** | load-bearing for §2.2's energy argument (D-5) | T2 abstract only; no PDF on disk. `TO_DOWNLOAD.md` row still valid |
| **Greening et al. 2022, *Nat Rev Microbiol* 20:513** | carries the "absent from budgets" negative (D-3) | T2 abstract only; paywalled, no PMCID |
| **D'Costa et al. 2011, *Nature* 477:457** | the pre-anthropogenic half of §2.4's resistome claim | T2 abstract only (not OA). Sufficient for the cite-only role it now has |
| **King & Weber 2007, *Nat Rev Microbiol* 5:107** | the CO-oxidiser breadth claim (D-4) | T2 abstract only; paywalled review |
| **Constant et al. 2010, *Environ Microbiol* 12:821** | "classic mesic-soil group 1h" (D-6) | T2 abstract only. ⚠ the on-disk `Constant_2016_streptomyces_h2.pdf` is **Liot & Constant 2016**, a *different* paper — do not mistake it for this one |
| **Parks et al. 2017, *Nat Microbiol* 2:1533** | optional second anchor for the L134 deep-lineage expectation | T2 abstract only, no PMCID. Not added — `liu2025qaidam` (T1, on disk) carries the claim |

Two useful routes, both confirmed working from this node and worth remembering: **Europe PMC render**
(`https://europepmc.org/articles/PMC<ID>?pdf=render`) reaches OA PMC PDFs even though `ncbi.nlm.nih.gov/pmc` is
blocked; and **OpenAlex** (`api.openalex.org`) resolves records Crossref does not — it supplied Holm 1979,
Mantel 1967 and the BWA-MEM preprint. arXiv's export API and Semantic Scholar were **not** reachable (empty
response / HTTP 429).

---

## 6. The 11 `\todo` markers deliberately left in place

Seen, not forgotten, not this session's scope:

| count | what | owner |
|---|---|---|
| 6 | `\todo{N}`, `{P}`, `{C}`, `{K}`, `{Q}`, `{R}` — §1.3 numeric placeholders | session 6, when the structural chain lands |
| 1 | "finalise magnitude + framing once 03_search + ESMFold-validation numbers land" | session 6 |
| 1 | "§2.2: promote to a 5,438 guild table and a HydDB/CODH placement figure" | session 9 (figures/tables) |
| 1 | companion ecology manuscript (Tawfiq et al.) — **unverifiable by definition**, it is unpublished | Alejandra: decide whether to cite as "in preparation" or drop |
| 1 | M.10 XRF workflow confirmation | blocked on Alejandra |
| 1 | M.14 ENA accessions / MAG deposition / archived-release DOI | blocked on Alejandra, before submission |

---

## 7. Reading records owed to `READING_NOTES.md`

Ten sources were read during this audit that have no entry yet. Worth distilling, especially the first two:
**Richardson 2023** (the biome-composition passage — it *disproves* a current claim, D-1), **Durairaj 2023**
(34% of UniRef50 clusters dark; 78% have an AFDB model; β-flower the only new fold; **explicitly no saturation
claim**), Aplakidou 2025, Bay 2021 aridity, Uritskiy 2019, Van Goethem 2018, D'Costa 2011, Pavlopoulos 2023
(scope = reference-invisibility), Coelho 2022 GMGC (303M genes, most habitat-specific — a §3.1 asset, not just
an Introduction name), Sun 2020 (the only soil-data source in the reference-bias bucket).

---

## 8. Summary

| | |
|---|---|
| distinct bibkeys cited in `main.tex` | 80 → **110** after this session (30 newly cited: the 21 new entries plus 9 already in the bib but never cited; none dropped) |
| existence verified | **80/80, zero fabricated or wrong-paper references** |
| entries in `references.bib` | 165 → **186** (21 added, 0 duplicates) |
| `\todo` blocks | 28 → **11** (all 17 citation todos resolved) |
| claim-support defects found | **13 still open (§4c) + 1 in the supplement (§9); 2 applied by a parallel session (D-7 fully, D-8 half); 3 fixed by adding a citation (§4b)** — *as found on 2026-08-18; see §4c Exception 2: D-9/D-10/D-11/D-13 went MOOT with the 2026-08-20 Introduction rewrite, so the live count is 9 + D-16* |
| bib metadata repaired | 2 (`grinter2023structural` byline 12→19; `uritskiy2019atacama` given names) |
| build | 42 pp, **0 undefined citations, 0 undefined references, 0 overfull boxes** |

**Bottom line: the bibliography is sound — nothing is invented.** The real finding is that a handful of claims,
mostly in the Introduction and §2.2, are attached to papers that are topically right but do not establish the
claim as worded, plus one factually inverted statement in §3.1 (D-1). None of it was rewritten; §4c is the
decision list.

---

## 9. `paper/supplementary.tex` — scanned separately, and it needs a decision

The supplement contains **zero `\cite` commands and zero `\todo` markers**, so a grep-based citation audit
sees nothing in it. That is misleading: it carries **three plain-text author-year citations** in the
trace-gas energy calculation (lines 78–81), because the file is standalone (`pdflatex supplementary.tex`)
and has no bibliography attached.

All three were verified this session and **all three are correct**:

| in-text | bibkey (exists in `references.bib`) | verified |
|---|---|---|
| "semi-arid chamber fluxes of Buzzard et al. 2022" | `buzzard2022h2uptake` — *PeerJ* 2022, "Sensitivity of soil hydrogen uptake to natural and managed moisture dynamics in a semiarid urban ecosystem" | Crossref |
| "global soil H₂ sink (60–80 Tg yr⁻¹; Ehhalt & Rohrer 2009)" | `ehhalt2009h2cycle` — *Tellus B* 61:500, "The tropospheric cycle of H₂: a critical review" | Crossref |
| "indirect GWP₁₀₀ ≈ 11 of H₂ (Ocko & Hamburg 2022)" | `ocko2022hydrogen` — *Atmos Chem Phys* 2022, "Climate consequences of hydrogen emissions" | Crossref |

The file's header comment claimed these bibkeys were **missing** and had to be added before conversion.
That comment was **stale** — all four keys it named (the three above plus `ortiz2021antarctic`) already
exist. The comment has been corrected in place.

**⚠ Decision D-16 (yours):** convert the supplement's three plain-text citations to `\citep`, which
requires attaching a bibliography to the standalone file, or leave them as author-year prose. This is a
structural change to how the supplement builds, so it was **not** done here.

