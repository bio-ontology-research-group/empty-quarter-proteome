# EQ-MAGs: open items

The one open-items register (rewritten 2026-09-24 after the paper evaluation and the L1/L2 literature
passes). Structure: `PAPER_STRUCTURE.md`; numbers: `CATALOGUE_NUMBERS.md`; what changed and why:
`CORRECTIONS.md`; literature: `../resources/literature/READING_NOTES.md`. Done items are removed
here and live in CORRECTIONS.

Rule: a proof-of-concept number never enters the registry or the manuscript until it has been rerun
in a workbench `_clean/` folder and registered with its basis.

## A. Decisions for Alejandra

| # | decision | bears on |
|---|---|---|
| A1 | Archaea in the title. A lineage-led paper is the only visible route to Nature Microbiology; the current decision keeps them out | title, abstract, Discussion |
| A2 | Two reagent contaminants confirmed by the no-template blank (*E. coli* `V42PRr2__single_SemiBin_33`, *C. acnes* `51Dr3__single_SemiBin_390`): keep with annotation (current draft) or remove from the 5,438. PI call. Never write anthrax | R1, R2, Discussion |
| A3 | R5 closing wording: "known parts in new arrangements" was withdrawn after Pfam found no domain in 91 % of members; current text reads "parts that resemble known folds structurally but are not recognisable sequence domains". Alejandra considers the structural result important; revisit the wording with her | R5 |
| A4 | Weight of the contamination material in the Discussion (competing-explanations paragraph and the Methods blank paragraph); candidate trim | Discussion, Methods |
| A5 | Target journal. As it stands Nature Communications (5,000 words, no trim). Nature Microbiology needs A1 and a 600-word trim (R3 to one paragraph, R5 to ~150 words, R6 prediction paragraph to one sentence). Decide at the end | word budget |
| A7 | Authorship and byline | front matter |
| A8 | Claim registry `reproducibility/claims.tsv`: no rows for the v2 sections, withdrawn values elsewhere. PI rebinds at end-game | Methods |
| A9 | XRF units and instrument; aridity climate source and years; fifth-expedition handheld XRF in or out | Methods (sampling), BLOCKED on her |
| A10 | Data availability: genome deposition target, archived-release DOI, companion KG citation | availability |

## B. Compute

| # | item | state | cost | unblocks |
|---|---|---|---|---|
| B1 | Known relatives of the 237: DONE and written 2026-09-24 (0 of 237 at species rank; best AAI median 64.8 %). Optional follow-up: the 49 SMAG genomes of the two classes once the SMAG archive from B2 is on disk | optional | minutes | R4 |
| B2 | SMAG in the R6 absence search: DONE, registered and written 2026-09-24 (40.47 % absent, gap +27.69; union with SMAG 26.83 / 0.36). SMAG JOINED the headline union 2026-09-27 (A11): 26.83 / 0.36 %, depth 5.4 / 35.8 / 25.37 %, written and registered. Occupancy on the new absent set registered (37 vs 38 unchanged). Remaining: Fig 5b SMAG point (see C5) | done | | Fig 5 |
| B5 | Gene-level check of the compartment contrasts on all assembled genes, unbinned half included | revision item | hours | R7 |
| B6 | Within-population pN/pS on dark and control member genes (bowtie2 + inStrain on 551 samples; 3,000–6,000 CPU-h) | revision item | days | R2 |

## C. Writing and display

| # | item | bears on |
|---|---|---|
| C1 | Commit the three Discussion fixes now in `main.tex`: Maeke wording (reference composition, not genome quality; "sister" removed), trace-gas comparator like-for-like (968 of 5,438 vs 58 %), Thermoplasmatota precedents (Zhou 2019, Sheridan 2022, Finstad 2017) | Discussion |
| C2 | Length as a fifth competing explanation: 76 % of dark representatives under 100 aa; answer with 67.1 % persistence at ≥250 aa, the frame-shift null, cultured dark median 59 aa; cite `duan2024gmsc` | Discussion |
| C3 | Intro paragraph 1: scale sentence naming Microflora Danica (10,683 samples, 5,518 species; `singleton2026danica`) and SMAG (3,304 metagenomes, 21,077 SGBs; `ma2023smag`). Intro paragraph 3: narrow the foil, global catalogues have the family layer and no context, regional catalogues the context and no controlled family layer | Intro |
| C4 | R4: never write TMEG, MBG-D, RC-III or Lutacidiplasmatales as former names of E2 or SW-10-69-26 (those are c__Thermoplasmata) | R4 |
| C5 | Figures: Fig 1 redraw in vector; Figs 2 and 3 unbuilt; Fig 6b shows the superseded geochemistry block (rebuild on the eleven-axis block, add the partition); Fig 6c and 6e still on the 548 cohort (rebuild from `results/v551/`, 45 sites, n 250); Fig 5b needs the SMAG point | display |
| C6 | Extended Data Table 1 (censuses) and Supplementary Table S3 (stress panels) are referenced with no file | display |
| C7 | Tawfiq companion and data-descriptor bibkeys (from Alejandra) | Intro, Methods |
| C8 | Re-check Hutchinson's 58 % (322/554) against the published ISME J version | Discussion |
| C9 | Trim to 3,500 only if A5 = Nature Microbiology | whole paper |

## D. Records to fix

- R1 order-gap sentence accounts for 8 of the 10 adjudicated candidates (six single-genome, two multi-genome clades); find the other two in the adjudication path and state them; the GUNC range 0.01–0.32 is not in the registry, register it.
- `CATALOGUE_NUMBERS.md` §2.3 basis line still says "44 site-paired"; the table says 45.
- `PAPER_STRUCTURE.md`: says seven sections with geochemistry as R7; the paper has eight, R1 and R6 headings differ. Align.
- `CATALOGUE_NUMBERS.md`: relabel §1.1–§3.1 to R1–R8 at the true-up; register the adjudication path; move 10.8x and "1.5x, p = 0.011" to the superseded map.
- `PIPELINE_GUIDE.md`: §8 provisional header, stale `refbias_2_magko.py:82` pointer, r226 warning, "GUNC never ran" false for the 10 candidates.
- `reproducibility/claims.tsv`, `REPORT.md`, `run_all.sh` on withdrawn numbers (PI, A8).
- `env_matrix.tsv` cannot be rebuilt; the availability statement ships the matrix.
- Supplementary Note 1 (energetics) survives only if the trace-gas ED item needs it; its "60–80 Tg yr⁻¹" conflates terms.

## E. Blocked

| item | on |
|---|---|
| Sampling-methods facts (XRF, aridity source) | Alejandra (A9) |
| Availability: reads = PRJEB125153 (554 runs, hold 2026-11-01); sample-level accessions from `code/ena_deposit/09_fold_accessions.py` once the deposit closes; genomes and code target | Alejandra / PI (A10) |
| Contaminant representatives in or out | PI (A2) |
| Plant-host genus axes (amplicon cache off-IBEX) | infrastructure |

## Not doing

- Nothing from the exploratory folders enters the paper without a registered rerun.
- Nothing from the `dark_function` annotation project enters this paper (2026-09-20).
- No further characterisation of the dark proteome beyond C2; no pangenome of the Thermoplasmatota on genomes below 90 % complete.
- The environment matrix is not rebuilt; the SMETANA cross-feeding pilot is cut.
- No geographic-endemism claim; database-relative absence under a named search rule only.
- No "radiation", no "genomes worldwide", no "hydrogenase-rich", no anthrax, no 10.8x, no 1.5x.
- Strain persistence across expeditions is the next study, not this paper.
