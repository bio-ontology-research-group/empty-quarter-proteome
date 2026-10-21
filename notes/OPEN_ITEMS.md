# EQ-MAGs: open items
Purpose: the open-items register for the PI's rewrite (15 and 20 Oct 2026: five Results sections, controls in Methods and supplement, 19 authors, Nature Communications as stretch target).
Pointers: structure `PAPER_STRUCTURE.md` (rewritten after A0); numbers `CATALOGUE_NUMBERS.md`; changes `CORRECTIONS.md`; manuscript = Overleaf main (`paper/` submodule in empty-quarter-mags; snapshot in `manuscript/`).
Rule: nothing unregistered enters the paper. M = main.tex line, S = supplementary_notes.tex line.

## A. Decisions for Alejandra
| # | item | recommendation | bears on |
|---|---|---|---|
| A0 | Section existence: his five sections vs the eight of 2026-09-20 (censuses to Note S3 and P3; R3, R8 into notes) | accept his five | structure, PAPER_STRUCTURE |
| A1 | Restore the matched-null sentence: 0.308 vs 0.048 (218 pairs, 6.5×), 0.310 vs 0.069 (237 pairs, 4.5×) | yes, no causal wording | M73 |
| A2 | Known relatives: 0 of 237 at species rank; best AAI median 64.8 % vs same-genus 75.4 %; nearest desert 73.1 %; skani step | yes | M75, S650 |
| A3 | Maeke clause: rarity and reference composition, not genome quality (3.9 % vs 22.5 %) | yes | M107 |
| A4 | Precedent sentence: Zhou 2019, Sheridan 2022, Finstad 2017 | yes | M85 or M111 |
| A5 | Microflora Danica beside SMAG (`singleton2026danica`) | yes | Intro M34 |
| A6 | Length clause: 76 % of dark representatives under 100 aa; 67.1 % of families ≥250 aa in all five expeditions | yes | M113, S898 |
| A7 | Structural floor beside raw counts: 1,031 of 2,838 (36.3 %) vs 101 of 3,013 (3.35 %) | yes; ratio waits on B4 | M213, S906 |
| A8 | Adopt his reruns: frame null 1.086 (n 9,029); root sets 20,136 / 41,046 | yes; mine to superseded | registry |
| A9 | 22 culture pairs (E. coli excluded) as the basis | yes | registry, S588–591 |
| A10 | "61 locations" → "60 sites" (settled with the companion papers 2026-08-17) | yes | M40, M46, M121, M176, M202, S408–417 |
| A11 | `tawfiq2026companion`, `rubalkhali_datapaper` instead of "in preparation" | yes if the bib entries are hers | S415 |
| A12 | Workbench STATUS lines in specificity_clean, uv_survival_clean, osmolyte_clean | yes (rob/ path, her word) | workbench |
| A13 | `setfacl` traverse on `/ibex/user/velazqal` for the PI (SMAG db) | yes; else 54 GB copy | B2 verification |
| A14 | One-line list to the PI: S380 path missing; natcomms_20261020 README/sbatch; "102 flagged", "6.31 → 6.34" without result file | yes | package |
| A15 | Remove 18 orphan figures, supplementary.tex, main_v1.tex from Overleaf | no, until submission | Overleaf |
| A16 | Title: his 20 Oct atlas, his 15 Oct archaea, or hers | his 15 Oct archaea title | title |
| A17 | XRF units and instrument (blocked since September) | answer in one line | Methods S440 |

## B. Decisions with the PI
| # | decision | options | recommendation |
|---|---|---|---|
| B1 | Where paper and record live | (a) text in Overleaf main, record in this repo, old repo archived; (b) edit only through him | (a); push to a branch plus one line to him |
| B2 | SMAG in the headline union | (a) 26.83 % / 0.36 %, 1.387 G; (b) his 30,479 / 1.258 G / 27.08 plus a with-SMAG sensitivity sentence; (c) drop | (a); his 18.54 matched test stays as the fairness test |
| B3 | Journal | Nature Communications vs decide after P3 | Nature Communications; no trim needed |
| B4 | Structural excess as a claim | (a) 6.6× / 8.1× standardised with his caveat; (b) crops-only framing | (b) in main, (a) in Note S2 |
| B5 | Authorship and contributions | 19 authors as set; contributions statement missing | accept; draft contributions from the record |

## C. Compute
| # | item | state |
|---|---|---|
| C1 | Gene-level compartment check on all assembled genes, unbinned half included | revision item |
| C2 | Within-population pN/pS on dark and control genes (bowtie2 + inStrain, 551 samples) | revision item |
| C3 | natcomms_20261020 folder needs README and sbatch | PI |
| C4 | Arsenic and climate gradients rerun on cohort 550 (resistome_clean, osmolyte_clean `results/v551`) | in progress 2026-10-21 |

## D. Writing and display
| # | item | waits on |
|---|---|---|
| D1 | SMAG union Methods sentence after S765; SPIRE +21.30 / SMAG +27.69 at M91; Fig. 4b SMAG point | B2 = yes |
| D2 | PAPER_STRUCTURE rewrite to his section table | A0 |
| D3 | Contributions statement | B5 |
| D4 | Phrase check: his `verify.py` reads only `manuscript/main.tex`; the 20 Oct split moved required phrases into `supplementary_notes.tex` | PI fix |

## E. Blocked
| item | on |
|---|---|
| XRF units and instrument | Alejandra (A17) |
| ENA PRJEB125153 public release 1 Nov 2026; MAG deposition and DOI; round7 submission prerequisites (raw paths, permissions, duplicate-run provenance, author contributions) | PI |

## Closed since 24 Sep
- A2 contaminants (PI: keep in the 5,438, exclude E. coli from culture pairs) · A3, A4 sections rewritten · A5 → B3 · A7 → B5 · A8 his harness replaces claims.tsv coordination · B1, B2 done
- C1–C3 done in her text, now transplants (A) · C5, C6, C8, C9 moot (new figure set, no ED table, trace gas in P3, no trim needed) · D records: registry updated 2026-10-21

## Not doing
- Nothing from the exploratory folders enters the paper without a registered rerun.
- Nothing from the `dark_function` annotation project enters this paper (2026-09-20).
- No further characterisation of the dark proteome beyond the length clause (A6); no pangenome of the Thermoplasmatota on genomes below 90 % complete.
- The environment matrix is not rebuilt; the SMETANA cross-feeding pilot is cut.
- No geographic-endemism claim; database-relative absence under a named search rule only.
- No "radiation", no "genomes worldwide", no "hydrogenase-rich", no anthrax, no 10.8x, no 1.5x.
- Strain persistence across expeditions is the next study, not this paper.
