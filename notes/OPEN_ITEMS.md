# EQ-MAGs: open items
Purpose: the open-items register for the PI's rewrite (15 and 20 Oct 2026: five Results sections, controls in Methods and supplement, 19 authors, Nature Communications as stretch target).
Pointers: structure `PAPER_STRUCTURE.md`; numbers `CATALOGUE_NUMBERS.md`; changes `CORRECTIONS.md`; manuscript = Overleaf main (`paper/` submodule in empty-quarter-mags; snapshot in `manuscript/`).
Rule: nothing unregistered enters the paper. M = main.tex line, S = supplementary_notes.tex line.

## A. Decisions for Alejandra
| # | item | recommendation | bears on |
|---|---|---|---|
| A15 | Remove 17 orphan figures (+ options/), supplementary.tex, main_v1.tex from Overleaf | no, until submission | Overleaf |

## B. Decisions with the PI
| # | decision | options | recommendation |
|---|---|---|---|
| B1 | Where paper and record live | (a) text in Overleaf main, record in this repo, old repo archived; (b) edit only through him | (a); push to main, one line to him |
| B2 | SMAG in the headline union, applied in the text (28,446 families, 26.83 % / 0.36 %, 1.387 G; depth check 25.37 %) | his confirmation; fallback = his 30,479 / 1.258 G / 27.08 % | one line to him |
| B3 | Journal | Nature Communications vs decide after P3 | Nature Communications; no trim needed |
| B4 | Structural excess as a claim | (a) 6.6× / 8.1× standardised with his caveat; (b) crops-only framing | (b) in main, (a) in Note S2 |
| B5 | Authorship and contributions | 19 authors as set; contributions statement missing | accept; draft contributions from the record |

## C. Compute
| # | item | state |
|---|---|---|
| C1 | Gene-level compartment check on all assembled genes, unbinned half included | revision item |
| C2 | Within-population pN/pS on all dark and control genes (bowtie2 + inStrain, 551 samples); only the Note S4 read-level pilot exists | revision item |

## D. Writing and display
| # | item | waits on |
|---|---|---|
| D1 | Fig. 4b: SMAG point drawn (verify after the figure rebuild); SPIRE +21.30 not in text, optional | Alejandra |
| D3 | Contributions statement | B5 |
| D4 | Phrase check: his `verify.py` reads only `manuscript/main.tex`; the 20 Oct split moved required phrases into `supplementary_notes.tex` | PI fix |

## E. Blocked
| item | on |
|---|---|
| ENA PRJEB125153 hold until 2026-11-25 (runs, assemblies, bins, MAGs submitted); archive DOI; umbrella project (PI) | PI |
| Fe > 45 in F17Sr1, S44PRr1, S56PRr1 (data-descriptor XRF tables) unexplained | PI / data owner |

## Not doing
- Nothing from the exploratory folders enters the paper without a registered rerun.
- Nothing from the `dark_function` annotation project enters this paper (2026-09-20).
- No further characterisation of the dark proteome beyond the length clause in main.tex; no pangenome of the Thermoplasmatota on genomes below 90 % complete.
- The environment matrix is not rebuilt; the SMETANA cross-feeding pilot is cut.
- No geographic-endemism claim; database-relative absence under a named search rule only.
- No "radiation", no "genomes worldwide", no "hydrogenase-rich", no anthrax, no 10.8x, no 1.5x.
- Strain persistence across expeditions is the next study, not this paper.
