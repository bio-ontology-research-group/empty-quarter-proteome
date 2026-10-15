# Paper structure — LOCKED 2026-08-09

Single source of truth for the empty-quarter-mags **Results** structure, decided after a full evidence
audit (branch `main`, PI's canonical rebuild). Supersedes the structure sketch formerly at `RQ/PAPER_NARRATIVE.md`, which was merged into this file
and `PLAN.md` and then deleted (see `notes/README.md`); its numbers were also dead. **Do not re-litigate without recording the change here.**

## Decision
- **3-part narrative, dark-matter first** — matches the title ("Known phyla, unknown biology"), the abstract,
  and the concept figure: Part 1 the dark matter → Part 2 function → Part 3 specificity.
- **One level of short descriptive subsection headings** (`\subsection*`); NO printed "Part N" labels, no
  nested sub-subtitles.
- The PI's **Part A/B/C scaffold is RETAINED** in `main.tex` as the frozen claim/number map — never deleted;
  fenced out of the PDF (`\iffalse…\fi`) once each subsection is ported. ([[feedback_keep_scaffold_guides]])
- **Current numbers only.** DEAD: 9,253 reps / R226 / 56 phyla / 10,007 folds / dN/dS-vs-0.29-control.

## Structure — each line is one `\subsection*`; [code] = main.tex scaffold source; status flag
### Part 1 — The catalogue and its dark matter
- **1.1 A taxonomically conventional catalogue of the world's largest sand sea** — [A0/H3] — ✅ current, done
  (14,426 bins → 5,438 reps · R232 · 51 phyla · first gaps 0/0/0 domain-phylum-class, 11 order, 142 family,
  1,343 genus, 3,755 species · 10 order-gaps adjudicated → 0 proposals).
- **1.2 A large fraction of the encoded proteome has no reference homolog** — [B1] — ✅ current
  (753M→449M reps → 69.9M/15.6% UniRef50-dark; confident ≥100aa 5.4%; genome-resolved 4.29M/9.89%;
  **length not completeness** discriminates darkness — 76% <100aa).
- **1.3 Structural novelty is bounded and calibrated** — [B1a] — ✅ DONE on the current dark set (1,466,943 queries → 68.0% triage no-hit → 17,221 fold families → **1,031 of 2,838 confident models (36.3%) dual no-hit under query-normalised TM vs a 3.35% control floor**; 95 (3.35%) under alignment-normalised TM vs its own 2.26% floor). Architectural, not local: 96.7% of confident models still contain a locally familiar region.
  (608/1,875 = 32.4% qTM dual no-hit; 97/1,875 = 5.2% aln-TM; known-fold control floor 3.35%. NOTE: on the
  historical tier1 ESMFold cohort; the current 4.29M-dark FoldSeek search is running, job 50253325).
- **1.4 The dark families are conserved and abundant, not spurious** — [B2 dN/dS + B3 abundance] — ✅ current, done
  (2026-08-14, rebuilt universe: 106,042 families / 536,437 member CDS from §1.2's 4.29M dark reps; median 12
  samples / 9 sites; dN/dS 0.425 (96.7%<1, n=97,200) vs matched annotated control 0.388 (99.6%<1); median
  98/550 samples, 4 single-sample, median peak 99 TPM; 79.4% all-5-trips (94.8% ≥4), 28.9% re-assembled all 5;
  97.9% single-phylum of the 74% attributable; 108,826→130,843 clusters at 30%. Written + Methods M.7 + Supp
  Table S5. June 23,382 / 0.391 / 52-sample numbers SUPERSEDED; old 0.29 comparator withdrawn.)

### Part 2 — What these organisms do, and why references miss it   (bias → known → dark)
- **2.1 Reference and marker-gene methods cannot predict this gene content** — [C1 PICRUSt] — ✅ DONE on
  5,438 (2026-08-16). NOT blocked: the gate was inherited/stale (CORRECTIONS 2026-08-16). On 448 matched
  samples: NSTI **0.159**, presence-recall **86.0%**, ρ **0.839** on shared KOs; ceiling **55.6% of
  protein-coding genes carry no KO** (58.4% CheckM2 denominator; 53.9% HQ-only — not a quality artefact).
  Second ceiling: **2,517/11,374 (22%)** encoded KOs outside the predictor's family set, 14.0% per sample.
  Quote the no-KO ceiling, NOT the over-prediction ratio (still withdrawn). Methods M.8.
- **2.2 The community runs on trace gases** — [A1 / RQ D1–D3] — ✅ DONE on 5,438 (2026-08-11) — **flagship known function** (968 reps 18% high-affinity [NiFe]-hydrogenase, 1,170 22% form-I CODH, 1,860 34% either, 278 5% both; MAG-direct abundance 18.4% median whole-community / 43.9% of mapped, 548 samples. (Recomputed on rep_ko_v2.) Energetic feasibility is the payoff; global sink <1% = one sentence. Capacity, not flux. Pending: 5,438 guild table + tree figure.)
  (22% MAGs high-affinity [NiFe]-hydrogenase, 25% form-I CODH, ~8% both; abundance-weighted ~13.5% of the
  community; climate-relevant H₂/CO sink. Capacity, not flux — no metatranscriptomics).
- **2.3 A genome pre-loaded for radiation, desiccation and salt** — [A2 + A3 / RQ D4–D5] — ✅ recomputed + drafted
  (2026-08-13, 5,438/rep_ko_v2, sample-level site-paired design. Floor: BER 99.6 / NER 99.5 / NHEJ 74.5 /
  photolyase 17.9% of species. Compartment: photolyase 24.3 vs 12.1 / carotenoid 70.3 vs 61.0 / NHEJ 80.0 vs 85.7
  % of mapped community, 44 paired sites, within-phylum replicated; MMR failed the phylum control → not claimed;
  within-species ecotype test null → species sorting. Osmolytes flat by compartment; salt = TURNOVER:
  Halobacteriota vs S rho +0.45 p=7e-14, CLR/detection-robust; haloarchaea lack ectoine/GG/MG, retain trehalose.
  CHANGE (2026-08-13): old honest-negative "osmolytes not enriched as a class, OR 1.10" REMOVED — phantom number;
  its source (results/darkmatter_ecology/OSMOLYTE.md) reports MWU p=0.11 and no OR, and the test belongs to the
  §2.5 dark-family decode, not A3. See CORRECTIONS 2026-08-13.)
- **2.4 A biosynthetic reservoir and an intrinsic resistome** — [A4 census + A5] — ✅ current, done
  (2026-08-15, 5,438 basis: 97.5% of species ≥1 BGC, terpene 62.1% / NRPS 25.1% / ectoine 16.5%
  [terpene-precursor 83.8%, separate]; AMR any 37.4%, β-lactam 24.6%; metal layer recovered, As 17.3% /
  Hg 7.5%, geogenic [Ca + aridity, geography-partialled]; biosafety control: 2 named pathogen species
  [0.04%], VFDB 5,435→202 at id≥40→90 = competition, not pathogenesis; ectoine antiSMASH×KO Jaccard
  0.92). **BGC novelty never computed (deferred); contig-edge 85.3% measured — census only.**
- **2.5 Most of the function is dark, and environmentally structured** — [B3 gene–environment] — ✅ current, done
  (2026-08-16, rerun on the 106,042 × 551 matrix: PERMANOVA geochem R² 0.123 / beyond-space 0.080 / climate|space
  0.030 / compartment 0.044, all p=0.001, vanadium null; family level S 54,227 / Ca 47,174 / Cl 46,972 / salinity
  37,280 / aridity 23,031 at q<0.05 of 94,568–100,293 tested, 87–96% negative = a FILTER; space-partial retention
  71–86%, salinity proxy 41%; V/Mo/Ni/Co/Mn/Zn = 0. Per-protein decoding → **dark-matter-annotator sibling
  paper**: one bridge sentence, nothing more. Display: rebuilt fig:environment + Supp Table S6).

### Part 3 — Specificity to the Rub' al-Khali
- **3.1 The dark proteome and the reference gap** — [B2/C3/H2] — ⏳ recompute in flight (2026-08-17)
  (**REFRAMED 2026-08-17, approved — recorded per the do-not-re-litigate rule.** The old heading,
  "biome-specific, *not generic*", asserted an exceptionalism we never measure, and biome restriction is
  the published norm — Aplakidou 2025: 63.4% of 608,258 novel families are single-ecosystem. New claim:
  absent from the named references · not an artefact of being short and rare [length/size-matched
  annotated control, new] nor of database size [e-value-corrected rarefaction] · because the reference
  set holds no desert like this one. **Still LOCKED: database-relative specificity under an explicit
  search rule, NOT geographic endemism**; provenance-level wording only — "no Arabian or Saharan site in
  any named desert comparator", never "absent from all catalogues". June 21.4/47.4/36.6/31.9 DEAD.
  **DISPLAY DECISION 2026-08-17 — supersedes "rebuild all three as-is"; Alejandra delegated the call.**
  Part 3 ships **one table + one figure + one supplementary table**, with no duplicated data:
  · `tab:specificity` — the nine catalogues, sizes, dark % and control %. A table presents nine rows of
    exact numbers better than bars; an earlier draft duplicated this as a figure panel.
  · `fig:specificity` → **`fig_specificity_controls.pdf`, 2 panels, both inferential**: (a) observed vs
    size-matched expected absence, where "below the diagonal" is the content-beyond-depth result and
    MGnify90 on the diagonal is the internal check; (b) median site occupancy of reference-absent vs
    reference-matched families **within length strata**. Only these two need a picture.
  · Supp Table **S7** — the full gap table with expectations.
  · **`fig:endemism` is CUT.** Of its three panels, one (occupancy) became `fig:specificity`b, one showed
    the occupancy × dN/dS analysis dropped as not reportable, and one duplicated `tab:specificity`. Its
    dangling `\ref` in the B3 scaffold was repointed; the scaffold claim is unchanged. Rationale: an
    unstratified occupancy panel would display a length artefact (median rep length falls 145 aa at 1–3
    sites to 37 aa at 56–61 sites). Generator `code/make_specificity_figures.py`; do NOT rerun
    `make_global_specificity_figure.py` or `make_darkmatter_figures.py`, both on the dead basis.
  See CORRECTIONS 2026-08-17.)

## CUT — not in this paper (be critical; RQ/PI deferrals)
- Per-protein EC decode (B5 / RQ A4–A5) → **dark-matter-annotator sibling paper**.
- HGT / mobilome / defensome / plasmid / CRISPR (RQ E3) → **empty-quarter-phages repo**.
- Open-pangenome mechanism (RQ C2) → dropped by the PI.
- Community-metabolism SMETANA (RQ G1 / main C2) → 16-sample pilot; ≤1 line or cut.
- BGC novelty % (RQ F2) → upper bound only; census yes, novelty no.
- Fast-evolving dN/dS>1 tail (B4 / RQ E2) → ≤1 hypothesis-generating sentence in 1.4.
- Rhodopsin (RQ D6) → one honest-negative sentence.
- Isolates (I4), predictor benchmark (I1), astrobiology (I3 → Discussion line only), paleoclimate,
  within-population pN/pS (inStrain) → out / deferred to revision.

## Status legend
✅ current & done · 🔧 recompute on the 5,438 catalogue (write structure now, numbers provisional) ·
⏳ blocked / job running · 🟡 basis-to-verify (23,382-orthogroup June tables) · ➕ new analysis to add.

## Number sources to trust (audit 2026-08-09)
`main.tex` Sections 1–4 + `CLAUDE.md` + `notes/CORRECTIONS.md` + `notes/CATALOGUE_NUMBERS.md` (current
registry; paper-section keys since 2026-08-14) + `reproducibility/claims.tsv` (PI mirror; may lag corrections).
Do NOT trust taxonomy/dN/dS/fold numbers in `AUDIT` or `RQ/SPINE_DOSSIER.md` — they carry the dead
9,253/R226/10,007 line. (`RQ/PAPER_NARRATIVE.md` no longer exists; it was merged here and deleted.)
Build: `pdflatex + bibtex + pdflatex×2` (natbib/bibtex, NOT biber).

## Display set (session 09, 2026-08-28/30)

**Figure rules (approved):** panel letters only, titles and ratios in the caption; descriptive
noun-phrase caption titles (comparator survey: 0 of 10 use a claim); no em-dashes; three text
tiers (8 / 7.5 / 6.8 pt) on a print-true canvas; palette roles in `code/figures/palette.py`
(teal focal, coral control, violet known reference, grey de-emphasis); generators assert every
plotted number against the registry before drawing.

**Built:** Fig 1 study design (AI-drafted raster, vector redraw owed before submission) ·
§1.3 structural search, calibration and candidate architectures (four panels; superposition slate
selected by a stated rule).

**Triage proposed 2026-08-30, PENDING RATIFICATION with Alejandra.** Test applied: does the claim
need a picture, or is it numbers a sentence already carries?
- Main figures: study design · §1.3 · §2.5 environment ordination · §3.1 specificity (+ its table).
- Conditional fifth: §2.2 trace gases, only if the HydDB/CODH placement tree is built with
  reference sequences; the counts alone are prose.
- Demote to supplementary: §2.1 reference bias, §2.3 stress modes, §2.4 BGC and resistome
  (currently main Figs 2 and 3).
- Prose only: §1.1 and §1.2 (covered by tab:scale and tab:ladder), §1.4 (four numbers; a
  distribution figure was built and rejected as a decorated table).
