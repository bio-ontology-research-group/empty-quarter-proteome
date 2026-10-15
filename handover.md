## Latest P3 state — 12 October 2026

The three-round scientific review/edit loop is complete: Astra rounds 1–2, Gemini
3.1 Pro High through agy round 3 at Robert's request after Opus failed. See
lifeonair/reviews/round3/response_round3.md and STATUS.md. Seven recovered expedition-4
blank candidates, expanded contamination sensitivity and corrected benchmark scope
are included. Group-1l relative depletion remains the main conclusion. The P3
verification passes; the repository-wide P4 harness retains its separate 30 existing
gate failures (round3/repository_harness.log). Public archive/accession confirmation,
portable public upstream reproduction and author/funding metadata remain required
before submission. Unrelated results/relic_climate/ was left untouched.

# Empty Quarter MAGs — handover

## Next agent: read first (9 Oct 2026)

Robert split the Empty Quarter work into four papers. The plan, data coverage and ordered task list are in
`plans/manuscript-split-2026-10.md`. P3 (life on air, ISME Journal) is in analysis in this repository. P4
(dark proteome, the rest of `manuscript/main.tex`) is Alejandra's.

**P3 review loop (12 Oct 2026).** A 3-round critical-review loop is running (Flash literature expansion via
`agy`, Astra review via `codex exec`, then revision). Round 1 is complete: review and response in
`lifeonair/reviews/round1/`, corrections in `notes/CORRECTIONS.md` (2026-10-12), current analyses in the
`*_v2` result folders (`lifeonair/results/README.md`, top section). H1 is now *relative* depletion of group 1l
beside roots, H2 is narrowed to group-1l specificity, and the bulk aridity association is reported as not
robust to spatial blocking. Rounds 2 and 3 follow. The H1–H3 list below is historical.
**P3 status (10 Oct 2026).** All P3 material is in `lifeonair/`; read `lifeonair/results/README.md` first.
The life-on-air passages of the P4 manuscript are commented out (`% P3 LIFE ON AIR`) and copied to
`lifeonair/manuscript/removed_from_P4_2026-10-10.tex`.
**P3 manuscript (10 Oct 2026).** Full ISME J draft at `lifeonair/manuscript/main.tex` (build with
`latexmk -pdf main.tex` in that folder; tables from `lifeonair/code/make_tables.py`). Figure 1
`lifeonair/figures/fig_h123.pdf` (`fig_h123.R`, hydB tree now rooted on group 2), Figure 2
`fig_mechanism.pdf` (`fig_mechanism.R`). Open placeholders: author list, data-descriptor and catalogue
citations, ENA release date, archived-release DOI, acknowledgements. Claims not yet in `reproducibility/claims.tsv`.
H1 mechanism proxies (10 Oct 2026, `lifeonair/results/h1_mechanism/`): roots shift communities to faster
growers and greener sites have less 1l, consistent with Li et al. 2023, but a 1l-specific loss remains beyond
growth rate and phylum (soil organic carbon was never measured). The plan's original hypothesis
(stable baseline despite strong turnover, compartment-independent) is contradicted by the data, so it was
replaced by three hypotheses:
- **H1:** roots selectively deplete group 1l H2 scavengers at the same 5–10 cm depth (about 4.7-fold,
  beyond control genes, in all 5 trips).
- **H2:** H2 and CO capacity are decoupled across microhabitats (form-I CoxL is invariant).
- **H3:** long-term aridity sets the spatial 1l baseline, but between-trip change does not follow rainfall.
  Caveat: the aridity effect is confounded with the main community-composition axis.

Analyses done:
- read-based guild genome equivalents for all 573 samples (IBEX array 53524441)
- MAG contamination screen
- stability and turnover models
- operon context and trees
- a draft figure

Open:
- register claims in `reproducibility/`
- write the manuscript
- decide the H3 framing
- optional gas-uptake assays on archived soil

## Start here

The public, buildable manuscript is `manuscript/main.tex`. The `paper/` directory is an optional
private Overleaf editing submodule; a clean clone does not need it. Run the complete lightweight
reproducibility package from the repository root:

```bash
bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

Do not run assembly, read mapping, large database searches, folding, or catalogue clustering on a
workstation. Those workflows belong on IBEX (CPU/high-memory) or ORIX (GPU). The local package is
deliberately limited to compact-table recomputation, figure generation, checksums, and LaTeX.

## Analysis universes

Several defensible results use different frozen workflow products. They must not be silently
joined or treated as nested denominators.

| Object | Size | Used for |
|---|---:|---|
| paired-read sample universe | 573 | canonical manifest and current mapping rebuild |
| metadata sample universe | 551 | ecology and abundance matrices |
| current refined-bin catalogue | 14,426 | current assembly inventory |
| current ANI95 HQ/MQ representatives | 5,438 | GTDB-Tk 2.7.2/R232 taxonomy |
| historical functional/Bakta MAG snapshot | 15,469 | 51.9-million-protein functional analyses |
| historical recurrent dark orthogroups | 23,382 | dN/dS, ecology, abundance, and specificity |
| historical assembly protein catalogue | 753,305,533 proteins | separate structural-tier screen |

The 15,469 and 14,426 counts come from different curation snapshots and inclusion rules. Their
difference is not a gain/loss statistic. The manuscript labels the historical snapshot explicitly.

## Current taxonomic result

GTDB-Tk 2.7.2 against R232 classifies 5,438 ANI95 HQ/MQ representatives across 51 named phyla.
First assignment gaps are mutually exclusive: zero at domain/phylum/class, 11 at order, 142 at
family, 1,343 at genus, and 3,755 at species. Only 187 representatives are assigned through
species; these contain 184 distinct named species labels. Ten order-gap representatives were
staged into de novo phylogenies and none was automatically advanced as a novel-rank proposal. One
low-MSA order-gap candidate remains unadjudicated.

This supports extensive lower-rank novelty and targeted recovery/coassembly work. It does not yet
support a new phylum or class. LQ bins, including 304 separately tracked coassembly-salvage bins,
are prospective recovery material and are not mixed into the current HQ/MQ taxonomy claims.

## Current dark-proteome result

The historical 51.9-million-protein ladder yields 23,416,070 90%-identity representatives and
2,839,030 UniRef50 no-hit representatives (12.1%). The predicted-3Di AFDB screen is retained as a
historical triage stage, not a count of proven folds.

The current structural sensitivity analysis uses full ESMFold coordinates, PDB through 30 July
2026, AFDB v6/UniProt 2025_03, two TM-score normalisations, explicit retrieval-failure accounting,
and known-fold controls. Among length/confidence-qualified candidates, query-normalised TM gives
608/1,875 dual no-hits and alignment-normalised TM gives 97/1,875. Identically gated known-fold
controls give 101/3,013 and 68/3,013, respectively. The query-normalised known-fold error rate is
3.35% (Wilson 95% CI 2.77–4.06%). These are reference-relative candidates with non-zero empirical
error floors, not validated new folds.

The historical recurrent-orthogroup object contains 23,382 families selected under the deposited
darkness rule: representative no UniRef50 match or less than 30% identity, with at least three
members. Of 22,767 evaluable families, the descriptive median pairwise NG86 dN/dS is 0.391 and
rises with the minimum family-size floor. This is evidence of coding-like constraint conditional
on recurrent alignable ORFs, not a catalogue-wide selection coefficient.

The historical dN/dS-greater-than-one tail contains 148 screened families; 123 pass a recorded
minimum-dS/pair-count heuristic. Eight have confident structural analogues and five analogue
descriptions motivate recognition/binding/surface/transport hypotheses. The recomputed pocket
tally is 23/332 nominal sites. This exploratory result does not establish positive selection,
interface diversification, or an arms race.

## Specificity and function

Specificity is always relative to named catalogues and the stated DIAMOND criterion. Of 23,382
historical dark orthogroups, 5,011 (21.4%) have no detected homolog in the named MGnify90, GEM, and
five-desert comparison. This is not evidence of universal absence or geographic endemism.

The functional snapshot supports an atmospheric trace-gas guild, near-universal repair systems,
a broad BGC census, a sparse named-pathogen signal, decoded enzyme candidates, and strong
database/reference ceilings. Capacity is not flux and structural/ligand annotations remain
hypotheses unless experimentally validated.

The historical PICRUSt2/MAG comparison remains labelled historical. A fail-closed current-catalogue
replacement is dependency-gated on all 573 raw-read CoverM profiles, matrix closure, and sample-
universe validation. Historical pangenome and BGC-family analyses are retained as exploratory
upper-bound records and are not promoted as current manuscript results.

## Reproducibility architecture

- `reproducibility/current_manuscript_artifacts.tsv` pins compact frozen inputs.
- `reproducibility/historical_provenance_ledger.tsv` pins heavy upstream counts to exact sources.
- `reproducibility/dnds_legacy_provenance.tsv` enforces the historical darkness/dN/dS scripts and
  result hashes, including withdrawn pre-guard code from immutable Git objects.
- `reproducibility/large_files.sha256` pins seven deposited large inputs mirrored on IBEX.
- `reproducibility/verify.py` recomputes or checks registered claims and scans active prose for
  superseded language.
- `reproducibility/regenerate_figures.sh` sets a fixed epoch/time zone and regenerates figures.
- `reproducibility/build_manuscript.sh` builds the public snapshot and rejects divergence from an
  initialised editing submodule.

Use `reproducibility/claims.tsv` and `reproducibility/REPORT.md` as the audit entry points. Plans in
`plans/` record historical intentions and are not evidence that an analysis was delivered.

## Remaining scientific work

1. Complete and validate the 573-sample current CoverM → MAG-KO → PICRUSt2 comparison on IBEX.
2. Finish AFDB/PDB/MGnify triage; stage only unresolved sequences to ORIX and run actual CUDA
   Foldseek/ESMFold workloads there.
3. Reassemble or coassemble samples carrying promising MQ/LQ and order-gap signals, then rerun bin
   refinement, contamination checks, ANI dereplication, GTDB R232 classification, read support,
   and phylogenetic adjudication. New high-order taxa are a hypothesis, not a current result.
4. Confirm the B4 candidate tail with recombination-aware codon models before making any selection
   claim.
5. Keep defensome/mobilome analyses in the sister phages repository and join them to MAG context
   only through versioned, checksum-pinned interface tables.
