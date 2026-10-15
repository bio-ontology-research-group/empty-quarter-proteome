# Empty Quarter manuscripts: split plan (decided 9 October 2026)

Follow-up for the next agent. Robert decided the venue strategy below on 9 October 2026 after a
read-only review of the ecology and data papers (PA record:
`~/Public/software/pa/journal/2026-10-09-empty-quarter-venues.md`). Read `handover.md` and
`CLAUDE.md` first; their claim rules still apply (capacity is not activity, no em-dashes in prose,
heavy compute only on IBEX/ORIX, nothing is sent, pushed or submitted without Robert).

## Decision

| Paper | Repository / source | Target | Status |
|---|---|---|---|
| P1 Data descriptor (knowledge base, amplicon + geochemistry) | `~/Public/software/empty-quarter-paper-repos/data-paper`, Overleaf clone `~/Documents/papers/empty-quarter-data-current` (`sn-article.tex`) | Scientific Data (Data Descriptor) | Blocked on public ENA runs |
| P2 Ecology (amplicon biogeography) | `~/Public/software/empty-quarter-paper-repos/ecology`, Overleaf clone `~/Documents/papers/empty-quarter-ecology-current` | Environmental Microbiome, nearly as is (fallback FEMS Microbiology Ecology) | Figures integrated and fixes done on main (9 Oct 2026, repo 116c073/c335bf7, Overleaf f42df80/c5a64da); waiting on public ENA runs |
| P3 "Life on air" (atmospheric H2/CO oxidation) | this repository, folder `lifeonair/` | ISME Journal, short, one hypothesis | New; analyses below |
| P4 Dark proteome (rest of current `manuscript/main.tex`) | this repository | Decide after P3 | Life-on-air text commented out 10 Oct 2026; Alejandra rewrites |

Order: P1 first (citable dataset, ENA public), P2 within weeks, P3 after the analyses, P4 later.
P2 and P3 cite P1. P2 does not cite or mention P3 (Robert, 9 Oct 2026).

## Metagenome coverage (checked 9 Oct 2026)

Source: `~/Public/software/empty-quarter-workbench/data/metadata/sample_metadata.tsv` (551 samples
with metadata; the paired-read universe is 573).

| Trip | Deep | PlantRoot | Surface | Total |
|---|---:|---:|---:|---:|
| 1 | 69 | 90 | 51 | 210 |
| 2 | 11 | 21 | 5 | 37 |
| 3 | 40 | 61 | 27 | 128 |
| 4 | 19 | 47 | 19 | 85 |
| 5 | 27 | 49 | 15 | 91 |
| All | 166 | 268 | 117 | 551 |

61 sites; sites sampled in 1/2/3/4/5 trips: 4/2/18/30/7. Time and compartment can therefore be
tested in the metagenomes themselves; the earlier note that only 150 Trip-5 libraries exist refers
to the subset used by the ecology paper (CoverM, 991 genomes), not to this catalogue.

## P3: Life on air (ISME Journal)

### Hypothesis
SUPERSEDED 10 Oct 2026. The original hypothesis (a stable, compartment-independent baseline despite strong
turnover, with functional redundancy) is contradicted:
- trip explains about 2% of composition
- compartment strongly affects hydB
- carriers do not compensate for each other

It is replaced by H1 (roots selectively deplete group 1l), H2 (H2 and CO capacity decoupled) and H3
(aridity sets the spatial 1l baseline, not between-year rainfall). See `lifeonair/results/README.md`.
Analyses 1-5 below are done; 6 (gas uptake) is open; 7 (rainfall) is covered by H3.

### Evidence already in this repository
- A1: HydDB placement, groups 1l (2,609 MAGs) and 1h (893); form-I CODH (3,828 MAGs); guild core with
  both H2 and CO = 1,204 MAGs (about 8%). Tree figure pending.
- D3 (`RQ/D3`): 16S-weighted capacity index across 1,227 amplicon samples, median 13.5% (IQR 9.6-18.0%).
  Uniform by compartment: medians Deep 0.142, PlantRoot 0.128, Surface 0.137.
- Capacity by trip (PA computation 9 Oct 2026, join of `RQ/D3/table_capacity_per_sample.tsv` with
  `~/Public/software/empty-quarter/data/processed/functional/picrust2/merged/sample_metadata.tsv`):
  trip 1 n=315 median 0.133; trip 2 n=24 0.171; trip 3 n=478 0.128; trip 4 n=177 0.146; trip 5 n=233 0.138.
  Not yet a formal test; reproduce in code before use.
- B8 (`RQ/B8`): dark-proteome composition turns over strongly between trips (trip R2 0.168, pseudo-F 27.3;
  site R2 0.130, pseudo-F 1.25). Needs repeating on whole-community and guild composition.
- Honest negatives to keep: no within-biome aridity gradient for group 1l (rho = -0.172, p = 0.196, n = 58
  sites) and no monotonic cross-biome gradient; trace-gas enzymes are not supplied by the dark proteome.

### Analyses required (in order)
1. **MAG-level decontamination screen** before any abundance claim. Flag MAGs whose GTDB genus matches
   genera enriched in the amplicon extraction blanks and PCR controls (control registry in the data
   repository), plus a frequency test (abundance against library DNA input) on CoverM abundances.
   Report guild results with and without flagged MAGs. Full read-level decontamination is not required.
2. **Read-based guild abundance per metagenome** (IBEX): reads against the confirmed 1h/1l hydrogenase and
   form-I CoxL protein set, normalised to single-copy marker genes, giving the fraction of genomes carrying
   each capacity per sample. Replaces 16S-weighted capacity as the primary measure (16S index stays as support).
3. **Stability test**: guild abundance by trip and compartment, mixed model with site as random effect,
   restricted to the 37 sites sampled in at least four trips for the paired temporal contrast.
4. **Turnover versus stability**: compositional turnover of the whole community and of guild carriers across
   trips (Bray-Curtis/Aitchison, PERMANOVA trip vs site as in B8); quantify redundancy as change in which
   genera carry the capacity while the total stays constant.
5. **Genomic context**: operon neighbourhoods of confirmed coxL and hydB (accessory/maturation genes) as in
   `plans/atmospheric-metabolism.md` Stage D3; finish the HydDB/form-I tree figure.
6. **Optional, highest value (needs Robert):** measured H2 and CO uptake by gas chromatography on archived
   soils from a subset of sites, compartments and trips. Possible collaboration with the Greening lab.
   Check whether archived soil exists and how it was stored before proposing it.
7. Exploratory only: guild fraction against antecedent rainfall (same caveats as the ecology paper).

### Format
ISME Journal: about 5,000 words and 8 display items (check whether Methods count). Draft outline:
1 design and guild phylogeny; 2 guild prevalence by phylum; 3 capacity by compartment and trip;
4 turnover versus functional stability; 5 genomic context; 6 activity (if measured). Methods and
negatives go to the supplement. Data availability through P1 (ENA, Zenodo, knowledge base).

## P2: Ecology paper to Environmental Microbiome (nearly as is)
- Done 9 Oct 2026 (on main, not on branch `rund-figures-2026-10-09`): Rund's five figures generated in the
  pinned environment and synced to Overleaf; abstract 1,227; broken refs and panel letters (4e, 5b, 5c)
  fixed; 5c caption; function stated plainly as predicted; data availability statement added
  (ENA PRJEB106069 within PRJEB104209, data repository, Zenodo KG, code repository).
- Trace-gas markers: Robert decided (9 Oct 2026) to de-emphasize, not remove. Done: the MAG coxL 38.5% and
  [NiFe]-hydrogenase 22.6% stay with the marker table, described as broad KO-level sets that do not
  identify high-affinity trace-gas oxidation; the Actinomycetota attribution and the "hydrogenase carriers
  most abundant in shallow-subsurface soil" statement were removed (it would contradict P3's
  compartment-independence hypothesis); the Discussion keeps one sentence on genetic capacity.
  Do NOT reference P3 or any in-preparation study in P2. Note for P3: P2's definitions differ (hydrogenase
  set lumps groups 1 and 3d; coxL KO does not separate CODH forms; 975 MAGs from 150 Trip 5 libraries), so
  P3 numbers (HydDB 1h/1l, form-I CODH, 15,469-MAG snapshot) will not match P2's and need not.
- Still open: shotgun/MAG accession for P2's data availability statement (only amplicon accession given).
- Note the contamination handling across campaigns (samples with high contaminant fractions).
- Owner of the text: Rund Tawfiq and Susana; Robert decides authorship order.

## P1: Data descriptor to Scientific Data
- Blocker: ENA PRJEB104209 and PRJEB106069 show 0 public runs (ENA filereport API, 9 Oct 2026). Release
  them or provide reviewer access.
- Decide shotgun scope: deposit reads/MAGs, or describe the shotgun data in P3 only and remove it here.
- Add Author Contributions and an ethics/sampling-permit statement (MEWA collection and export, Nagoya);
  fix the contradictory Funding text; state MIxS checklist compliance.
- Expand Data Records (currently 361 words) with a per-file table: content, format, size, accession.
- Align the Zenodo record title (10.5281/zenodo.23134168) with the paper; remove "staging snapshot, no DOI"
  from the data repository README; move OWL/ShEx detail from Methods to the supplement or usage notes.

## P4: Dark proteome (remaining manuscript)
- Done 10 Oct 2026: life-on-air passages commented out of the Overleaf source (`% P3 LIFE ON AIR`; Results
  census sentence and hydrogenase root-loss clause, Discussion comparison, Methods typing). Cite P3 once it exists.
- Venue decided after P3; candidates depend on the structural results (keep the current honest framing:
  reference-relative candidates, known-fold controls, both TM normalisations).

## Open questions for Robert
1. Authorship and lead for P3 (and P2 order: Rund first?).
2. Gas-uptake measurements: archived soil available? collaborator?
3. Shotgun data in P1 or only in P3/P4?
4. ~~Whether P2 should mention the trace-gas result at all before P3 is out.~~ Decided 9 Oct 2026: keep, de-emphasized, no reference to P3.
