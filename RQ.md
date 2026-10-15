# Research questions — Empty Quarter MAGs and proteome

This register separates current, manuscript-supported results from historical exploratory products
and planned work. “Done” means that the bounded statement is represented in deposited evidence;
it does not mean every upstream heavy computation can be rerun locally.

| Status | Meaning |
|---|---|
| **Current** | Included in the current manuscript with an executable or checksum-pinned claim path. |
| **Historical** | Reproducible snapshot from an older analysis universe; explicitly labelled and not substituted for the current catalogue. |
| **Exploratory** | Useful hypothesis/upper bound, but not promoted as a manuscript conclusion. |
| **Planned** | Data or design exists, but the requested result has not been delivered. |
| **Wet lab** | Requires new experimental evidence. |

The current taxonomy universe is 14,426 refined bins and 5,438 ANI95 HQ/MQ representatives under
GTDB-Tk 2.7.2/R232. The 15,469-MAG functional/Bakta snapshot, 23,382 historical dark orthogroups,
and 753,305,533-protein assembly catalogue are separate workflow products. See `handover.md` and
`reproducibility/claims.tsv` before combining any result.

## A. Dark-proteome scale and interpretation

| RQ | Status | Bounded answer |
|---|---|---|
| A1 | Historical + current calibration | The historical ladder contains 2,839,030 UniRef50 no-hit representatives (12.1%). Predicted-3Di no-hits are triage candidates; the current full-coordinate PDB/AFDB analysis reports both TM normalisations and empirical known-fold error floors. |
| A2 | Current, descriptive | Recurrent alignable dark ORFs show a median pairwise NG86 dN/dS of 0.391 and recurrent abundance. This supports coding-like constraint conditional on the filters; it is not a catalogue-wide purifying-selection estimate. |
| A3 | Current | At length at least 100 aa and mean pLDDT at least 0.70, qTM yields 608/1,875 candidate dual no-hits versus 101/3,013 known-fold controls; alignment-normalised TM yields 97/1,875 versus 68/3,013. Survivors require experimental confirmation. |
| A4 | Historical, bounded | A calibrated structure/context pipeline prioritises ten specific-EC candidates and three pocketed structural targets, with a low stated correctness ceiling. |
| A5 | Historical, negative/bounded | Decoded candidates map mainly to known reaction space. No genuinely new enzymatic reaction is established. |

## B. Ecology and biogeography of the historical 23,382-orthogroup object

| RQ | Status | Bounded answer |
|---|---|---|
| B1 | Historical | Dark-family composition is associated with environment after the deposited spatial controls; the result is conditional on the historical abundance matrix. |
| B2 | Historical, bounded | 5,011/23,382 (21.4%) have no detected homolog in the named MGnify90, GEM, and five-desert comparison under the stated rule. This is database specificity, not universal absence or geographic endemism. |
| B3 | Historical | Deposited association tables identify water/salt/geochemical structure; trace-metal nulls are retained alongside positive associations. |
| B4 | Historical, negative | The family-level vanadium signal is null; no vanadium–alternative-nitrogenase coupling is claimed. |
| B5 | Historical | About 30.6% of families meet the deposited compartment-specificity rule, with rhizosphere and deep-soil signals prominent. |
| B6 | Historical, descriptive | Site-restricted and EQ-wide families have descriptive median NG86 values 0.332 and 0.420. Clustered/phylogenetic inference is pending. |
| B7 | Historical | Deposited variation-partitioning results favour environmental filtering over pure space, within their model and sample universe. |
| B8 | Historical | Most families recur across expeditions; trip effects are inseparable from sequencing/batch effects in the current design. |

## C. Reference bias and predictability

| RQ | Status | Bounded answer |
|---|---|---|
| C1 | Historical; current rebuild running | The old PICRUSt2/MAG comparison shows a 10,543-KO constant PICRUSt2 support and median 4,482 MAG KOs, but richness ratios depend strongly on threshold symmetry. A 573-sample raw-read CoverM rebuild must validate before replacement. |
| C2 | Exploratory, not promoted | A historical pangenome analysis exists, but its catalogue/provenance is not part of the current manuscript evidence chain. Treat the proposed open-pangenome mechanism as unconfirmed here. |
| C3 | Historical | KO coverage and taxonomy-novelty associations are descriptive properties of the historical functional snapshot. |
| C4 | Historical + current caveat | Structure can supply analogue hypotheses beyond sequence, but reference version, TM normalisation, retrieval failure, model confidence, and the known-fold error floor materially change the apparent rescue/no-hit rate. |

## D. Extreme-survival biology

| RQ | Status | Bounded answer |
|---|---|---|
| D1 | Historical | The functional snapshot contains widespread high-affinity H2/CO marker capacity and a 1,204-MAG dual-marker core. Capacity is not activity or flux. |
| D2 | Historical | HydDB/phylogenetic classification supports group-1l/1h hydrogenases and form-I CODH calls under the deposited pipeline. |
| D3 | Historical | CoverM-weighted marker abundance provides a relative trace-gas capacity index, not a biogeochemical flux measurement. |
| D4 | Current (2026-08-13) | On the 5,438 catalogue with a site-paired sample-level design: NER/BER near-universal (99.5–99.6%); light-dependent defences enrich at the surface (photolyase 24.3 vs 12.1% of the mapped community) and NHEJ enriches below (80.0 vs 85.7%), replicated within phyla; mismatch repair fails the phylum control; within-species ecotype test null (species sorting). |
| D5 | Current (2026-08-13) | First proper synthesis census: 79.5% of species carry ≥1 compatible-solute pathway, flat across compartments; salinity is answered by turnover toward Halobacteriota (S rho +0.45, p=7e-14, composition-robust), which lack ectoine/GG/MG. The old "OR=1.10" class-enrichment negative is a phantom number (source reports MWU p=0.11) and belongs to the §2.5 dark-family decode. |
| D6 | Historical, negative | Rhodopsin calls are uncommon and do not show the proposed compartment gradient. |

## E. Evolution and genome dynamics

| RQ | Status | Bounded answer |
|---|---|---|
| E1 | Current, descriptive | The filtered recurrent-ORF dN/dS distribution is below one in median, with a small historical greater-than-one screen tail. Neither quantity alone establishes a selection regime. |
| E2 | Exploratory | The 148-family historical screen has 123 heuristic survivors, eight confident structural analogues, and five broad recognition/binding/surface/transport hypotheses. The 23/332 pocket tally is descriptive; no positive-selection, interface, or arms-race conclusion follows. |
| E3 | Planned | Test dark-family enrichment on genomic islands, plasmids, and phage-associated regions using versioned interfaces to the defensome/mobilome repository. |
| E4 | Historical | Most historical dark families are taxonomically narrow under the deposited family-to-taxon table; this is conditioned on catalogue construction and detection. |

## F. Natural products and biotechnology

| RQ | Status | Bounded answer |
|---|---|---|
| F1 | Historical | The antiSMASH census reports terpene, NRPS, and ectoine capacity; it is a broad BGC census, not proof of product expression. |
| F2 | Exploratory, not promoted | The historical BiG-SLICE/MIBiG analysis is an upper bound because about 90% of called BGCs touch contig edges. A complete-BGC, current-catalogue rerun is required before a novelty conclusion. |
| F3 | Historical candidates + wet lab | Decoded enzymes and pocketed proteins are expression candidates; activity and extremotolerance need experiments. |

## G. Community metabolism

| RQ | Status | Bounded answer |
|---|---|---|
| G1 | Historical pilot | Sixteen-sample gapseq/SMETANA results are competition-leaning and do not preferentially anchor exchange to the trace-gas guild. |
| G2 | Planned | Biome-scale cross-feeding remains a high-cost cluster analysis. |
| G3 | Partial | Only decoded dark reactions can enter current GEMs; broader participation remains unresolved. |

## H. Comparative context and rigor

| RQ | Status | Bounded answer |
|---|---|---|
| H1 | Partly addressed | Recurrence, abundance, named-pathogen screening, and cross-catalogue comparison argue against a trivial contamination tail; a formal negative-control/decontamination analysis remains desirable. |
| H2 | Historical, bounded | The five named desert comparisons retain substantial no-hit fractions. This establishes comparator-relative difference, not geographic uniqueness. |
| H3 | Current | R232 shows no first assignment gaps at phylum or class, but 11 order-, 142 family-, 1,343 genus-, and 3,755 species-level gaps. Ten tree-adjudicated order candidates yield no automatic new-rank proposal; one low-MSA candidate remains. |

## I. Synthesis and validation

| RQ | Status | Bounded answer |
|---|---|---|
| I1 | Planned | A held-out sequence/structure/context/ecology predictor benchmark has not been delivered. |
| I2 | Planned | A current-catalogue minimum survival-genome analysis has not been delivered. |
| I3 | Framing | Astrobiology relevance is discussion context, not a result established by these data alone. |
| I4 | Wet lab | Candidate expression, ligand binding, enzyme activity, and trace-gas measurements require new experiments. |

No aggregate “answered” tally is reported because current, historical, exploratory, planned, and
wet-lab statuses are not equivalent. Per-question evidence and scripts are indexed in `RQ/INDEX.md`.
