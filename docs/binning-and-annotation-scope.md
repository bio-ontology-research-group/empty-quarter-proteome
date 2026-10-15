# Binning and annotation scope — current and historical universes

This document prevents counts from different pipeline generations from being combined. The active manuscript
uses a current refined/dereplicated catalogue for taxonomy and recurrent-family analyses, while several functional
censuses remain explicitly labelled historical snapshots. It does not claim that all functions were rederived on
one common catalogue.

IBEX roots used by the historical workflows are:

```text
M=/ibex/scratch/projects/c2014/rob/metagenomics
P=/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites
```

## Current canonical sample and MAG catalogue

The checksum-pinned source of truth is `reproducibility/canonical_catalogue_summary.json`, mirrored from the
workbench canonical manifest. Its distinct universes are:

| Universe | Count | Role |
|---|---:|---|
| paired-read samples | 573 | raw paired-read inventory |
| metadata samples | 551 | metadata-covered universe |
| binnable samples | 550 | inputs eligible for binning |
| coassembly samples / groups | 527 / 126 | staged coassembly salvage |
| DAS_Tool-refined bins | 14,426 | current bin universe |
| refined HQ / MQ / LQ bins | 4,225 / 8,209 / 1,992 | CheckM2/MIMAG quality strata |
| ANI95 HQ/MQ representatives | 5,438 | current nonredundant taxonomy/family universe |

The current HQ+MQ arithmetic is **4,225 + 8,209 = 12,434**. The frozen historical dRep membership table contains
12,433 inputs and yields 5,438 representatives. Thus one current HQ/MQ bin is outside that historical membership
snapshot. This one-bin non-closure is disclosed rather than silently changing either source. The representative
set is used as the frozen ANI95 product; any future rebuild must identify and adjudicate the missing membership
before replacing it. The exact current-only identifier is `25PRr2__self_SemiBin_2_sub`; the reconciliation,
source paths and SHA-256 digests are frozen in `reproducibility/drep_membership_reconciliation.tsv`.

GTDB-Tk 2.7.2/R232 taxonomy and higher-rank novelty claims refer to the 5,438 representatives. The active R232
census assigns all representatives through class; lower-rank gaps are reported separately. An unplaced result is
not, by itself, a new higher-order taxon.

## Historical binning and functional snapshots

Two SemiBin2 strategies were run and annotated in the processed tree:

| Historical strategy | HQ+MQ MAGs | All bins including LQ | Driver list |
|---|---:|---:|---|
| `single_sample` pretrained model | 15,469 | 45,869 | `M/03_annotation/mag_list.txt` |
| `self_supervised` model | 21,616 | 86,079 | `M/03_annotation/mag_list_self_supervised.txt` |

A raw-union dRep95 comparison of 37,085 HQ/MQ bins (`M/04_drep_self_vs_single/`) flagged 3,386
self-supervised bins as nonredundant with the single-sample set. This is a historical cross-binner comparison, not
the active 5,438-representative catalogue and not evidence for 3,386 novel taxa.

Per-bin Bakta and EggNOG outputs for both strategies live under:

```text
P/<trip>/<site>/<sample>/wgs/annotation/
  single_sample/{bakta,eggnog}/SemiBin_<N>/
  self_supervised/{bakta,eggnog}/SemiBin_<N>/
```

Historical GTDB-Tk and CheckM2 outputs also exist for both strategies. Their presence shows that the bins were
processed; it does not make their release, taxonomy or catalogue membership current.

## Protein universes

| Protein universe | Count | Scope and manuscript role |
|---|---:|---|
| historical MAG-protein snapshot | 51,920,778 | Bakta proteins from 15,469 single-sample MAGs; source of the historical 51.9 M → 23.4 M → 2.84 M sequence-dark ladder and several explicitly labelled functional censuses |
| recurrent sequence-dark orthogroups | 23,382 | pyrodigal calls clustered across the frozen 5,438 representatives; used for ecology/specificity analyses and kept distinct from the historical ladder families |
| assembly-level catalogue | 753,305,533 | all predicted proteins on binned and unbinned contigs across 552 assemblies; a separate “annotate everything” arm, not a denominator for current manuscript percentages |

The assembly-level catalogue is stored on IBEX at `M/04_protein_catalog/`. Its ProstT5/3Di/FoldSeek processing
uses GPU jobs on Orix under `/mnt/data/u/hohndor/eq-prostt5/`. Those outputs must be promoted through a manifest
before they can replace any manuscript number.

## Interpretation and update rules

- Label every statistic with its catalogue and database release; do not combine the 15,469 historical functional
  denominator with the 5,438 current representative denominator.
- Treat the 753 M assembly catalogue as a separate discovery universe until its membership, searches and output
  tables enter the checksum-closed claim ledger.
- A coassembly or LQ-bin candidate may motivate recovery work, but higher-rank novelty requires source-read
  support, cellularity/quality checks, dereplication context, GTDB placement and de novo phylogeny/RED adjudication.
- Rebuilding the representative set or rederiving all functions on it is a new analysis version. Preserve the old
  universe and hashes rather than overwriting historical evidence.
