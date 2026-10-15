# Current literature synthesis and claim decisions

Date: 2026-08-03

Scope: the MAG/taxonomic-novelty, protein-dark-matter, and structural-novelty
parts of the Empty Quarter manuscript. This supplements the 2026-07-29
internal reviews and uses literature available through 2026-08-03.

## Executive decision

The strongest defensible discovery programme is not a raw catalogue-size or
raw no-hit story. It is a prospective, falsifiable recovery-and-calibration
story:

1. recover candidate lineages with candidate-targeted sample selection and
   independent split reconstructions;
2. require GUNC/marker/coverage coherence and stable marker-jackknife
   placement before calling a high-rank lineage;
3. describe protein darkness relative to named, versioned references and
   separate coding support from function;
4. call current structure-reference survivors *calibrated candidate folds*,
   not new folds, until experimental validation.

The current catalogue already supports substantial novelty at order and below,
but no new phylum or class. A new high-rank result can only come from the
prospective recovery arm.

## 1. Coassembly can recover deep novelty, but generic pooling is not enough

[Bin Chicken](https://www.nature.com/articles/s41592-025-02901-1) is the most
directly relevant recent result. It selects metagenomes for coassembly using
shared marker sequences in raw reads, can prioritize markers divergent from
references, and recovered first genomic representatives of six phyla and 41
classes. This supports Robert's premise that deeper Empty Quarter lineages may
be recoverable when evidence is distributed across samples.

Two recent benchmarks impose important constraints:

- A 2025 real-data benchmark found that short-read coassembly recovered fewer
  MQ/near-complete/HQ MAGs than single- or multi-sample strategies and warned
  that coassembly can create inter-sample chimeric contigs and erase
  sample-specific variation
  ([Han et al. 2025](https://www.nature.com/articles/s41467-025-57957-6)).
- A 2026 benchmark found wide tool-to-tool variation in chimerism, noted that
  CheckM2 does not reliably detect non-redundant contamination, and found that
  multi-sample binning was best around 20 related samples, with both smaller
  and larger sets sometimes performing worse
  ([Kim et al. 2026](https://www.nature.com/articles/s41467-026-71521-w)).

### Consequence for the running analysis

The current 126 site-by-compartment assemblies remain a useful broad rescue
screen, but they are not by themselves a phylum-discovery design. After the
screen completes, the high-rank follow-up should use candidate-marker/read
recruitment to select a bounded set of samples, cap the primary multi-sample
arm near 20 when enough positive samples exist, and retain ecological grouping
as a sensitivity arm. Every candidate must be reconstructed in two
non-overlapping sample sets without seeding bins from the discovery assembly.
COMEBin is justified as one high-performing neural binner, but independent
SemiBin2/VAMB or ensemble agreement and GUNC are required because binner choice
changes chimerism.

This design can produce either a real discovery or a strong negative. A branch
that disappears under split reconstruction or marker jackknifing is an
assembly/binning artefact, not a borderline new phylum.

## 2. An unassigned rank is not a taxon, and a taxon is not yet a name

The modern phylum concept is phylogenetic and GTDB rank normalization is based
on concatenated marker trees and relative evolutionary divergence, not on an
empty classifier field
([Pallen 2024](https://doi.org/10.1099/ijsem.0.006508)). GTDB itself uses ANI
for species and RED on marker-gene trees for higher ranks
([GTDB release 10](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807784/)).

If naming is eventually warranted, the current
[SeqCode](https://registry.seqco.de/page/seqcode) requires an INSDC-deposited
type sequence; its type-MAG quality floor is greater than 90% completeness and
less than 5% contamination, with genome/16S agreement. These are naming
requirements, not sufficient evidence that a branch deserves phylum rank.

### Consequence for the manuscript

The R232 census should report the operational rank gaps exactly: no
domain/phylum/class gap among 5,438 HQ/MQ representatives, then 11 order, 142
family, 1,343 genus, and 3,755 species gaps. The 11 order gaps are candidates
for placement, not 11 new orders. A high-rank claim requires a coherent,
independently recovered HQ genome, stable concatenated-marker placement,
rank-normalized divergence, contamination/strain-mixture exclusion, and a
clearly separated nomenclatural step.

## 3. Protein darkness remains large, but short ORFs need positive evidence

The 2025 microbial-dark-matter review emphasizes that MAGs/SAGs expose large
uncultured functional space but do not functionally characterize it
([Iqbal et al.](https://pubmed.ncbi.nlm.nih.gov/38385313/)). Global primary
studies make the comparator much harder than UniRef alone. NMPF reported
106,198 large reference/Pfam-negative families across global metagenomes
([Pavlopoulos et al. 2023](https://www.nature.com/articles/s41586-023-06622-3)),
and GMSC contains roughly 965 million non-redundant small ORFs across 63,410
metagenomes. Crucially, only 43.6 million GMSC entries (4.5%) passed all
in-silico quality tests plus transcriptional or translational evidence, and
short-sequence homology recovery deteriorates sharply as sequences approach
20 aa
([Duan et al. 2024](https://www.nature.com/articles/s41467-024-51894-6)).

The very recent microbial-unknome synthesis likewise argues for expression,
spatial/coexpression context, and reverse genetics rather than treating an
unknown label as a function
([Vanni et al. 2026](https://www.nature.com/articles/s41437-026-00870-5)).

### Consequence for the manuscript

The current GMSC refresh and matched genomic-decoy work are essential rather
than optional. Recurrence across independently assembled samples, intact
coding boundaries, AntiFam/low-complexity checks, and codon-level constraint
support a predicted coding ORF; none proves translation or function. MGnify,
GMSC, UniRef and structure searches must retain separate versioned
denominators. “No homolog anywhere” is not supportable.

The corrected B2 dN/dS result now fits this literature: it is bounded to
22,767 recurrent, alignable dark orthogroups and is consistent with purifying
selection. It is not a calibrated test of all predicted dark ORFs.

## 4. Structural novelty is now a calibration problem

AFESM combines more than 800 million AFDB/ESM Atlas models. Its current
versions find few genuinely new domain folds relative to the immense scale,
while discovering about 11,941 previously unseen multi-domain architectures
([Yeo et al.](https://doi.org/10.1101/2025.04.23.650224)). The precise novel
fold count changed materially with re-prediction and preprint version, which is
itself evidence that model quality and version freezing are decisive. A 2026
review concludes that fold space appears near saturation while remote
homology, architecture and functional assignment remain open
([Sahakyan et al. 2026](https://doi.org/10.1002/pro.70397)).

### Consequence for the running analysis

PDB/AFDB/MGnify triage before Orix folding is correct. The current historical
cohort result—608 dual current-reference no-hits among 1,875 length- and
confidence-qualified models—must remain a *candidate* result and be shown with
the 101/3,013 known-fold false-no-hit calibration. The recurrent-family gate
must finish before freezing the residual Orix manifest. Orix should fold only
that checksum-closed residue, with GPU utilization monitored; no raw
“new-fold” denominator should be promoted.

## Manuscript promotion checklist

- Replace all R226/15,469-MAG census text with the R232/14,426-bin and
  5,438-representative universes only after the final generated tables exist.
- Make lower-rank novelty the current taxonomic result; make high-rank recovery
  a prospective, separately validated result.
- Add the 2025 and 2026 binning benchmarks to Methods/Discussion when the
  recovery output is promoted.
- Report reference-relative no-hit counts with database release, threshold,
  query universe and length/coding-quality stratum.
- Prefer “candidate structural novelty” or “structure-reference survivor” to
  “novel fold”. Report the known-fold miss curve beside every candidate count.
- Treat novel domain architecture, conserved genomic neighborhood and
  experimentally testable function as more plausible near-term headlines than
  a large uncalibrated fold count.
