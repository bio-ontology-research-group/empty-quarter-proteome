# DESIGN: Conserved Gene Neighbourhoods of Low-Query-TM vs Fold-Matched Sequence-Dark Families

**Hypothesis 2:** Conserved gene neighbourhoods of low-query-TM versus fold-matched sequence-dark families.  
**Author:** Gemini Flash High (acting on behalf of Robert Hoehndorf)  
**Date:** 2026-10-13  
**Status:** Pre-computation specification (mandatory prior to computation)  
**Local Workspace:** `reproducibility/darkmatter/round6_hypotheses/h2_context/`  
**Remote Workbench:** `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h2_context/`

---

## 1. Scientific Aims & Exact Hypothesis

### 1.1 Context and Scope
In the Empty Quarter metagenomic catalogue, thousands of gene families lack detectable sequence homology to known reference databases (sequence-dark). A subset of these dark families have confidence-qualified ESMFold structure predictions (mean pLDDT $\geq 0.70$, length $\geq 100$ aa; $n = 2{,}838$). In the manuscript (§2.4), these structural models split into two arms:
1. **Low-query-TM arm ($K$):** 1,031 families whose representatives have query-normalised TM-scores $< 0.5$ against both PDB and AFDB (dual no-hits; candidate structurally novel folds).
2. **Fold-matched sequence-dark arm ($Q \setminus K$):** 1,807 families whose representatives have query-normalised TM-scores $\geq 0.5$ against at least one reference (familiar/matched fold despite lacking sequence similarity).

The biological function of low-query-TM families is completely unknown. Genomic context—specifically whether a gene consistently co-occurs beside specific functionally annotated neighbours across phylogenetically distinct species—can generate bounded, falsifiable functional hypotheses without requiring ungrounded expression or wet-lab claims.

### 1.2 The Formal Hypothesis
* **Research Question:** Do low-query-TM families repeatedly occur beside the same annotated neighbours across distinct species, compared with fold-matched sequence-dark families, after controlling for sequence length, member count, species coverage, contig fragmentation, and host lineage?
* **Null Hypothesis ($H_0$):** The rate of conserved, multi-species functional gene neighbourhood embedding does not differ between low-query-TM dark families and fold-matched sequence-dark families ($P(C_{\text{fam}} = 1 \mid K) = P(C_{\text{fam}} = 1 \mid Q \setminus K)$) after adjusting for family size, species breadth, contig edge truncation, and taxonomic lineage.
* **Alternative Hypothesis ($H_1$):** Low-query-TM dark families exhibit differential gene neighbourhood conservation relative to fold-matched sequence-dark families—either depleted (e.g., if concentrated in fast-evolving, mobile, or orphan loci) or enriched in specific conserved cellular metabolic/repair operons.

---

## 2. Primary Estimand & Unit of Analysis

### 2.1 Unit of Analysis
The **gene family** is the primary unit of analysis ($n = 1{,}031$ for $K$; $n = 1{,}807$ for $Q \setminus K$).
* *Anti-pseudoreplication rule:* Member genes within a family cannot be treated as independent observations because they share evolutionary history and genomic context. All flanking gene observations are aggregated to family-level summary statistics.

### 2.2 Population and Admission Criteria
A family is admitted to the inferential test if and only if it satisfies:
1. **Confidence-qualified structure model:** mean pLDDT $\geq 0.70$, length $\geq 100$ aa in `04_structure/00_structure_clean/05_validate/per_rep_structures.tsv`.
2. **Multi-species presence:** Admitted members reside on contigs binned into at least two distinct species clusters (dRep 95% ANI secondary clusters from `Cdb.csv`).
3. **Callable flanks:** At least two distinct species have non-truncated contig flanks (the focal gene does not terminate at a contig edge).

### 2.3 Primary Estimand
Let $F$ be an admitted family with members across $S_F \geq 2$ distinct species.
For each member $m \in F$ on contig $c$:
* Define the flanking window $W = \pm 2$ coding sequences (CDSs) immediately upstream and downstream along the contig.
* Record the annotation of each flank $j \in \{-2, -1, +1, +2\}$:
  * eggNOG KEGG Orthologue ($KO$)
  * COG functional category
  * Bakta product string
  * Sequence darkness indicator (via 90% cluster representative in `dark_rep_ids.txt.gz`)
  * Truncation status (flagged if contig edge terminates within window)
  * Orientation relative to focal gene (same strand vs opposite strand)

**Family-level Conserved Neighbourhood Indicator ($C_{\text{fam}}$):**
$$C_{\text{fam}} = \mathbb{I}\left(\exists \, k \in \text{Annotated KOs} \text{ such that } k \text{ is present in } W \text{ across } \geq 2 \text{ distinct species}\right)$$

**Continuous Neighbourhood Recurrence Score ($R_{\text{fam}}$):**
$$R_{\text{fam}} = \max_{k \in \text{KOs}} \frac{\text{Number of distinct species with } k \text{ in flank } W}{S_F}$$

The primary estimand is the risk difference:
$$\Delta = \mathbb{E}[C_{\text{fam}} \mid K] - \mathbb{E}[C_{\text{fam}} \mid Q \setminus K]$$
and the odds ratio adjusted for covariates via logistic regression.

---

## 3. Controls, Confounds & Boundary Enforcement

### 3.1 Essential Confounds and Mitigation
1. **Family Prevalence & Species Breadth:** Families with more members across more species have more combinatorial opportunities to observe recurrent neighbours by chance.
   * *Mitigation:* Match or condition on $S_F$ (number of distinct species) and member count $N_F$.
2. **Contig Fragmentation and Edge Truncation:** Metagenomic contigs frequently terminate within a few hundred base pairs of a gene. A missing neighbour caused by a contig edge must NOT be conflated with biological absence of an operon.
   * *Mitigation:* Explicitly record contig length, distance to 5' edge, and distance to 3' edge. Flag any window truncated by a contig boundary as `edge_truncated`. Exclude truncated flanks from the denominator of callable context.
3. **Annotation Ascertainment Bias:** A neighbour may be a real protein that simply lacks a KO assignment (e.g., hypothetical protein). Missing annotation must NOT be coded as absence of gene context.
   * *Mitigation:* Separate unannotated/hypothetical neighbours from annotated neighbours. Evaluate both KO recurrence and 90% cluster recurrence.
4. **Strand and Spacing:** Same-strand proximity ($<100$ bp intergenic gap) is consistent with co-transcription (polycistronic operons), but does NOT prove operon expression or shared regulation without transcriptomic data.
   * *Mitigation:* Report strand orientation and intergenic distance as structural context descriptors, never as proof of an operon or enzymatic pathway.
5. **Mobilome / Phage Boundary:** Transposases, integrases, and phage structural proteins frequently flank inserted genetic elements and can exhibit recurrent context across unrelated taxa.
   * *Mitigation:* Flag neighbours annotated with mobilome/phage Pfams/KOs (IS elements, transposases, integrases, terminases, portal proteins). Respect companion phage/mobilome repository boundaries; do not reclassify mobile elements as novel cellular pathways.
6. **No Invented Expression Cohort:** No RNA sequencing or expression validation exists in this paper. Reject any language implying 136 expressed targets or transcriptomic confirmation.

### 3.2 Statistical Null Model
To determine whether observed neighbourhood recurrence exceeds background genomic expectation:
* **Within-lineage permutation null:** For each contig harbouring a focal family member, permute the position of the focal gene among the CDS positions of the contig, or shuffle gene labels across contigs belonging to the same GTDB order/class while preserving contig length, gene count, strand distribution, and annotation status.
* Compute $C_{\text{fam}}^{\text{null}}$ over 1,000 permutations to generate an empirical null distribution and two-sided p-value.

---

## 4. Upstream Data Architecture & Join Verification

All data joins are anchored in verified, immutable project tables:

```
[per_rep_structures.tsv] (04_structure/00_structure_clean/05_validate/)
   |-- id (modelled representative)
   |-- confidence_qualified == 1 (2,838 models)
   |-- conf_qualified_dual_novel == 1 (1,031 K) vs == 0 (1,807 QnotK)
   v
[fam_reps.tsv] (04_structure/00_structure_clean/05_validate/)
   |-- chosen_rep -> family
   v
[clusters.tsv] (04_structure/00_structure_clean/04_cluster/)
   |-- family -> member proteins (4,880 in K; 7,676 in QnotK)
   v
[assembly_proteins/{sample}/{sample}.faa.gz] (04_protein_catalog/)
   |-- Pyrodigal headers: >sample::contig_gene # start # end # strand # partial=...
   v
[contig_to_bin_refined.tsv.gz] (04_protein_catalog/)
   |-- sample::contig -> refined MAG bin
   v
[Cdb.csv / Wdb.csv] (05_refine/_drep/out/data_tables/)
   |-- MAG bin -> species cluster (secondary_cluster, 95% ANI)
   v
[EmptyQuarter_Data/.../annotation/.../{bakta,eggnog}/]
   |-- Bakta .tsv: contig_N, cds, start, stop, strand, locus_tag, product
   |-- eggNOG .emapper.annotations: locus_tag -> KEGG_ko, COG, Description
```

### Join Audit Results (Verified on IBEX 2026-10-13):
* **Pyrodigal-to-Bakta Coordinates:** Identical 1-based start and stop coordinates verified between `04_protein_catalog/assembly_proteins` and Bakta `.tsv` files on refined MAG contigs.
* **Refined Binning Rate:** 52.9% of K member contigs (2,565 / 4,845) and 47.7% of QnotK contigs (3,627 / 7,602) reside in refined MAGs.
* **Multi-Species Representation:** 608 K families (59.0%) and 960 QnotK families (53.1%) have admitted members binned into $\geq 2$ distinct species clusters, confirming adequate statistical support for multi-species neighbourhood comparison.

---

## 5. Feasibility Pilot & Scaling Strategy

### 5.1 Prototype/Pilot Scope
To establish feasibility and test the complete pipeline under modest resource limits ($\leq 16$ CPUs, $\leq 64$ GB RAM, $\leq 30$ min wallclock on IBEX):
* Sample a deterministic, stratified cohort of 100 families (50 $K$ and 50 $Q \setminus K$) with members in $\geq 2$ distinct species, stratified by member count and phylum.
* Execute end-to-end join:
  1. Extract exact contig coordinates and flanking genes ($\pm 2$ CDSs) from assembly FASTAs.
  2. Map contigs to refined MAGs and dRep species clusters.
  3. Resolve Bakta products and eggNOG KOs for all flanking genes.
  4. Evaluate contig edge truncation.
  5. Compute family-level conserved neighbourhood statistic $C_{\text{fam}}$ and continuous recurrence $R_{\text{fam}}$.
  6. Run within-lineage gene shuffle null (100 permutations).
* Verify pipeline integrity, compute runtimes, and document edge cases.

### 5.2 Go / No-Go Decision Thresholds
* **GO Criterion:**
  - $\geq 90\%$ successful coordinate matching between Pyrodigal headers and Bakta annotations on binned contigs.
  - Callable flanks (non-edge) available for $\geq 50\%$ of member genes.
  - Complete traceability of all input hashes and execution logs.
* **NO-GO / Refusal Criterion:**
  - If coordinate joins fail across callers ($<80\%$ match), or
  - If $>80\%$ of candidate flanks are truncated at contig boundaries, or
  - If annotation availability is $<10\%$, preventing meaningful functional inference.
  Under any NO-GO condition, the hypothesis must be reported as non-viable due to technical limitations rather than claiming biological absence of neighbours.

---

## 6. Structure of Artifacts

* `DESIGN.md` (this file): Pre-computation study design and statistical protocol.
* `INPUTS.tsv`: Pinned paths, sizes, access status, and SHA-256 checksums of all input tables.
* `scripts/`:
  - `01_sample_pilot_cohort.py`: Deterministic cohort selection.
  - `02_extract_neighbourhoods.py`: Flank coordinate extraction and annotation join.
  - `03_analyse_conservation.py`: Family-level recurrence statistics and null permutation.
  - `run_pilot.sbatch`: Slurm submission wrapper for IBEX.
* `results/`: Compact TSV and JSON summaries of pilot joins, recurrence rates, and null comparisons.
* `DOSSIER.md`: Comprehensive evaluation of hypothesis viability, pilot results, confounds, and full execution recipe.
