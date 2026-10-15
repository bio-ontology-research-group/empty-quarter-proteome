# DESIGN.md — Hypothesis 4: Within-Model Recovery of Confident Structural Regions Under Reference-Conditional Cropped Search

**Analysis Directory**: `reproducibility/darkmatter/round6_hypotheses/h4_structure`  
**Remote Workbench**: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h4_structure`  
**Author**: Gemini Flash High (acting for Robert Hoehndorf)  
**Date**: 2026-10-13  
**Status**: Revised Post-Audit Specification (V2)  
**Superseded Version**: `scripts/v1_superseded/`, `results/v1_superseded/`  

---

## 1. Scientific Context and Scope Narrowing

### 1.1 Context in Current Manuscript
In the Empty Quarter MAGs dark-matter analysis (§1.3, lines 307–351), Foldseek structural searching of 2,838 confidence-qualified representative models ($\ge 100$ aa, mean pLDDT $\ge 0.70$) against the Protein Data Bank (PDB) and AlphaFold Database (AFDB) identified:
- **1,031 dual no-hit models ($K = 1,031$, 36.3%)** under query-normalised TM score ($qTM < 0.5$).
- In contrast, under alignment-normalised TM score ($alnTM < 0.5$), only **95 models (3.35%)** were dual no-hits.
- Co-searched with them were 3,013 gate-passing PDB-derived positive controls ($KF...$), reporting an apparent error floor of 3.35% (101/3,013) under $qTM < 0.5$ and 2.26% (68/3,013) under $alnTM < 0.5$.

The current manuscript already cautions:
> *"The query and alignment scores are independent maxima across reported hits and may come from different reference targets... Domain arrangement, extensions and uncertain inter-domain placement can affect this comparison. Confidence-qualified predictions alone do not validate a novel fold (orr2026spurious)... The controls derive from PDB chains present in the references and serve as positive retrieval controls. Homology-excluded controls and domain-resolved analysis would be needed to calibrate remote homology and delimit unmatched structural regions."*

### 1.2 Scope Narrowing from Original "Domain" Objective
Root's audit (`reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`) identified that earlier exploratory drafts overreached by:
1. Labeling confidence-thresholded residue runs "bona fide novel structural domains".
2. Claiming biological disorder, flexible linkers, and quaternary reorientation without physical/biochemical testing.
3. Estimating region coverage by overlapping whole-chain hits rather than actually isolating and searching the regions.
4. Calling 4-letter source-PDB identifier exclusion "homology-free calibration".

**Explicit Narrowing**:
- Contiguous residue segments passing confidence criteria ($\ge 50$ aa, smoothed pLDDT $\ge 0.70$) are operational **confident structural regions**, not confirmed evolutionary or functional domains.
- Low pLDDT is an index of algorithmic model uncertainty and low predicted alignment error, **not established biological intrinsic disorder**.
- The scientific objective is strictly narrowed to: **testing whether physical cropping of outcome-independent confident structural regions recovers structural matches ($qTM_{crop} \ge 0.5$) in PDB among models that previously scored $qTM_{wc} < 0.5$ at the whole-chain level**, evaluated under multi-tier sequence identity exclusion controls and an internal fold-matched dark comparator.

---

## 2. Identifiable Hypothesis and Estimands

### 2.1 Formal Hypotheses
- **Alternative Hypothesis ($H_1$)**: Among confidence-qualified dark models exhibiting whole-query divergence ($qTM_{wc} < 0.5$) but possessing local structural hits ($alnTM_{wc} \ge 0.5$; stratum $K_{local}$), physical cropping and independent Foldseek search of confident regions ($\ge 50$ aa, smoothed pLDDT $\ge 0.70$) recovers structural matches ($qTM_{crop} \ge 0.5$) against PDB in a substantial fraction of models ($\ge 50\%$), demonstrating that their whole-chain $qTM < 0.5$ was driven by query-length normalization dilution over multi-region architectures or flanking segments. In contrast, models lacking any whole-chain hit ($K_{dual\_none}$) exhibit significantly lower region recovery ($\le 35\%$) and elevated zero-retrieval rates.
- **Null Hypothesis ($H_0$)**: Within-model region recovery after cropping among $K_{local}$ models does not exceed that of $K_{dual\_none}$ models, or fails to exceed the recovery rate expected when regions remain fundamentally unmatched in PDB ($qTM_{crop} < 0.5$).

### 2.2 Primary Statistical Unit and Estimands
- **Primary Statistical Unit**: The representative protein model (and its corresponding fold family cluster). Region-level evaluations are aggregated to the parent model to account for intra-protein regional dependence.
- **Outcome-Blind Boundary Definition**:
  - Contiguous segments of length $\ge 50$ aa with 15-aa rolling-window mean pLDDT $\ge 0.70$.
  - Gaps $\le 5$ aa below threshold are bridged to prevent fragmenting single secondary structure elements.
  - Boundaries are defined strictly from the intrinsic per-residue pLDDT vector of the predicted coordinate file, with zero knowledge of Foldseek search hits.
- **Primary Estimands**:
  1. **Within-Model Region Recovery Rate ($P_{rec}$)**:
     $$P_{rec} = \mathbb{P}\left(\exists r \in \mathcal{R}(M) : qTM_{crop}(r) \ge 0.5 \;\middle|\; qTM_{wc}(M) < 0.5\right)$$
     where $\mathcal{R}(M)$ is the set of confident cropped regions in model $M$.
  2. **Model-Level Outcome Decomposition**:
     - *Recovered after Cropping*: $qTM_{wc} < 0.5$ and $\max_{r} qTM_{crop}(r) \ge 0.5$.
     - *Persistent Low Score*: $qTM_{wc} < 0.5$ and $\max_{r} qTM_{crop}(r) < 0.5$ (with at least one region retrieved).
     - *Zero Retrieval*: No hits returned by Foldseek for any cropped region (or model possesses 0 confident regions $\ge 50$ aa).
     - *Whole-Chain Matched*: $qTM_{wc} \ge 0.5$ (baseline status for controls and fold-matched dark models).
  3. **Multi-Tier Sequence-Identity Exclusion Sensitivity**:
     For each cropped region, evaluate top PDB hit across exclusion tiers:
     - **Tier 0**: Unfiltered top hit.
     - **Tier 1 (Source-PDB Excluded)**: Exclude hits where the 4-letter target PDB code matches the control's known source PDB accession (`pdb_structure_id` in manifest).
     - **Tier 2 (SeqId < 30%)**: Exclude hits with sequence identity $\ge 30\%$ ($\text{fident} \ge 0.30$).
     - **Tier 3 (SeqId < 25%)**: Exclude hits with sequence identity $\ge 25\%$ ($\text{fident} \ge 0.25$).
     - **Tier 4 (SeqId < 20%)**: Exclude hits with sequence identity $\ge 20\%$ ($\text{fident} \ge 0.20$).

---

## 3. Study Design and Stratified Cohorts

### 3.1 Stratified Analysis Cohorts
To evaluate recovery and control sensitivity, the study stratifies the confidence-qualified population into four mutually exclusive groups:
1. **$K_{local}$ ($N = 100$ pilot; full population $n = 936$)**: Candidates with whole-chain $qTM < 0.5$ but local whole-chain hit ($alnTM \ge 0.5$).
2. **$K_{dual\_none}$ ($N = 50$ pilot; full population $n = 95$)**: Candidates with whole-chain $qTM < 0.5$ and no whole-chain hit in PDB or AFDB.
3. **$Q \setminus K$ ($N = 75$ pilot; full population $n = 1,807$)**: Internal fold-matched dark comparator models with whole-chain $qTM \ge 0.5$.
4. **$Controls$ ($N = 75$ pilot; full population $n = 3,013$)**: Blinded known-fold PDB-derived controls passing confidence qualification.

### 3.2 Two-Stage Execution Strategy
1. **Stage 1 (Completed Full-Cohort Audit, $N = 5,851$)**:
   - Audit all 2,838 dark models and 3,013 controls from raw hit files (`pdb_hits.raw.tsv`, `afdb_hits.raw.tsv`).
   - Quantify disconnected maxima between $qTM_{max}$ and $alnTM_{max}$.
   - Separate zero retrieval from score failure across the entire population.
   - Quantify exact source-PDB reference leakage in controls.
2. **Stage 2 (Completed Cropped Search Pilot, $N = 300$ models, $413$ regions)**:
   - Deterministic sample (seed `20261013`) across all four strata.
   - Physically crop PDB models into individual PDB files (`crops/{model_id}_reg{reg_id}.pdb`), preserving exact residue numbering, coordinates, and chain IDs.
   - Execute Foldseek search (`easy-search --alignment-type 1 -e 10 -s 9.5 --max-seqs 50`) directly on cropped regions against PDB (`foldseek_pdb`).
   - Compute single-hit exact metrics ($qTM_{crop}, alnTM_{crop}, \text{fident}, \text{alnlen}, \text{qcov}$).
   - Evaluate multi-tier sequence-identity sensitivity and within-model recovery.

---

## 4. Threats to Validity and Explicit Mitigations

| Threat to Validity | Nature of Bias | Methodological Mitigation |
|---|---|---|
| **Disconnected Maxima** | Taking $qTM_{max}$ and $alnTM_{max}$ across different targets distorts coverage | Hits evaluated strictly per-row; region-level $qTM_{crop}$ and $alnTM_{crop}$ refer to the identical aligned target. |
| **Circular Boundaries** | Defining regions from search hits inflates apparent coverage | Boundaries defined strictly from per-residue pLDDT profiles prior to and independent of search. |
| **Overlapping vs Cropped Search** | Whole-chain search coordinates do not reflect region-isolated alignment behavior | Residues are physically cropped into independent PDBs and searched de novo with Foldseek. |
| **Control Reference Leakage** | Controls derived from PDB structures hit their own crystal structure (94.4% self-hit rate) | Tier 1 excludes exact 4-letter source PDB IDs; Tiers 2–4 impose sequence identity thresholds (<30%, <25%, <20%). |
| **Residual Homology in Controls** | PDB contains close paralogs and homologous families not removed by 4-letter ID | Evaluated across decreasing sequence identity tiers down to 20% to document sensitivity and lower bounds. |
| **Search Truncation** | Foldseek heuristic prefilters (`-s 9.5`, `--max-seqs 50`) may miss remote structural hits | Zero retrieval is tracked separately from score divergence; retrieval rate is explicitly reported for every stratum. |
| **Conflating Low pLDDT with Disorder** | Unfolded or poorly predicted regions are not confirmed biological disorder | Labeled strictly as "low-confidence regions"; no claims of intrinsic disorder or flexible linkers are made. |
| **Claims of Novel Folds** | Persistent low score ($qTM < 0.5$) does not prove fold novelty | Labeled strictly as "persistent structural divergence under reference-conditional Foldseek search"; no fold novelty claimed. |

---

## 5. Success Criteria and Informative Negative Boundaries

### 5.1 Criteria for Methodological Success
1. **Verification of Cropping Pipeline**: Successful physical extraction and Foldseek search of $\ge 400$ cropped regions with zero corrupted structures.
2. **Separation of Outcomes**: Complete resolution into zero retrieval, persistent low score, and recovered match for all models.
3. **Control Sensitivity Quantification**: Precise quantification of how sequence identity filtering affects control match rates relative to dark candidates.

### 5.2 Informative Negative Boundaries (Falsification Criteria)
- If within-model region recovery in $K_{local}$ is $\le 20\%$ (indistinguishable from $K_{dual\_none}$), the hypothesis that whole-chain divergence is driven by query-length normalization dilution over multi-region architectures is refuted.
- If control match rates under sequence identity $< 25\%$ fall to or below the candidate match rate, the apparent divergence between candidates and controls is attributable to database representation bias rather than structural uniqueness.
