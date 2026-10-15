# Full Biological Discovery Census Report: Tracks F3, F4, and F5

**Date:** 2026-10-15 07:36:34 UTC  
**Execution Status:** COMPLETED_AUDIT_VERIFIED (Strict Provenance, Zero Guesswork)  
**Cohorts Audited:**
1. **Empty Quarter Reps:** 5,438 MAGs (20,451,986 proteins from Root full cohort `cohort.faa`)
2. **SMAG Harmonized Environmental Cohort:** `smag_pool_canonical.faa` (4,210,592 proteins)
3. **GEM Harmonized Environmental Cohort:** `gem_harmonized_pool_canonical_seed20261014.faa` (4,210,592 proteins)

---

## 1. Executive Summary Table

| Track | Family / Class | Total Candidates | Complete CDS (`partial=00`) | EQ MAGs | SMAG Envr. | GEM Envr. | Unique Seq SHA256 |
|---|---|---|---|---|---|---|---|
| **F3** | Type-1 Microbial Rhodopsin | 319 | 312 | 144 | 12 | 163 | 305 |
| **F3** | Heliorhodopsin Type-3 | 186 | 178 | 67 | 39 | 80 | 183 |
| **F4** | Ectoine Synthase ($EctC$) | 1085 | 1059 | 966 | 56 | 63 | 1045 |
| **F5** | $\\gamma$-Carbonic Anhydrase | 6668 | 6402 | 4381 | 1051 | 1236 | 6398 |
| **F5** | $\\beta$-Carbonic Anhydrase | 6351 | 6034 | 4301 | 1011 | 1039 | 6122 |
| **F5** | Carboxysome Shell CA ($CsoSCA$) | 59 | 48 | 42 | 3 | 14 | 54 |
| **Total** | **All Families** | **14668** | **14033** | **9901** | **2172** | **2595** | **14107** |

---

## 2. Track F4 Synteny & Pathway Context

Granular analysis of $EctC$ genomic neighborhood (same-contig $\\le 5,000$ bp distance and strand orientation):

| Pathway Synteny Organization | Candidate Count | Percentage | Description |
|---|---|---|---|
| `Complete_ABCD_Operon` | 499 | 46.0% | Co-localized $EctA + EctB + EctC + EctD$ on same contig |
| `Complete_ABC_Operon` | 256 | 23.6% | Co-localized $EctA + EctB + EctC$ without nearby $EctD$ |
| `Partial_Operon_BCD` | 51 | 4.7% | Co-localized $EctB + EctC + EctD$ (missing nearby $EctA$) |
| `Partial_Operon_BC` | 30 | 2.8% | Co-localized $EctB + EctC$ (missing nearby $EctA$) |
| `Partial_Operon_AC` | 8 | 0.7% | Co-localized $EctA + EctC$ (missing nearby $EctB$) |
| `Partial_Operon_CD` | 38 | 3.5% | Co-localized $EctC + EctD$ (missing nearby $EctA/B$) |
| `Solitary_sameORF_EctCD_fusion` | 2 | 0.2% | Single ORF encoding non-overlapping EctC + EctD domains (NOT an operon) |
| `Isolated_EctC_Candidate` | 197 | 18.2% | No pathway markers within 5,000 bp on same contig |

**Assembly Limitation & Completeness Context:**
- **Zero Self-Gene Operon Artifacts:** Fixed self-exclusion ($mpid \\neq pid$) in neighbor scanning; self-hits are strictly prevented from declaring false operons.
- **Profile-Supported Same-ORF EctC+EctD Fusions:** Detected via non-overlapping domain coordinate intervals ($\min(c_{to}, d_{to}) < \max(c_{from}, d_{from})$). Pilot example `V60PRr2__single_SemiBin_127_sub::V60PRr2_211258_4` confirmed as solitary fusion (EctC 1..125, EctD 147..422).
- **CheckM Absence is Conditional:** Completeness estimates cannot prove absence. Candidates in MAGs with CheckM $\ge 90\\%$ are classified as conditional absence (`Unlinked_EctC_in_high_completeness_MAG_absence_conditional`), never claiming "true complete genome" or "true orphan".
- **External Cohorts:** SMAG and GEM cohorts record completeness as `NA` (`External_Cohort_CheckM_Not_Loaded_cannot_exclude_unlinked_pathway`).

---

## 3. Sequence Deduplication & Cross-Cohort Sharing Matrix

| Cohort Sharing Category | Unique Sequences (SHA256) | Notes |
|---|---|---|
| **EQ Rep MAGs Only** | 9775 | Endemic to Empty Quarter MAG cohort |
| **SMAG Environmental Only** | 2066 | Found only in SMAG cohort |
| **GEM Environmental Only** | 2241 | Found only in GEM cohort |
| **Shared: EQ MAG + SMAG** | 2 | Exact identical 100% protein sequence |
| **Shared: EQ MAG + GEM** | 2 | Exact identical 100% protein sequence |
| **Shared: SMAG + GEM** | 21 | Environmental cross-pool overlap |
| **Shared Across All Three** | 0 | Ubiquitous across EQ, SMAG, and GEM |

*Cross-source memberships are preserved; shared sequences are never arbitrarily claimed as EQ-only.*

---

## 4. ProstT5 Input Packages for Downstream Structure Embedding

Primary packages are strictly restricted to verified complete candidates (`partial=00`, `is_complete_gene=1`) and source-checked curated references. Incomplete candidates are retained in separate partial tables and never silently dropped:

| Family Package | Curated Refs | Complete Candidates | Retained Partials | Primary FASTA |
|---|---|---|---|---|
| **Type-1 Microbial Rhodopsin** | 30 | 298 | 7 | `prostt5_type-1_microbial_rhodopsin.faa` |
| **Heliorhodopsin Type-3** | 26 | 175 | 8 | `prostt5_heliorhodopsin_type-3.faa` |
| **Ectoine Synthase ($EctC$)** | 30 | 1019 | 26 | `prostt5_ectc.faa` |
| **$\\gamma$-Carbonic Anhydrase** | 28 | 6132 | 266 | `prostt5_gamma_ca.faa` |
| **$\\beta$-Carbonic Anhydrase** | 12 | 5805 | 317 | `prostt5_beta_ca.faa` |
| **Carboxysome Shell CA ($CsoSCA$)** | 3 | 43 | 11 | `prostt5_csosca.faa` |

---

## 5. Strict Negative Assertions & Scientific Guardrails
1. **Biosynthetic Context $\\neq$ In Vitro Activity:** Co-localized $EctABC$ / $EctABCD$ gene neighborhoods predict biosynthetic pathway capacity based on domain homology, not confirmed metabolic flux or biochemical turnover.
2. **Absence is Conditional:** Even in high-completeness MAGs ($\ge 90\\%$) or CheckM 100% estimates, statistical absence cannot prove complete physical absence of pathway genes. Never claim "true complete genome" or "true orphan".
3. **Sensory Rhodopsins $\\neq$ Novel Phototrophy:** Detection of microbial and heliorhodopsin domain homologs reflects phylogenetic sequence diversity and putative sensory/ion-translocating capacity, not novel phototrophic energy generation.
4. **Enzymatic Activity Plausibility Only:** $\\gamma$-Carbonic anhydrase identification establishes sequence and fold homology; active-site coordination supports plausibility only and does not establish enzymatic catalysis.
5. **Strict Pre-Execution Code Immutability:** All executed scripts were frozen to `results/full_v1/frozen_code/` before execution and verified by SHA256 checksums.
