# Cross-Paper Shared Foundation: Rub' al Khali Multi-Omics

## 1. Overview & Scope

This module establishes the frozen, cross-paper shared foundation across the three companion Rub' al Khali (Empty Quarter) desert manuscripts:
1. **Data Manuscript** (`empty-quarter-data-paper` / `empty-quarter-data-current`): Foundational catalogue, specimen custody, sequencing, metadata representation, and knowledge graph provenance.
2. **Ecology Manuscript** (`empty-quarter-ecology` / `empty-quarter-ecology-current`): Landscape biogeography, environmental/community turnover, and MAG-supported functional potential across landforms and soil depth compartments.
3. **Proteome / Mechanism Manuscript** (`empty-quarter-proteome`): Mechanistic investigations of archaeal adaptation, metabolic requirements and deficiencies, viral host associations, and genomic neighbourhoods of dark exported proteins.

This module unifies the exact genome identifiers, canonical GTDB release 232 taxonomy reconciliation, and metabolic/viral cohort memberships so that all three manuscripts share an identical, verifiable foundation.

---

## 2. Frozen Foundation Cohorts

### A. Total Genome Catalogue (5,438 MAGs)
- **Source**: `inputs/genome_manifest.tsv` (SHA256: `835ded5ac9aee7e50be1405ebf5e9794a4382157b07d30a157d593e94b004bb5`)
- Exactly **5,438 unique MAG identifiers** (`10Dr2__single_SemiBin_0` .. `V60SRr2__single_CONCOCT_10_sub`).
- Zero duplicate identifiers.

### B. Canonical GTDB Release 232 Taxonomy Reconciliation
- **Source**: `inputs/taxonomy_reconciliation_all5438.tsv.gz` (SHA256: `1bcbdbedec9b16e9f901a989778a5f831d2f475ef4c0c3e559a47568c2f9c5a4`)
- Exactly **5,438 rows**, matching the genome manifest 1:1 in both set membership and row order.
- Reconciles earlier analytical taxonomy snapshots against canonical `canonical_gtdb_r232_classification`.
- Domain breakdown:
  - **4,847 Bacteria** (`canonical_domain == "d__Bacteria"`)
  - **591 Archaea** (`canonical_domain == "d__Archaea"`)

### C. Archaeal Partitions & Thermoplasmatota Cohort (237 MAGs)
- **Source**: `inputs/taxonomy_reconciliation_thermo237.tsv` (SHA256: `13214500eea0fc3516b3262a21044344f550c87d4fa7ff3869fc9d958f9aaee8`)
- Within the 591 canonical Archaea:
  - Exactly **237 Thermoplasmatota** (`canonical_phylum == "p__Thermoplasmatota"`).
  - Exactly **354 other Archaea** across seven verified GTDB r232 phyla:
    - `p__Halobacteriota`: 269
    - `p__Thermoproteota`: 56
    - `p__Nanohalarchaeota`: 14
    - `p__Nanobdellota`: 7
    - `p__EX4484-52`: 4
    - `p__Aenigmatarchaeota`: 3
    - `p__Methanobacteriota`: 1
- The 237 Thermoplasmatota MAGs form the focal archaeal cohort for metabolic modeling and surface protein neighborhood analysis, joining 1:1 between `taxonomy_reconciliation_all5438.tsv.gz` and `taxonomy_reconciliation_thermo237.tsv` in ID set, row counts, and full canonical GTDB r232 taxonomy string.

### D. Metabolic 484 Model Inventory & Exact Cohort Join
- **Source**: `inputs/exp1_model_inventory_and_qc_484.tsv` (SHA256: `044014df6aa1f7bf9b5d4fcaea3a2daea0ef62f3dfcb3c79ef08e54b0347eb2b`)
- Exactly **484 unique models**:
  - **237 EQ Desert MAGs** (`cohort == "EQ_desert"`, source `EQ_MAG`): exactly equal to the canonical 237 Thermoplasmatota cohort (1:1 join).
  - **247 External Comparator Genomes** (`cohort == "external_comparator"`): reference genomes with declared inventory source labels from `GTDB_r232_rep` (n=218) and `MGnify` (n=29).
- **Strict Disjointness**:
  - 0 overlap between EQ desert genomes and external comparator genomes.
  - 0 overlap between external comparator genomes and the full 5,438 Empty Quarter MAG catalogue.

---

## 3. Directory Layout & Manifests

```
shared_foundation/
├── README.md                                         # This documentation
├── verify_shared_foundation.py                       # 6-gate standalone lightweight exact-join verifier (<1s)
├── inputs/                                           # Compact frozen source tables (1.6 MB total)
│   ├── exp1_model_inventory_and_qc_484.tsv          # 484 metabolic model inventory (237 EQ + 247 ext)
│   ├── genome_manifest.tsv                           # 5,438 Empty Quarter MAG manifest
│   ├── taxonomy_reconciliation_all5438.tsv.gz        # Canonical GTDB r232 reconciliation for all 5,438 MAGs
│   └── taxonomy_reconciliation_thermo237.tsv         # Canonical 237 Thermoplasmatota focal subset
├── manifests/
│   └── sha256_manifest.tsv                           # Byte-exact SHA256 checksums and file sizes
├── provenance/                                       # Independent root audit proofs
│   ├── exp1_inventory_exact_cohort_join_audit.json   # Prior root proof: 484 inventory cohort join
│   └── shared_taxonomy_input_join_audit.json         # Prior root proof: 5,438 manifest & taxonomy join
└── tests/
    └── test_shared_foundation.py                     # Pytest suite with adversarial / mutation tests
```

---

## 4. Verification Procedures

### A. Standalone Lightweight Verifier (~0.5 seconds)
Execute the 6-gate verifier directly:
```bash
cd /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/shared_foundation
uv run python verify_shared_foundation.py
```

#### Gates Tested:
1. `[1/6]` **SHA256 Manifest Integrity**: Byte-exact checksum and file size verification of all 6 foundation assets.
2. `[2/6]` **Genome Manifest ID Count & Uniqueness**: Exactly 5,438 unique MAG IDs with 0 duplicates.
3. `[3/6]` **Canonical GTDB r232 Taxonomy Syntax & Rank Consistency**: Exact seven rank prefixes (`d__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`) and strict consistency between canonical rank columns and classification strings across all 5,438 genomes.
4. `[4/6]` **Archaea (591) / Thermoplasmatota (237) Subsets & Cross-Table Concordance**: Exact 591 Archaea, 237 Thermoplasmatota, and 354 other Archaea across the 7 verified GTDB phyla; 237 Thermoplasmatota subset matches `taxonomy_reconciliation_thermo237.tsv` 1:1 in ID set and full canonical GTDB r232 taxonomy string.
5. `[5/6]` **Metabolic Cohort Exact Join & Explicit External Cohort**: Exactly 484 unique models; 237 EQ desert IDs (`cohort == "EQ_desert"`) equal Thermoplasmatota 237; 247 external IDs (`cohort == "external_comparator"` explicitly verified); strict disjointness (0 overlap EQ vs ext, 0 overlap ext vs all 5,438 MAGs).
6. `[6/6]` **External Comparator Source Labels & Harmonization Boundary**: Verifies declared source labels (`GTDB_r232_rep` [218], `MGnify` [29]); documents that while source labels are verified as declared in inventory, cross-release genus/order harmonization between external comparator genomes and EQ MAGs is not independently audited by this join.

#### Informational Notice (Not a verification gate):
- **Scientific Interpretation Boundary & Pending Status**: Frozen inputs form the common cross-paper foundation, not proof that downstream analyses finished. Root owns scientific interpretation; ongoing metabolic modeling (Exp 1) and viral searches (Exp 3) remain pending independent review.

### B. Comprehensive Test Suite (Pytest)
Run the unit, regression, and adversarial mutation test suite:
```bash
cd /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/shared_foundation
uv run pytest tests/test_shared_foundation.py
```
The test suite proves that the verifier successfully catches:
- Substituted wrong IDs preserving the exact count of 5,438.
- Duplicate IDs in the genome manifest, taxonomy table, or metabolic inventory.
- Invalid Thermoplasmatota subsets (dropped member or substituted non-Thermo archaeon).
- Full taxonomy string or rank column mismatch between `taxonomy_reconciliation_thermo237.tsv` and `taxonomy_reconciliation_all5438.tsv.gz`.
- Invalid GTDB rank prefixes (e.g. `phylum__` instead of `p__`).
- Lineages missing ranks (e.g. 6 ranks instead of 7).
- Inconsistency between individual canonical rank columns and classification strings.
- Unexpected cohort values (must strictly be `'EQ_desert'` or `'external_comparator'`).
- Altered cohort membership (EQ vs external flipping, or ID collisions with Empty Quarter MAGs).
- Unexpected external source labels.
- Corrupted SHA256 checksums in manifest.

---

## 5. Explicit Limitations & Evidence Boundaries

1. **Taxonomy Harmonization Limitation**:
   - The 237 Empty Quarter Thermoplasmatota MAGs were reconciled to GTDB release 232 from earlier analytical taxonomy snapshots.
   - The 247 external comparator genomes carry declared inventory source labels from `GTDB_r232_rep` (n=218) and `MGnify` (n=29).
   - Full cross-release genus/order harmonization between external comparator genomes and Empty Quarter MAGs is not verified or harmonized by this join; documented uncertainty is preserved.

2. **Foundation vs Completed Downstream Analyses**:
   - These frozen tables represent the **input foundation** across manuscripts. They are **not proof that ongoing downstream mechanistic analyses have finished**.
   - Ongoing metabolic modeling (Experiment 1: 484-model FBA/FVA ensembles, nutrient deficiencies, and masking controls) and viral host-association analyses (Experiment 3: 591 archaeal genome CRISPR/prophage screen) remain under active review.
   - Confirmed biochemical functions are unassigned for dark exported protein families (Experiment 4); lead genomic neighborhood associations do not survive genus-dereplicated joint FDR.

3. **Interpretation Ownership**:
   - The scientific root owns interpretation.
   - No tentative viral or metabolic outcomes are presented as completed or integrated into companion manuscripts until independent review concludes.
