# Exp 4 Surface Genomic Module Synteny Figures

This self-contained reproducibility module generates publication-quality main and supplementary figures illustrating the syntenic gene architecture of two candidate surface-associated dark protein families identified in the genome-resolved Rub' al-Khali archaeal unknown-proteome study (Experiment 4).

---

## 1. Reproduction Commands

### Generate Figures
```bash
# Build main figure (fig_surface_modules.pdf/png) and supplementary atlas (fig_surface_modules_all_loci.pdf/png)
uv run python3 build_surface_figure.py
```

### Run Test Suite
```bash
# Run all 11 unit and regression tests
uv run pytest -v test_surface_figure.py
```

---

## 2. Input Provenance & Compact Data

1. **`exp4_lead_flank_coordinates_independent.tsv`** (SHA256: `75c51a9009ae46dbde21dc96a4932cf3d540772b67066e45f31b91a8c4923866`):
   - 24 verified independent locus pairs for the two surface-associated leads:
     * **Lead A** (`22PRr1::22PRr1_2531955_3`): 13 recovered pairs across 6 genus-or-singleton blocks (8 in genus `g__JAJPHT01`, 5 unclassified singletons); gaps range from -4 to 37 bp.
     * **Lead B** (`V60PRr2::V60PRr2_1314979_4`): 11 recovered pairs across 5 genus-or-singleton blocks (7 in genus `g__JAJPHT01`, 4 unclassified singletons); gaps range from -4 to 6 bp.
   - Source provenance: identical bit-for-bit to audited evidence in `reproducibility/mechanism_20261109/exp4_surface/compact_evidence/exp4_lead_flank_coordinates_independent.tsv`.
   - Core structural invariants: All 24 loci are strictly co-oriented on the same contig strand (`same_strand == True`), with the candidate gene positioned 5' (upstream) to the partner in the direction of transcription. Negative gaps denote translational/genomic overlaps under 1-based inclusive coordinates (`flank_start - focal_end - 1` for strand +1, `focal_start - flank_end - 1` for strand -1).

2. **`exp4_surface_final_summary_v2_audited.json`** (SHA256: `18327f271eb58d53f1410a2cc9a96113b108a1fc1511c8dcb0e0b446d9b4430c`):
   - Audited summary from the 47-family deep dive across 428 loci (326 verified complete interior), recording the final dispositions:
     * 38 uncertain / insufficient evidence (80.9%)
     * 9 exploratory fold or motif candidate without operon convergence (19.1%)
     * 0 biochemical function assignments validated (0.0%)
   - Full 1,460-family / 8,760-test context screen joint FDR: minimum genus-restricted joint $q = 1.0$ (no significant genus-level associations). Primary genome-level joint FDR: $q = 0.00206$ (Lead A) and $q = 0.00035$ (Lead B).
   - Source provenance: identical bit-for-bit to audited evidence in `reproducibility/mechanism_20261109/exp4_surface/compact_evidence/exp4_surface_final_summary_v2_audited.json`.

3. **`locus_id_mapping.tsv`** (SHA256: `e54afec45859a1f20b61b6c8994d5361071f00da235dc124f3e512d2a2717109`):
   - Complete crosswalk linking short unique IDs `L01`–`L24` to sample IDs, genome names, genus/singleton blocks, coordinates, strands, and exemplar representative flags.

---

## 3. Output Figures

- **`fig_surface_modules.pdf`** & **`fig_surface_modules.png`** (Main Figure):
  - Compact two-panel layout ($11.0\text{ in} \times 5.5\text{ in}$, 300 DPI) displaying verified complete interior representatives:
    * **Panel a (Lead A: PilJ-profile neighbour)**: exemplar `1PRr2__self_SemiBin_81` (`L01`, strand -1), Candidate (1,245 bp) and Partner (471 bp) separated by a +5 bp gap.
    * **Panel b (Lead B: envelope-profile neighbour)**: exemplar `V47PRr1__self_SemiBin_699_sub` (`L14`, strand +1), Candidate (684 bp) and Partner (339 bp) separated by a +5 bp gap.
  - Common horizontal axis: `Position relative to pair left edge (bp)` ($0$ to $1{,}800\text{ bp}$).
  - Single common top legend for Candidate (teal `#0E7490`) and Partner (orange `#D97706`).

- **`fig_surface_modules_all_loci.pdf`** & **`fig_surface_modules_all_loci.png`** (Supplementary Figure):
  - Comprehensive synteny atlas displaying all 24 verified loci with short IDs `L01`–`L24`, sample IDs, genus classifications, strand orientations, and exact intergenic gaps (including -4 bp and -1 bp overlaps).

---

## 4. Figure Captions & Audited Disclosures

### Main Figure: Conserved syntenic gene-neighbourhood modules of two candidate dark protein families across desert archaeal genomes.

**a, Lead A: PilJ-profile neighbour.** Syntenic architecture of the exemplar verified interior locus 1PRr2__self_SemiBin_81 (L01; sample 1PRr2, genus g__JAJPHT01). The candidate dark family 22PRr1::22PRr1_2531955_3 (1,245 bp, teal arrow) is immediately adjacent to a 'PilJ-profile partner' (1PRr2_468829_2; 471 bp, orange arrow; bacterial PilJ PF04964 profile match) separated by an intergenic gap of +5 bp. Both genes reside on the reverse contig strand (strand -1; native orientation pointing left, transcription candidate to partner). This locus is the verified interior representative for n = 13 verified context pairs recovered across six genus-or-singleton blocks (eight loci in shared genus g__JAJPHT01 and five independent unclassified singletons), with intergenic gaps ranging from -4 to 37 bp (including negative translational overlaps of -4 bp and -1 bp; see Supplementary Figure S_ for all 13 loci).

**b, Lead B: envelope-profile neighbour.** Syntenic architecture of the exemplar verified interior locus V47PRr1__self_SemiBin_699_sub (L14; sample V47PRr1, genus g__JAJPHT01). The candidate dark family V60PRr2::V60PRr2_1314979_4 (684 bp, teal arrow) is immediately adjacent to an 'Envelope-profile partner' (V47PRr1_105096_19; 339 bp, orange arrow; Pfam PF18204/PF10731 profile match) separated by an intergenic gap of +5 bp. Both genes reside on the forward contig strand (strand +1; native orientation pointing right, transcription candidate to partner). This locus is the verified interior representative for n = 11 verified context pairs recovered across five genus-or-singleton blocks (seven loci in shared genus g__JAJPHT01 and four independent unclassified singletons), with intergenic gaps ranging from -4 to 6 bp (including negative overlaps of -4 bp and -1 bp; see Supplementary Figure S_ for all 11 loci).

**Orientation and coordinate representation:**
Coordinates (bp) are plotted on a common scale relative to the pair left edge on the contig (x-axis: 'Position relative to pair left edge (bp)'; pair left edge = 0 bp). Arrow lengths correspond strictly to verified coding sequence lengths (actual start and end coordinates; no imputed or invented flanks). Arrow directions indicate native contig strand (+1: forward / right; -1: reverse / left). Across all 24 loci, candidate and partner genes are strictly co-oriented on the same contig strand, and the candidate is consistently positioned 5' (upstream) to the partner in the direction of transcription.

**Evidence boundaries and scientific interpretation:**
Candidate biochemical and cellular functions remain experimentally unresolved. Partner annotations ('PilJ-profile partner' and 'Envelope-profile partner') reflect sequence and HMM profile matches only and do not imply validated partner function or assign candidate function. Structural fold matching against curated reference panels yielded no qualifying fold matches for either candidate (max qTM = 0.2898 to 8RH5 ArlB for Lead A; no alignment to positive S-layer 7PTR for Lead B), and verified models exist only for contig-edge seed members with exact sequence matches rather than complete interior representatives.
Statistical significance reflects primary genome-level context associations that survive joint Benjamini-Hochberg FDR correction across the 1,460-family screen (Lead A joint FDR q = 0.00206; Lead B joint FDR q = 0.00035). However, genus-restricted joint FDR is 1.0 throughout (min genus joint q = 1.0 across all tested families), indicating that local synteny supports conserved module hypotheses across multiple genus-or-singleton blocks but provides genus-limited statistical support and does not establish independent genus-level replication.

---

### Supplementary Figure S_: Complete synteny alignment of all 24 verified loci across desert archaeal genomes for Lead A and Lead B.
Panels a and b display all 13 and 11 recovered loci, respectively, identified by short unique identifiers L01–L24 alongside sample names and genus-or-singleton block classifications (full mapping provided in locus_id_mapping.tsv). Verified interior representatives (L01 and L14) are highlighted. Arrow directions encode contig strand (+1 right, -1 left) and lengths reflect exact verified coordinates relative to the pair left edge. Intergenic gaps (ranging from -4 to 37 bp in Lead A and -4 to 6 bp in Lead B, including -4 bp and -1 bp translational overlaps) are annotated per locus.
