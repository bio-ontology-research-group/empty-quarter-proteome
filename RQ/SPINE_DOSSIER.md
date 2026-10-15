> ARCHIVED July 2026 planning document. Its biological interpretations and proposed novelty floors are superseded by manuscript/main.tex and the October round7 audits. Historical numbers below are retained for traceability, not current claims.

# SPINE DOSSIER — Empty Quarter MAGs Paper

**Purpose:** decisive, paper-facing consolidation of the 10 spine RQs. Every number is flagged **OLD** (pre-audit single-15,469 'cluster90' proteome) or **MERGED-VALID** (post-audit 2026-07-02 catalogue). The merged catalogue = single HQ+MQ **16,333** (was 15,469); self HQ+MQ **21,616**; merged dRep **9,253** species (composition: 5,715 recovered by both binnings, 2,886 self-supervised-only, 652 single-only — merging adds 3,538 either alone would miss; the 6,198/3,055 self/single figure is only which strategy won the *representative* slot, not who found the species).

**The pivot in one line:** the Empty Quarter's novelty is **functional and structural, not phylogenetic** — known organisms (0 novel phyla/classes) running unknown biology (vast endemic dark fold-space + community-wide life-on-air).

---

## 1. Summary table

| RQ | One-line answer | Catalogue status | Paper statement (headline number) |
|----|-----------------|------------------|-----------------------------------|
| **H3** Taxonomic novelty | None — the pivot. Conventional arid-soil tree. | **DONE-on-merged** | 0 novel phyla / 0 novel classes / 11 novel orders on merged 16,333 MAGs (~44 phyla). |
| **A1** Dark-proteome ladder | Vast unseen fold space organised at family scale. | **needs-recompute** (scales/holds) | 2.84M dark reps → 334,810 historical candidate novel-fold families → floor **10,007** with no AFDB/PDB analog. *OLD.* |
| **A2** Is dark matter real? | Yes — purifying selection + abundant recruitment. | **needs-recompute** (scales/holds) | dark median dN/dS **0.39** (99.3% <1) vs 0.29 control; median peak cov ~585×, only 1.3% single-sample. *OLD.* |
| **A3** Are novel folds genuine? | Raw 68,808 inflated ~7×; defensible floor ~10k. | **needs-recompute** (ratios stable) | 68,808 candidates → **10,007 (16.1%)** dual-DB TM<0.5 floor; mean pLDDT 0.51. *OLD.* |
| **C1** Does PICRUSt2 predict gene content? | No — good placement, wrong genome. | **scales-holds** (recompute optional) | ~**2.4×** overprediction (10,543 vs 4,482 KOs); 56.7% of proteins carry no KO. NSTI side catalogue-independent. *MAG side OLD.* |
| **C2** Open pangenomes? | Yes — mechanism behind C1; dark lives in accessory. | **ready-now** (rerunnable on merged dRep) | 97.7% open (median Heaps α **0.555**, ~6.5% core); dark **1.79×** enriched in accessory (p<1e-300). *pre-audit dRep.* |
| **D1** Living on air? | Yes, community-wide genomic capacity. | **needs-recompute** (fractions stable) | 22% H2 / 25% CO / **8% both**; 13.5% abundance-weighted, uniform across compartments. *OLD.* |
| **D2** Bona-fide high-affinity enzymes? | Yes — within known atmospheric clades, not wastebasket. | **needs-recompute** (verdict final; counts only) | 95% of hydB atmospheric (1l+1h); 5,645 form-I CODH (motif-concordant). *OLD.* |
| **H2** EQ-specific vs generic desert? | EQ-specific; endemism plateaus vs deep refs. | **needs-recompute** (scales/strengthens) | 31.9% absent from 5-desert union; 47.4% GEM; 36.6% MGnify90 — plateaus, not decays. *OLD.* |
| **B2** Endemic vs undersampled? | Genuinely endemic. | **needs-recompute** (scales/strengthens) | **21.4%** (5,011/23,382) no homolog in MGnify ∪ GEM ∪ deserts; plateaus at 32–47%. *OLD.* |

---

## 2. Per-RQ sections

### H3 — Taxonomic novelty: none (the pivot) — **MERGED-VALID**
- **Answer:** No novel deep lineages. GTDB-Tk r226 recovers a conventional arid-soil community.
- **Key numbers (MERGED, 16,333 single HQ+MQ):** 0 novel phyla, 0 novel classes, **11 novel orders**. OLD-catalogue detail not yet re-tallied on merged: ~44 phyla (37 bacterial + 7 archaeal), Actinomycetota 6,788 / Pseudomonadota 3,740 dominant. **Note:** CLAUDE.md "38 phyla" is stale; curated count is **44** (verify.py `scale.nphyla`).
- **Evidence:** `results/soil_gtdbtk/gtdbtk.{bac120,ar53}.summary.tsv`; `RQ/H3/result.md`+`methods.md`; verify.py `scale.*`. Merged job ID not in curated repo — pull from IBEX_SOIL_PROGRESS before quoting a job ID.
- **Reviewer risk:** GTDB-Tk under-calls deep novelty by construction (RED+ANI force-fit; MQ genomes may lack markers). "0 novel phyla" = no *classifiable* deep novelty. **Backstop to have in hand:** de novo concatenated-marker phylogeny of the 11 order-level outliers vs GTDB refs + RED/completeness-vs-placement check.
- **Phase-3 action:** **Done.** Optional light polish only — re-tally phylum count (44) and dominant-phylum counts on merged 16,333 for exact Table-1 numbers. Headline does not move.

### A1 — Dark-proteome novelty ladder — **OLD (needs-recompute)**
- **Answer:** The dark matter is unseen fold space, organised at family scale.
- **Key numbers (OLD):** 51.9M proteins → 23.4M mmseqs reps → **2.84M** UniRef50-dark reps → 2.53M structurally novel (89.2%) → **334,810** novel multi-member families → 68,808 candidate novel-fold families → validated floor **10,007** (16.1%; dual-DB TM<0.5).
- **Evidence:** `RQ/A1/result.md`+`methods.md`; verify.py `b1.*` (PROVENANCE), `retired legacy floor counter=10007` (COMPUTE via `RQ/A3/best_*_tm.tsv`, job 47844272). Ladder rungs are DOC/PROVENANCE — not re-derivable from committed tables.
- **Reviewer risk:** later rungs prediction-confidence-limited (mean pLDDT 0.51; only 9.4% ≥0.70) — hence retreat to 10,007 floor. AFDB near-saturation ("no analog" = undersampling?); MAG-called ORFs (rebutted by A2 selection + abundance). Every number on superseded catalogue.
- **Phase-3 action:** full end-to-end rerun **after merged proteome is built** — merged proteins → mmseqs reps → DIAMOND vs UniRef50 → ProstT5/ESMFold + FoldSeek vs AFDB+PDB → 30% family clustering → dual-DB TM floor. All six headline numbers provisional; conclusion expected to hold/strengthen.

### A2 — Is the dark proteome real? — **OLD (needs-recompute)**
- **Answer:** Real, not spurious — two independent lines.
- **Key numbers (OLD):** (1) **Selection** — dark median dN/dS **0.39** vs 0.29 size-matched control; 22,767 multi-member dark families (97.4%), 99.35% median <1; mean within-family id ~0.72 (dS not saturated). (2) **Abundance** — 23,382 families over 551 samples: only **1.3%** single-sample, 99.5% peak ≥10×, 89.4% ≥100×, median peak cov **584.8**, median family in 52/551 samples.
- **Evidence:** `results/darkmatter_ecology/dnds.tsv`, `darkfam_abundance.tsv`, `ABUNDANCE.md`; `code/darkmatter_abundance.py`; verify.py `b2.*`, `b3d.*`. Fully materialized — no new job for OLD result.
- **Reviewer risk:** dark dN/dS (0.39) is *higher* than control (0.29) — honest reading is "purifying but less constrained," not more. Skeptic: elevated ratio = young/fast-evolving or partly mis-called. Defenses: 99.3% <1, id ~0.72, abundance signature rules out contamination tail. NG86 crude; selection test covers only multi-member families (single-member reps rest on abundance alone).
- **Phase-3 action:** rerun **after merged proteome** — regenerate `dnds.tsv` (NG86 within-family + control) and CoverM `ABUNDANCE.md` over merged dark families. Distributional stats; qualitative conclusion expected to hold, point estimates move.

### A3 — Are the novel folds genuine? — **OLD (needs-recompute; ratios stable)**
- **Answer:** Raw "68,808 novel folds" inflated ~7×; replace with confidence-aware dual-DB floor.
- **Key numbers (OLD):** 68,808 candidates (68,527 ESMFold-modelled). pLDDT: mean **0.507**, only 9.4% ≥0.70, 52.9% ≥0.50. Real-coordinate reconfirm: PDB-novel 23,142 (33.6%); among high-confidence reps only 27.6% stay PDB-novel (~72% already have a PDB analog). Dual-DB census (job 47844272, 62,018 reps): AFDB TM<0.5 = 23.8%, PDB = 26.4%, **novel floor (both) = 10,007 (16.1%)**.
- **Evidence:** job 47844272 (`recover_afdb_convertalis.sh`) recovering timed-out job 47788341; `RQ/A3/best_{afdb,pdb}_tm.tsv`, `finalize_tm.py`, `fold_validation_summary.tsv`; verify.py `retired legacy floor counter` recomputes 10,007.
- **Reviewer risk:** the 10,007 floor is **not** itself confidence-filtered — most models low-confidence, TMalign unreliable at low pLDDT, so a "no-analog" call can be a poor-model artifact. Truly defensible HC core is smaller (only 1,651 HC reps PDB-novel; HC-AND-dual-DB-novel intersection not separately reported, likely low thousands). Per-family recurrence not attachable (novelfold rep-ids vs family-ids disjoint; only 46/68,808 match).
- **Phase-3 action:** rerun full A1 chain then dual-DB FoldSeek/TMalign census on merged proteome. Proportions (16.1%, 9.4% HC) near-invariant; only absolute ~10k floor scales up. **Recommended add:** report the HC-AND-dual-DB-novel intersection explicitly.

### C1 — Does PICRUSt2/16S predict gene content? — **MAG side OLD; NSTI side catalogue-independent**
- **Answer:** No — classification ≠ functional predictability. Good placement, wrong genome.
- **Key numbers:** 436 dual-data samples. Placement fine (NSTI median **0.159**, p95 0.322; housekeeping recall 93.3%). Failures: (1) **2.35×** overprediction (median 10,543 predicted vs 4,482 encoded KOs; Jaccard 0.384; ~6,363 phantom KOs/sample; marquee hoxH 99.1% pred/0.2% enc, anfG 99.1%/0.0%); (2) decoupled signature-function abundances (marquee median Spearman ~0.09; phr −0.123); (3) **KO-space ceiling** — 56.7% of proteins no KO, 19.1% no eggNOG, 3,376 encoded KOs never predicted.
- **Evidence:** `results/picrust2_bias/` (`per_sample_concordance.tsv`, `dark_fraction.tsv`, `marquee_functions.tsv`); `code/picrust2_bias_analysis.py`; verify.py `c1.*`; MAG-KO job 47804253; PICRUSt2 v2.4.1 job 45925397.
- **Reviewer risk:** overprediction magnitude/encoded denominator deflated by MAG incompleteness (MQ ~50% floor) + unbinned taxa. Rebuttals: NSTI good; marquee ~0% encoded at 99% predicted can't be manufactured by incompleteness; 56.7% KO-invisible is a fraction of *observed* proteins. Two different denominators — robust claims are direction + KO ceiling.
- **Phase-3 action:** **not in mandatory recompute set.** NSTI/amplicon side needs nothing. MAG side (cheap, mechanical): rerun `130_mag_ko_metagenome.py` + `picrust2_bias_analysis.py` once merged-proteome EggNOG annotations exist. Expected: ratio nudges below 2.35×, 56.7% shifts modestly, direction unchanged.

### C2 — Open pangenomes (the mechanism) — **pre-audit dRep; READY-NOW to rerun**
- **Answer:** Pangenomes wide open; dark matter sits in the open, unpredictable accessory — upgrades C1 to a first-principles mechanism.
- **Key numbers:** 345 non-redundant species (dRep ≥95% ANI, ≥10 genomes each; 8,635 genomes). **97.7% open** (Heaps α <1); median α **0.555**; median core/pan **0.065** (~6.5% core). Dark 8.24% of accessory vs 4.79% of core → **OR 1.786** (Fisher p<1e-300; 2×2 [8245/164589/20488/228945]).
- **Evidence:** `RQ/C2/140_pangenome_openness.{py,sbatch}`, `finalize.py`; job 47843640; `enrichment_summary.json`, `species_openness.tsv`; inputs `Cdb.csv` (col `secondary_cluster`, **May-13 pre-audit**, 37,085 members), Prodigal proteomes, `is_dark.tsv`. Params MIN_GENOMES=10 CORE_FRAC=0.95 …; mmseqs easy-cluster --min-seq-id 0.5 -c 0.8.
- **Reviewer risk:** MAG incompleteness confound — MQ genomes miss genes, mechanically shrinking core and inflating accessory → α biased "open." Needs completeness-stratified / HQ-only sensitivity check (merged adds more MQ). ≥10-genome threshold restricts to well-sampled species (345/600 raw); "dark" is DB-relative; enrichment on rep proteomes only.
- **Phase-3 action:** **Runnable now** on the post-audit merged dRep Cdb.csv (9,253) once python+mmseqs2 env reconstructed. **Two halves:** (a) **openness** (Heaps α, core/pan) runs immediately from Cdb.csv + Prodigal proteomes — can only reinforce openness; (b) **dark-enrichment OR** depends on `is_dark.tsv`, which needs the merged-proteome rebuild first. **Blocker note:** the *script default* Cdb.csv is the May-13 pre-audit one — point it at our post-audit merged re-dRep (`04_drep_merged/out/data_tables/Cdb.csv`, 9,253).

### D1 — Living on air (capacity) — **OLD (needs-recompute; fractions stable)**
- **Answer:** Yes, community-wide genomic capacity for atmospheric H2/CO oxidation.
- **Key numbers (OLD):** **3,386 MAGs (22%)** high-affinity [NiFe]-hydrogenase (HydDB 1l+1h; 95% of hydB high-affinity; desert Group 1l dominant 5,292 prot/2,609 MAGs vs Huc 1h 986/893). **3,828 (25%)** form-I CODH (5,645 proteins; 8,003 ambiguous excluded). Dual-marker guild core **1,204 (~8%)**. Abundance-weighted (D3, 1,227 amplicon samples): median **13.5%** (IQR 9.6–18.0), uniform across Surface/Deep/PlantRoot. gapsmith FBA: 616/648 (95%) guild MAGs carry a hydrogenase reaction.
- **Evidence:** jobs 47554112/47554035/47554113; `results/atmospheric/` (`hydB_class*.tsv`, `coxL_*`, `guild_core_mags.txt`, `atmo_guild_by_phylum.tsv`, trees); verify.py `a1.*`; D3 `results/atmospheric/abundance/`. Culture cross-check: Group 1l absent from 1,221 cultures (dominant scavenger uncultured).
- **Reviewer risk:** genomic **capacity, not flux** — no transcriptomics/rates ("living on air" inferred from markers; mitigated by D2 clades + FBA + abundance, true flux stays open). 8,003 ambiguous coxL excluded (conservative but boundary is a judgement call). Counts on pre-audit catalogue.
- **Phase-3 action:** DIAMOND census is light; **blocker is the merged proteome build**, not compute. Rerun hydB/coxL DIAMOND+HydDB+motif on merged proteins; re-derive per-MAG counts + guild core; re-map D3 16S weighting via GTDB. Absolute counts shift; fractions (~8%/22%/25%/13.5%) and conclusion expected stable.

### D2 — Bona-fide high-affinity enzymes? — **OLD counts; verdict FINAL**
- **Answer:** Bona-fide, not wastebasket — markers classify to high-affinity forms and fall within known atmospheric clades.
- **Key numbers (OLD):** hydB 6,278/6,602 (**95%**) atmospheric (1l 5,292/2,609 MAGs; 1h 986/893); minor tail across 10 groups (2a 246, 1d 32…; median id ~58%) — not a one-reference artifact. coxL: 5,645 proteins/3,828 MAGs meet form-I (>55% id, ≥90% qcov + AYRGAGR motif); 8,003 ambiguous excluded; only 27 form-II/other. Prevalence 22%/25%/8%. FBA corroboration 616/648.
- **Evidence:** `results/atmospheric/PHYLOGENY.md`, `hydB.nwk`, `coxL.nwk` (FastTree, EQ reps + refs); jobs 47554112/47554035/47554113; verify.py `a1.hydB1l.prot=5292` etc. Refs sondergaard2016hyddb, ji2017atmospheric, cordero2019.
- **Reviewer risk:** capacity not activity. Form calls rest on DIAMOND best-hit + identity/motif threshold; 8,003 boundary proteins "ambiguous." Trees are the anti-wastebasket armor but the **annotated publication tree figure is not yet rendered** — reviewer currently sees best-hit tables.
- **Phase-3 action:** **not in mandatory dark-matter recompute set** — bona-fide verdict is a per-protein clade/form property independent of catalogue size. Merged rerun is **count-level only** (refresh protein/MAG counts + prevalence %, whose 15,469 denominator is shared with D1). **Recommended:** render the annotated hydB/coxL tree figure.

### H2 — EQ-specific vs generic desert dark matter? — **OLD (needs-recompute; strengthens)**
- **Answer:** Substantially EQ-specific, not shared generic desert dark matter.
- **Key numbers (OLD):** same 23,382 UniRef50-dark reps, DIAMOND blastp --sensitive e<1e-3, no-hit = endemic. Absent per desert: arid 41.7%, Mojave 51.9%, Antarctic 55.9%, autotroph 62.8%, Atacama 89.1%. Absent from **5-desert union (147M) 31.9%** (7,469); **GEM (123M) 47.4%** (11,075); **MGnify90 (718M) 36.6%** (8,552). Plateaus at ~32–47% at deepest refs. Cross-catalogue floor (=B2) **21.4%** (5,011).
- **Evidence:** `results/global_specificity/MULTICATALOG.md`, `desert_specificity_summary.tsv`, `global_specificity_summary.tsv`; query `endemism/query.faa`; scripts `131_global_specificity.sh` (job 47810126, SPIRE still finishing), `131b_desert_specificity.sh` (job 47810171); verify.py `c3.*`; `paper/main.tex` L216/L318.
- **Reviewer risk:** "absent" bounded by comparator completeness + DIAMOND sensitivity (misses remote homologs, overstates endemism). Comparators uneven — shallow deserts (Atacama 6M) mechanically inflate absence; cross-desert numbers partly track depth. Defense = plateau + deep anchors. SPIRE not yet folded in. All on OLD proteome.
- **Phase-3 action:** re-derive dark query set from merged proteome (query.faa size/composition changes), then DIAMOND vs the **fixed external comparators** (MGnify90/GEM/5-deserts unchanged, MGnify row cached) via scripts 131/131b. Only the query side moves; conclusion expected to survive/strengthen. **Fold in SPIRE (job 47810126).**

### B2 — Genuinely endemic vs undersampled? — **OLD (needs-recompute; strengthens)**
- **Answer:** Genuinely endemic — endemism plateaus rather than decaying with reference depth.
- **Key numbers (OLD):** same 23,382 reps. Absent MGnify90 36.6%, GEM 47.4%, 5-desert union 31.9%. **Intersection: 5,011/23,382 = 21.4%** with no homolog in MGnify ∪ GEM ∪ deserts (with-hit anywhere 18,371). Plateau at 32–47% = endemism signature, not depth artifact. Corroboration (B6): endemic families under stronger purifying selection (≤3-site median dN/dS 0.332 vs EQ-wide 0.420, MWU p<1e-67).
- **Evidence:** scripts 131/131b (jobs 47810126/47810171); `global_specificity_summary.tsv`, `desert_specificity_summary.tsv`, `MULTICATALOG.md`; MGnify count in `darkfam_meta.tsv`/GLOBAL_ENDEMISM; verify.py `c3.intersection` (21.4, 5011/23382); B6 dN/dS in `ENDEMISM.md`.
- **Reviewer risk:** "no homolog anywhere" is only an **upper bound** on novelty (DIAMOND misses remote homologs; DBs keep growing) — reviewer: 21.4% shrinks toward 0 with deeper refs. Defense = plateau. SPIRE unfinished → quote as **"≤21.4% pending SPIRE."** Spurious-ORF concern countered by H1 controls + purifying selection.
- **Phase-3 action:** same as H2 — rebuild dark query from merged proteome → DIAMOND vs unchanged external refs → recompute 4 headline fractions + intersection. Adding self-sup EQ genomes should add EQ-specific families → conclusion holds/strengthens. SPIRE only tightens (lowers) 21.4%.

---

## 3. Phase-3 recompute checklist (ordered by dependency)

### Tier 0 — Ready now (no merged-proteome build required)
1. **H3 — essentially DONE.** 0 novel phyla/classes / 11 orders already merged-valid on 16,333. *Optional:* re-tally 44-phylum count + Actinomycetota/Pseudomonadota dominance on merged 16,333 for exact Table-1 numbers. Correct the stale "38 phyla" in CLAUDE.md → 44. **Reviewer backstop:** de novo marker phylogeny of the 11 order-level outliers.
2. **C2 — openness half, runnable now.** Reconstruct python + mmseqs2 env, point at the **post-audit merged re-dRep** (`04_drep_merged/out/data_tables/Cdb.csv`, 9,253 reps), rerun `140_pangenome_openness.{py,sbatch}` (32cpu/256G/24h) for Heaps α + core/pan. Can only reinforce openness. *(Dark-enrichment OR half is gated on Tier 1.)*

### Tier 1 — Critical blocker: **build the merged proteome catalogue**
> **Step 0 (gates A1, A2, A3, H2, B2, C1-MAG, C2-dark-half, D1, D2-counts):** assemble/predict proteins over the merged MAG set (9,253-rep merged dRep / 16,333 merged single + 21,616 self) → mmseqs reps → DIAMOND vs UniRef50 → **new dark-family reps** (replaces the 23,382 query set) + merged-proteome EggNOG annotations + rebuilt `is_dark.tsv`.

Once Step 0 lands, in dependency order:

3. **A1** — full novelty ladder: reps → UniRef50-dark → ProstT5/ESMFold + FoldSeek vs AFDB+PDB → 30% family clustering → dual-DB TM floor. *(Feeds A3.)*
4. **A3** — dual-DB FoldSeek/TMalign census (job-47844272 style) on A1's ESMFold models. Also add the HC-AND-dual-DB-novel intersection.
5. **A2** — regenerate `dnds.tsv` (NG86 within-family + size-matched control) and CoverM `ABUNDANCE.md` over the new dark families.
6. **C2 (dark-enrichment half)** — recompute the accessory-vs-core dark OR (was 1.786) using the rebuilt `is_dark.tsv`. *(Openness half already done in Tier 0.)*
7. **H2 / B2** — re-derive dark query.faa, DIAMOND blastp --sensitive e<1e-3 vs the **unchanged** external refs (MGnify90/GEM/5-deserts, MGnify cached) via 131/131b; recompute fractions + 21.4% intersection. **Fold in SPIRE (job 47810126).**
8. **C1 (MAG side)** — rerun `130_mag_ko_metagenome.py` + `picrust2_bias_analysis.py` on merged-proteome EggNOG annotations. *(PICRUSt2/NSTI amplicon side needs nothing.)*
9. **D1 / D2 (counts)** — rerun hydB/coxL DIAMOND+HydDB+motif census on merged proteins; re-derive guild core + prevalence %; re-map D3 16S weighting via GTDB. **D2 verdict is final** — counts only.

### Independent loose end
- Render the annotated **hydB/coxL publication tree figure** (D2 armor) — not catalogue-dependent, do anytime.

---

**Bottom line for the manuscript:** H3 is the one spine result already clean on the merged catalogue and anchors the pivot. C2's openness signal is one env-reconstruction away from merged-valid. Everything else (A1/A2/A3/H2/B2 mandatory; C1-MAG, D1, D2-counts derived) is gated on building the merged proteome catalogue — **no dark-matter/protein number should be quoted for the merged catalogue until Step 0 lands.** In all cases the *direction* of the conclusion is expected to hold or strengthen; only point estimates move.

_Generated 2026-07-02 (spine-dossier workflow, 10 per-RQ distillers + synthesis)._
