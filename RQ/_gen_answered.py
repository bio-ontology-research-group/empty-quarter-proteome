#!/usr/bin/env python3
"""Generate bounded methods/result snapshots for curated RQs and stubs for not-now RQs.

Statuses distinguish current verified claims from historical exploratory outputs; the presence of a generated
folder is not evidence that a result is promoted to the manuscript. Registered numbers are checked by
../reproducibility/verify.py where explicitly stated.
Run from repo root: python3 RQ/_gen_answered.py
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RQ = os.path.join(ROOT, "RQ")

# id: (title, status, answer_md, evidence_md, methods_md, refs)
ANS = {
"A1": ("How much of the EQ proteome is dark, and is it organised into families/folds?", "✅ answered",
 "**A historical 51.9 M-protein snapshot contains 2.84 M UniRef50-dark representatives and 334,810 recurrent "
 "sequence-dark families.** Its 68,808 pre-correction structural candidates are not fold calls. Under the current "
 "full-coordinate confidence/length gate, 608/1,875 are qTM dual no-hits and 97/1,875 are alignment-TM dual "
 "no-hits, against error floors of 101/3,013 and 68/3,013 in known-fold controls.",
 "`reproducibility/verify.py` (claims `b1.*`, PROVENANCE); source counts in the workbench darkmatter pipeline.",
 "Historical mmseqs2 clustering and UniRef50/AFDB screens are retained as provenance. The current headline uses "
 "full ESMFold coordinates, current PDB and AFDB v6 searches, both qTM and alignment-TM definitions, an identical "
 "length/confidence gate, and blinded known-fold controls. Heavy computation runs on IBEX/ORIX.",
 "pavlopoulos2023darkmatter; vanni2022agnostos; durairaj2023folds; yeo2025afesm; rodriguezdelrio2024fesnov"),
"A2": ("Is the dark proteome real (selected genes) or spurious ORFs?", "✅ answered",
 "**The recurrent, alignable subset carries coding-like evidence.** Its filtered pairwise NG86 median is "
 "**0.391**, and read recruitment finds only 1.3% in a single sample while 99.5% reach peak abundance ≥10. "
 "Pairwise within-population dN/dS is descriptive and is not a direct estimate of long-term selection; these "
 "results do not authenticate every predicted ORF.",
 "`results/darkmatter_ecology/dnds.tsv`; `results/darkmatter_ecology/ABUNDANCE.md`; reproduced by "
 "`code/darkmatter_abundance.py` and `reproducibility/verify.py` (`b2.*`, `b3d.*`).",
 "Per-family NG86 dN/dS over within-family pairs vs a size-matched control; CoverM read-recruitment per family "
 "across 551 samples for the abundance/prevalence statistics.",
 "nei1986dnds; kryazhimskiy2008dnds; rodriguezdelrio2024fesnov; aroney2025coverm"),
"A4": ("What concrete functions can be assigned to dark proteins, and what is the ceiling?", "✅ answered",
 "**10 well-supported specific-EC dark enzymes + 3 high-priority structural targets** were decoded (e.g. "
 "glucosylglycerol-phosphate phosphatase EC 3.1.3.69; uroporphyrinogen-III synthase EC 4.2.1.75). The decode is "
 "calibrated: P(correct) ≈ 0.34–0.48, so most of the dark proteome stays undecodable — an honest ceiling.",
 "`results/darkmatter_ecology/family_ec.tsv`; workbench `darkmatter_function_assignment` (C·V·T engine).",
 "Context (C) + Validation (V) + Templates (T) scoring with LRO calibration; FoldSeek structural evidence + "
 "genomic context + pathway gap-filling. (Provenance: workbench.)",
 "heinzinger2024prostt5; vankempen2024foldseek; durairaj2023folds"),
"B1": ("Is dark-proteome composition environmentally structured beyond space?", "✅ answered",
 "**Yes.** PERMANOVA shows compartment, climate|space and geochem|space each explain significant variance in "
 "dark-family composition *beyond* the spatial term; distance-decay Mantel r = +0.16 (p=0.001).",
 "`results/darkmatter_ecology/phase_b_summary.txt`; `code/phase_b_ecology.py`. (See also RQ-B7.)",
 "adonis/PERMANOVA per environmental block with conditioning on space; Mantel test community vs geography.",
 "fierer2017soil; maestre2015aridity"),
"B2": ("Is EQ dark matter endemic/desert-specific, or just globally undersampled?", "✅ answered",
 "**Database-relative specificity remains substantial.** Under one search rule, 36.6% of dark families have no "
 "detected MGnify90 homolog, 47.4% have no GEM homolog, 31.9% have no homolog in the union of five desert "
 "catalogues, and **21.4% have no detected homolog in the named MGnify90, GEM and five-desert comparison**. "
 "These no-hits prioritise broader searches but do not establish geographic absence or endemism.",
 "`results/global_specificity/`; `code/make_global_specificity_figure.py`; `reproducibility/verify.py` (`c3.*`).",
 "DIAMOND blastp `--sensitive` e<1e-3 of the 23,382 dark family reps vs each catalogue under one hit criterion; "
 "catalogue-specificity = % with no detected hit; the named-catalogue intersection is explicitly bounded.",
 "nayfach2021gem; richardson2023mgnify; uritskiy2019atacama"),
"B6": ("Are endemic genes more or less constrained than cosmopolitan ones?",
 "historical descriptive screen; not promoted as a constraint comparison",
 "In the archived, unclustered diagnostic, site-restricted families (≤3 of 60 sites) had median pairwise NG86 "
 "dN/dS 0.332, versus 0.420 for EQ-wide families (≥30 sites); Spearman(n_sites, dN/dS)=+0.18. The nominal "
 "Mann--Whitney p-value is retained in the source table, but orthogroup/pair dependence and unequal support are "
 "not modelled. The result therefore does **not** show that one class is under stronger evolutionary constraint. "
 "The occupancy labels are within-EQ descriptions, not evidence that either class is endemic or cosmopolitan.",
 "`results/darkmatter_ecology/ENDEMISM.md`; `code/phase_d_endemism.py`.",
 "Historical site occupancy per family (60 sites) joined to pairwise NG86 summaries; unclustered diagnostics only.",
 "nei1986dnds"),
"C1": ("Does 16S/PICRUSt2 predict the gene content of these classified taxa?", "✅ answered",
 "**Only partly.** Weighted NSTI is 0.16 and median KO presence-recall is 93%, yet 56.7% of MAG proteins lack a "
 "KO and signature-function abundance agreement is heterogeneous. The historical 2.40× richness ratio uses an "
 "asymmetric denominator; the shared relative-abundance rule gives 31.41×, so neither is retained as biological "
 "enrichment evidence.",
 "`results/picrust2_bias/`; `code/picrust2_bias_analysis.py`; `reproducibility/verify.py` (`c1.*`). Fig "
 "`figures/fig_picrust2_bias.pdf`.",
 "Historical PICRUSt2 v2.4.1 KO predictions from 16S (436 matched amplicon samples) vs MAG-encoded EggNOG KOs; "
 "the threshold sensitivity reports the asymmetric PIC-threshold/MAG-any-presence ratio separately from the "
 "shared relative-abundance threshold ratio. A current-catalogue raw-read remap replacement is dependency-gated.",
 "douglas2020picrust2; matchado2024limits; sun2020inference"),
"C3": ("Is predictability worse for more phylogenetically novel taxa?", "✅ answered",
 "**Yes.** The KO-space ceiling (fraction of proteins that any 16S predictor could see) rises with GTDB novelty — "
 "the more novel the lineage, the larger the unpredictable gene fraction.",
 "`figures/fig_picrust2_bias.pdf` panel c; `code/make_picrust2_bias_figure.py`.",
 "MAG no-KO fraction stratified by GTDB-Tk novelty (RED / lowest assigned rank).",
 "douglas2020picrust2; matchado2024limits"),
"D1": ("Is the community metabolically 'living on air' (atmospheric H2/CO)?", "✅ genomic-capacity census",
 "**Trace-gas oxidation capacity is widespread.** 3,386 MAGs carry high-affinity atmospheric [NiFe]-hydrogenase "
 "markers (groups 1l+1h) and 3,828 carry form-I CO dehydrogenase markers; **1,204 MAGs (~8%)** carry both. These "
 "genomic markers support metabolic potential, not a direct rate or activity measurement.",
 "`results/atmospheric/PHYLOGENY.md`; `reproducibility/verify.py` (`a1.*`). Fig `figures/fig_lifeonair_by_phylum.pdf`.",
 "HMM census (hydB/coxL) + HydDB typing; per-sample CoverM abundance of the guild.",
 "ji2017atmospheric; bay2021trace; greening2016genomic; sondergaard2016hyddb"),
"D2": ("Are the atmospheric enzymes bona-fide high-affinity forms, not wastebasket homologs?", "✅ answered",
 "**Bona-fide.** EQ hits classify to HydDB groups **1l (5,292 proteins) and 1h (986)** and CODH **form-I (5,645)**, "
 "and fall *within* the known atmospheric clades on the marker trees — not in wastebasket/heterotrophic groups.",
 "`results/atmospheric/{PHYLOGENY.md,hydB.nwk,coxL.nwk}`; `reproducibility/verify.py` (`a1.*`).",
 "Classification vs the HydDB reference set (Søndergaard 2016) + form-I CODH references; tree placement (FastTree).",
 "sondergaard2016hyddb; ji2017atmospheric"),
"D4": ("How is UV/desiccation survival encoded, and does it partition by sun exposure?", "✅ prevalence census",
 "**Repair-marker capacity is widespread and photolyase prevalence differs by compartment.** NER 99.4%, BER "
 "99.5% and NHEJ 78.4% of MAGs; photolyase-marker prevalence is 22.4% in surface versus 13.1% in deep samples. "
 "The observational contrast does not by itself establish sun exposure as the cause.",
 "`results/uv_survival/`; `code/uv_survival_analysis.py`; `reproducibility/verify.py` (`a2.*`). Fig amr/uv panels.",
 "KO/pathway presence across MAGs (group_prevalence) + compartment-stratified prevalence.",
 "makhalanyane2015microbial"),
"D5": ("Are compatible-solute (osmolyte) strategies convergent and salinity-linked?",
 "historical exploratory association; not a selection result",
 "Trehalose-, glucosylglycerol- and ectoine-linked annotations recur across multiple dark families, and the "
 "archived partial-correlation screen placed some families preferentially in saline/rhizosphere samples. "
 "Osmolyte-linked families were not enriched as a class (OR 1.10, nominal p=0.45). This exploratory association "
 "does not establish convergence, causation, or evolutionary selection; a current-catalogue, family-cluster-aware "
 "dose--response analysis would be required for those claims.",
 "`results/darkmatter_ecology/OSMOLYTE.md`; `code/phase_e_partial_osmolyte.py`. (Dose-response deepening: see RQ-B3.)",
 "Osmolyte-pathway family annotation + partial correlation of family abundance with the salinity proxy; class-level "
 "enrichment test.",
 "galinski1995osmoadaptation"),
"E1": ("Purifying vs diversifying selection across the dark proteome?",
 "partially answered; historical screens are descriptive",
 "The retained result is a filtered pairwise NG86 median of 0.391 for the recurrent, alignable dark-family subset. "
 "Pairwise within-population dN/dS is descriptive and is not a direct estimate of long-term selection. A historical "
 "FUBAR candidate-site screen is preserved under `positive_selection/`, but it did not provide a calibrated, "
 "family-aware basis for an interface or adaptive-selection conclusion. RQ-E2 records the bounded structural "
 "follow-up and its negative inferential boundary.",
 "`results/darkmatter_ecology/dnds.tsv` + `positive_selection/`; `code/phase_c_positive.py`.",
 "Historical per-family NG86 summaries and HyPhy-FUBAR candidate-site screening on aligned families; retained as "
 "exploratory provenance rather than a promoted selection test.",
 "nei1986dnds; murrell2013fubar"),
"F1": ("What is the secondary-metabolite (BGC) repertoire?", "✅ answered",
 "**Pervasive.** terpene in 64.5% of MAGs, NRPS 25.3%, ectoine 17.0% (+ T1PKS, RiPPs, etc.) — a broad natural-product "
 "reservoir. (Novelty vs MIBiG/BiG-FAM is RQ-F2.)",
 "`results/bgc_product_summary.tsv`; `code/make_figures.py`; `reproducibility/verify.py` (`a4.*`). Fig "
 "`figures/fig_bgc_product.pdf`.",
 "antiSMASH region calls aggregated by product type across MAGs.",
 "blin2023antismash; kautsar2021bigfam"),
"G1": ("Is the community competition- or cooperation-structured?", "✅ answered",
 "**Competition-leaning, function broadly distributed.** SMETANA global MRO ≈ 0.7 (metabolic-resource overlap = "
 "competition signal); cross-feeding is present but **not guild-anchored** — the atmospheric guild are not "
 "preferential metabolite donors (guild ≈ non-guild, p=0.69, robust to a completeness control).",
 "`results/smetana/CROSSFEED.md`; `code/smetana_analysis.py` (16-sample pilot). (Scaling = RQ-G2.)",
 "gapseq GEMs + SMETANA MRO/MIP + FBA cross-feeding on selected-sample HQ MAGs.",
 "zelezniak2015smetana; zimmermann2021gapseq"),
"H2": ("How does EQ compare to other hot/cold deserts?", "✅ answered relative to named comparison catalogues",
 "Under one historical DIAMOND rule, no-hit fractions against individual named desert catalogues ranged from "
 "42--89% (Atacama 89%, Antarctic 56%, Mojave 52%, arid-soil comparator 42%), and 31.9% of EQ dark families had "
 "no detected homolog in the union of the five named desert catalogues. These are database-relative detection "
 "results; they do not establish universal desert uniqueness or geographic absence.",
 "`results/global_specificity/MULTICATALOG.md`; `code/make_global_specificity_figure.py`.",
 "Same historical DIAMOND criterion against SPIRE per-study desert protein sets; per-catalogue and named-union "
 "no-hit fractions.",
 "uritskiy2019atacama; bay2021desert"),
"H3": ("Is there taxonomic novelty (new deep lineages)?", "✅ answered (bounded at high ranks; lower ranks open)",
 "**No current phylum- or class-gap candidate, but substantial lower-rank novelty.** GTDB-Tk 2.7.2/R232 assigns "
 "all 5,438 ANI95 HQ/MQ representatives through class across 51 named phyla. First assignment gaps occur for 11 "
 "representatives at order, 142 at family, 1,343 at genus and 3,755 at species. Ten order gaps entered de novo "
 "phylogenies and none advanced to proposal; one 46.49%-MSA candidate remains unadjudicated and in recovery.",
 "`reproducibility/verify.py` (`scale.*`); checksum-pinned R232 census and phylogeny artifacts under "
 "`reproducibility/`.",
 "GTDB-Tk 2.7.2/R232 first-empty-rank census over the canonical ANI95 HQ/MQ catalogue; separate LQ salvage "
 "sensitivity; de novo phylogeny/RED adjudication for staged order gaps.",
 "parks2018gtdb; chaumeil2022gtdbtk"),
}

# not-now RQs: brief status + plan (so every folder has content)
NOTNOW = {
"F2": ("How novel are the BGCs vs MIBiG/BiG-FAM?", "historical exploratory output; current analysis incomplete",
 "An archived BiG-SLICE/MIBiG run exists, but 90.1% of input antiSMASH regions were contig-edge fragments, so its "
 "99.3% no-reference-GCF result is an upper bound and is not promoted. A complete-region sensitivity analysis and "
 "BiG-FAM query are still required.", "kautsar2021bigfam; blin2023antismash"),
"E3": ("Is the dark proteome mobile (HGT / genomic islands / plasmid-borne)?", "🟡 answerable now (scope caveat)",
 "Operon-clustering was a clean negative (`results/atmospheric/integrate/OPERON_DARKMATTER.md`); dark proteins are "
 "NOT genomically isolated (`FUNCTIONAL_NOVELTY.md`). A dedicated genomic-island/MGE detection + dark-family-on-island "
 "enrichment is not yet run. **Coordinate with the sister phages/mobilome repo to avoid overlap.** Effort ~3–5 d.",
 ""),
"G2": ("Who feeds whom at biome scale (scale the cross-feeding)?", "🟡 answerable now (compute-heavy)",
 "Pilot = 16 samples (RQ-G1). Scaling gapseq+SMETANA to many more samples/MAGs is design-proven but needs weeks of "
 "IBEX compute. Deferred pending priority. ", "zelezniak2015smetana; zimmermann2021gapseq"),
"B4": ("Is the vanadium → alternative-nitrogenase link real?", "🟡 answerable now (leaning negative)",
 "Vanadium was *null* in the family-level PERMANOVA, and the V indicator table can be checked for nitrogenase "
 "families. A clean test needs per-sample vnf/anf marker abundance vs measured `V` across the 251 XRF sites "
 "(marker-abundance extraction from the workbench). Honest expected result: no V–vnf coupling. Effort ~2 d.", ""),
"I1": ("Can dark-gene function be predicted by integrating sequence + structure + context + ecology?", "🟡 answerable now",
 "The C·V·T decode already integrates these layers (RQ-A4). A held-out benchmark of the integrated predictor vs a "
 "sequence-only baseline would make it a methods contribution; not yet run. Effort ~3–5 d.", ""),
"I2": ("What is the minimal 'survival genome' for hyperarid oligotrophy?", "🟡 current analysis not run",
 "A current-catalogue core/accessory analysis joined to survival-function annotations is required. The archived "
 "RQ-C2 pangenome output is exploratory and is not a checksum-closed substitute.", ""),
"I3": ("Does the EQ proteome constitute a Mars/astrobiology analogue signature?", "🔴 framing / Discussion",
 "Not a data analysis — position the survival toolkit + dark folds as an arid-analogue reference in the Discussion. "
 "No new computation answers it.", ""),
"I4": ("Can EQ enrichment isolates validate dark-gene predictions in vitro?", "🔴 needs wet-lab",
 "Only 13 usable enrichment strains (the 185 isolates are aquatic — separate paper). Targeted expression / "
 "transcriptomics of dual-marker or decoded-enzyme strains is the natural wet-lab lever.", ""),
"F3": ("Can the dark proteome supply biotech-relevant (extremo)enzymes?", "✅ candidates / 🔴 validation",
 "**Candidates yes; validation no.** Decoded osmolyte phosphatases + 3 pocketed structural targets are expression "
 "candidates (RQ-A4/B5). Wet-lab heterologous expression/assay is new work (🔴).", ""),
"G3": ("Does the dark proteome participate in community metabolism (dark reactions)?", "🔴/🟡 partly addressable",
 "GEMs by construction only see decodable reactions. Partial route: map the decoded dark enzymes (RQ-A4) into GEM "
 "gap-fills and test whether they resolve network gaps; full answer needs broader decoding (🔴).", ""),
"H1": ("Is the database-relative dark signal contamination, kit-ome, or relic DNA?", "🟡 partial controls only",
 "A historical named-pathogen census and read-recruitment summaries are reassuring but do not exclude reagent "
 "contamination or relic DNA. A formal negative-control/decontam analysis and the staged relic-DNA workflow are "
 "required before a contamination-exclusion claim.", "salter2014contamination; davis2018decontam"),
}

def write(folder, title, status, answer, evidence, methods, refs):
    d = os.path.join(RQ, folder); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "result.md"), "w") as f:
        f.write(f"# RQ-{folder} — {title}\n\n**Status:** {status}\n\n## Answer\n{answer}\n")
        if evidence: f.write(f"\n**Evidence / reproduce:** {evidence}\n")
        if refs: f.write(f"\n**References:** {refs} (keys in `manuscript/references.bib`).\n")
    if methods is not None:
        with open(os.path.join(d, "methods.md"), "w") as f:
            f.write(f"# RQ-{folder} — methods\n\n{methods}\n")

n=0
for fid,(t,s,a,e,m,r) in ANS.items():
    write(fid,t,s,a,e,m,r); n+=1
skipped = 0
for fid,(t,s,a,r) in NOTNOW.items():
    # A formerly pending RQ may have been answered manually. Never overwrite that
    # richer evidence with an old planning stub; create stubs only when absent.
    if os.path.exists(os.path.join(RQ, fid, "result.md")):
        skipped += 1
        continue
    write(fid,t,s,a,None,"See `result.md` (status + plan).",r); n+=1
print(f"wrote {n} answered/new RQ folders; preserved {skipped} existing pending/manual folders")
