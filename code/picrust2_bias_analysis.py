#!/usr/bin/env python3
"""C1 reference-bias analysis: PICRUSt2 16S-predicted KO content vs MAG-encoded reality.

Even though GTDB classifies the Empty Quarter taxa (known phyla/genera), is their *gene content*
predictable from 16S? PICRUSt2 (the standard tool for amplicon-only studies) places ASVs on a
reference tree and copies the reference genomes' KO content. We test how well that recovers what the
genome-resolved MAGs from the SAME samples actually encode.

Inputs
  MAG side  (pulled from IBEX 05_picrust_validation/ -> results/picrust2_bias/mag/):
    mag_ko_relab_by_sample.tsv     KO x sample, CoverM-relab-weighted copy number (PICRUSt2-comparable)
    mag_ko_count_by_sample.tsv     KO x sample, unweighted pangenome copy number
    mag_ko_presence_by_sample.tsv  KO x sample, # MAGs encoding KO
    mag_annotation_stats.tsv       per-MAG: sample,bin,quality,n_proteins_eggnog,n_proteins_ko,coverm_relab,eggnog_present
  PICRUSt2 (amplicon repo, already computed):
    merged/ko_metagenome_unstrat.tsv   KO x <picrust2_col> predicted abundance
    merged/weighted_nsti.tsv           per-sample weighted NSTI
    merged/sample_metadata.tsv         picrust2_col -> sample_id (+ compartment etc.)
  Local: workbench/results/soil_gtdbtk/checkm2_all_samples.tsv (Total_Coding_Sequences -> dark denom)
         workbench/results/soil_gtdbtk/gtdbtk.bac120.summary.tsv (phylum, optional)

Out: results/picrust2_bias/{per_sample_concordance.tsv, marquee_functions.tsv, dark_fraction.tsv, SUMMARY.md}
Additional correction input: results/picrust2_bias/rebuild/threshold_sweep.tsv
"""
import os, csv, math, statistics as st
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
OUT = os.path.join(REPO, "results", "picrust2_bias")
MAG = os.path.join(OUT, "mag")
AMP = os.path.expanduser("~/Public/software/empty-quarter/data/processed/functional/picrust2/merged")
CK = os.path.join(REPO, "workbench", "results", "soil_gtdbtk", "checkm2_all_samples.tsv")
GT = os.path.join(REPO, "workbench", "results", "soil_gtdbtk", "gtdbtk.bac120.summary.tsv")
SWEEP = os.path.join(OUT, "rebuild", "threshold_sweep.tsv")
REBUILT = os.path.join(OUT, "rebuild", "per_sample_rebuilt.tsv")

MARQUEE = {  # KO -> (gene, function)  the EQ "story" functions
    "K03520": ("coxL", "CO dehydrogenase form-I large subunit (atmospheric CO)"),
    "K03519": ("coxM", "CO dehydrogenase medium subunit"),
    "K03518": ("coxS", "CO dehydrogenase small subunit"),
    "K06281": ("hyaB", "[NiFe]-hydrogenase large subunit (H2 oxidation)"),
    "K00437": ("hoxH", "[NiFe]-hydrogenase"),
    "K06720": ("ectC", "ectoine synthase (osmolyte)"),
    "K00836": ("ectB", "ectoine biosynthesis aminotransferase"),
    "K06718": ("ectA", "ectoine biosynthesis acetyltransferase"),
    "K00697": ("otsA", "trehalose-6-P synthase (osmolyte)"),
    "K01087": ("otsB", "trehalose-6-P phosphatase"),
    "K01669": ("phr",  "deoxyribodipyrimidine photolyase (UV repair)"),
    "K00531": ("anfG", "nitrogenase (alternative/Fe-only)"),
}

def spearman(x, y):
    pts = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    n = len(pts)
    if n < 4:
        return float("nan"), n
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(v):
            j = i
            while j+1 < len(v) and v[order[j+1]] == v[order[i]]:
                j += 1
            for k in range(i, j+1):
                r[order[k]] = (i+j)/2.0+1
            i = j+1
        return r
    rx, ry = rank([p[0] for p in pts]), rank([p[1] for p in pts])
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = (sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry)) ** 0.5
    return (num/den if den else float("nan")), n

def read_matrix(path):
    """KO x sample -> (cols, {KO: {sample: val}})"""
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        cols = hdr[1:]
        m = {}
        for line in f:
            c = line.rstrip("\n").split("\t")
            m[c[0]] = {cols[i]: float(c[i+1]) for i in range(len(cols)) if c[i+1] != "0"}
    return cols, m

def main():
    # ---- sample_id <-> picrust2_col ----
    col2sid, sid2comp = {}, {}
    with open(os.path.join(AMP, "sample_metadata.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            col2sid[r["picrust2_col"]] = r["sample_id"]
            sid2comp[r["sample_id"]] = r.get("compartment", "")
    # ---- PICRUSt2 KO metagenome (KO x picrust2_col) -> reindex to sample_id ----
    pic_cols, pic = read_matrix(os.path.join(AMP, "ko_metagenome_unstrat.tsv"))
    pic_samples = {col2sid.get(c, c) for c in pic_cols}
    pic_by_sid = defaultdict(dict)  # KO -> sid -> val (collapse dup runs by max)
    for ko, d in pic.items():
        for c, v in d.items():
            sid = col2sid.get(c, c)
            pic_by_sid[ko][sid] = max(pic_by_sid[ko].get(sid, 0), v)
    # ---- NSTI ----
    nsti = {}
    with open(os.path.join(AMP, "weighted_nsti.tsv")) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            nsti[col2sid.get(c[0], c[0])] = float(c[1])
    # ---- MAG matrices ----
    _, mag_relab = read_matrix(os.path.join(MAG, "mag_ko_relab_by_sample.tsv"))
    mag_cols, mag_count = read_matrix(os.path.join(MAG, "mag_ko_count_by_sample.tsv"))
    mag_samples = set(mag_cols)

    matched = sorted(mag_samples & pic_samples)
    # invert MAG matrices to per-sample dicts
    def per_sample(mat):
        out = defaultdict(dict)
        for ko, d in mat.items():
            for s, v in d.items():
                out[s][ko] = v
        return out
    pic_s = defaultdict(dict)
    for ko, d in pic_by_sid.items():
        for s, v in d.items():
            pic_s[s][ko] = v
    magrelab_s = per_sample(mag_relab)
    magcount_s = per_sample(mag_count)

    # KOs the MAGs encode that PICRUSt2 NEVER predicts in ANY sample
    mag_all_kos = set(mag_count)
    pic_all_kos = set(pic_by_sid)
    mag_only_global = mag_all_kos - pic_all_kos

    # ---- per-sample concordance ----
    rows = []
    for s in matched:
        Mset = set(magcount_s[s]); Pset = set(pic_s[s])
        if not Mset:
            continue
        inter = Mset & Pset
        recall = len(inter)/len(Mset)
        union = Mset | Pset
        jacc = len(inter)/len(union) if union else float("nan")
        # abundance Spearman over union (zeros filled), MAG relab vs PICRUSt2
        xs = [pic_s[s].get(k, 0.0) for k in union]
        ys = [magrelab_s[s].get(k, 0.0) for k in union]
        rho, n = spearman(xs, ys)
        rows.append(dict(sample=s, compartment=sid2comp.get(s, ""), nsti=nsti.get(s, float("nan")),
                         n_mag_ko=len(Mset), n_pic_ko=len(Pset), n_shared=len(inter),
                         recall=recall, mag_only=len(Mset-Pset), pic_only=len(Pset-Mset),
                         jaccard=jacc, spearman=rho))
    with open(os.path.join(OUT, "per_sample_concordance.tsv"), "w") as o:
        flds = ["sample","compartment","nsti","n_mag_ko","n_pic_ko","n_shared","recall",
                "mag_only","pic_only","jaccard","spearman"]
        o.write("\t".join(flds)+"\n")
        for r in rows:
            o.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in flds)+"\n")

    med = lambda key: st.median([r[key] for r in rows if r[key] == r[key]])

    # ---- dark fraction (PICRUSt2/KO-space invisibility ceiling) ----
    tcds = {}
    with open(CK) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                tcds[f"{r['sample']}__{r['Name']}"] = int(r["Total_Coding_Sequences"])
            except (ValueError, KeyError):
                pass
    dark_rows = []
    sum_cds = sum_ko = sum_eggnog = 0
    with open(os.path.join(MAG, "mag_annotation_stats.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["eggnog_present"] != "1":
                continue
            mid = f"{r['sample']}__{r['bin']}"
            cds = tcds.get(mid)
            nko = int(r["n_proteins_ko"]); negg = int(r["n_proteins_eggnog"])
            if not cds or cds <= 0:
                continue
            frac_no_ko = 1 - nko/cds
            frac_no_eggnog = 1 - negg/cds
            dark_rows.append((mid, r["sample"], r["quality"], cds, nko, negg, frac_no_ko, frac_no_eggnog))
            sum_cds += cds; sum_ko += nko; sum_eggnog += negg
    with open(os.path.join(OUT, "dark_fraction.tsv"), "w") as o:
        o.write("mag\tsample\tquality\ttotal_cds\tn_ko\tn_eggnog\tfrac_no_ko\tfrac_no_eggnog\n")
        for row in dark_rows:
            o.write("\t".join(str(x) for x in row[:6]) + f"\t{row[6]:.4f}\t{row[7]:.4f}\n")
    prot_frac_no_ko = 1 - sum_ko/sum_cds if sum_cds else float("nan")
    prot_frac_no_egg = 1 - sum_eggnog/sum_cds if sum_cds else float("nan")
    med_mag_no_ko = st.median([row[6] for row in dark_rows]) if dark_rows else float("nan")

    # ---- marquee functions: PICRUSt2 vs MAG prevalence/abundance across matched samples ----
    marq = []
    for ko, (gene, desc) in MARQUEE.items():
        pic_prev = sum(1 for s in matched if pic_s[s].get(ko, 0) > 0)/len(matched)
        mag_prev = sum(1 for s in matched if magcount_s[s].get(ko, 0) > 0)/len(matched)
        pic_ab = [pic_s[s].get(ko, 0) for s in matched]
        mag_ab = [magcount_s[s].get(ko, 0) for s in matched]
        rho, _ = spearman(pic_ab, mag_ab)
        marq.append((ko, gene, pic_prev, mag_prev, st.median([v for v in pic_ab]),
                     st.median([v for v in mag_ab]), rho, desc))
    with open(os.path.join(OUT, "marquee_functions.tsv"), "w") as o:
        o.write("KO\tgene\tpicrust2_prevalence\tmag_prevalence\tpic_median_pred\tmag_median_copies\tspearman_abund\tdescription\n")
        for row in marq:
            o.write(f"{row[0]}\t{row[1]}\t{row[2]:.3f}\t{row[3]:.3f}\t{row[4]:.1f}\t{row[5]:.1f}\t{row[6]:.3f}\t{row[7]}\n")

    # ---- raw richness counts + function-level abundance agreement ----
    # PICRUSt2 emits a non-zero value for the same 10,543 KOs in every matched
    # sample under this rule. Treat it as a checked constant, not a median.
    pic_ko_values = {int(r["n_pic_ko"]) for r in rows}
    if len(pic_ko_values) != 1:
        raise ValueError(f"n_pic_ko is not constant: {len(pic_ko_values)} values")
    pic_ko_constant = pic_ko_values.pop()
    med_mag_ko = int(st.median([r["n_mag_ko"] for r in rows]))
    # marquee abundance agreement (signature CO/osmolyte/UV functions encoded widely)
    marq_real = [row for row in marq if row[3] >= 0.5]
    marq_rhos = [row[6] for row in marq_real if row[6] == row[6]]
    marq_rho_min = min(marq_rhos) if marq_rhos else float("nan")
    marq_rho_med = st.median(marq_rhos) if marq_rhos else float("nan")
    marq_rho_max = max(marq_rhos) if marq_rhos else float("nan")

    # The original zero-threshold richness contrast is not interpretable
    # without the shared-threshold rebuild.
    with open(SWEEP) as f:
        sweep = list(csv.DictReader(f, delimiter="\t"))
    with open(REBUILT) as f:
        rebuild_n = sum(1 for _ in csv.DictReader(f, delimiter="\t"))
    asymmetric_ratio = {
        float(r["threshold"]): float(r["ratio_picrust_over_mag_presence"])
        for r in sweep
    }
    shared_relab_ratio = {
        float(r["threshold"]): float(r["ratio_picrust_over_mag_relab"])
        for r in sweep
    }
    asymmetric_zero = asymmetric_ratio[0.0]
    asymmetric_1e5 = asymmetric_ratio[1e-5]
    shared_zero = shared_relab_ratio[0.0]
    shared_1e5 = shared_relab_ratio[1e-5]

    # ---- NSTI summary on matched set ----
    nv = sorted(nsti[s] for s in matched if s in nsti)
    nsti_med = st.median(nv); nsti_p95 = nv[int(len(nv)*0.95)]

    # ---- SUMMARY ----
    L = []
    L.append("# C1 — Reference bias: PICRUSt2 (16S) vs MAG-encoded gene content\n")
    L.append(f"Matched samples (both 16S-PICRUSt2 and genome-resolved MAGs): **{len(matched)}**. "
             f"GTDB classifies these taxa, so this tests *gene-content* predictability, not placement.\n")
    L.append("## Phylogenetic distance is NOT the failure mode")
    L.append(f"Weighted NSTI median **{nsti_med:.3f}** (p95 {nsti_p95:.3f}) — within the typical-soil range "
             "and well under PICRUSt2's NSTI<2 cutoff. ASVs place fine; the failure is downstream, in gene content.\n")
    L.append("## What PICRUSt2 gets right — and the remaining limits")
    L.append(f"PICRUSt2 *does* recover the housekeeping core: median **{med('recall')*100:.1f}%** presence-recall "
             "of MAG-encoded KOs.\n")
    L.append("**(1) KO-richness differences are threshold-sensitive, not evidence of 2.4-fold enrichment.**")
    L.append(f"The original presence comparison counted every non-zero value: PICRUSt2 has exactly "
             f"**{pic_ko_constant:,} predicted KOs in every matched sample**, whereas MAG-encoded richness has "
             f"median **{med_mag_ko:,}**. In the rebuilt {rebuild_n}-sample sweep, the historical asymmetric "
             f"PIC-threshold/MAG-any-presence ratio is **{asymmetric_zero:.2f}** at zero and "
             f"**{asymmetric_1e5:.2f}** at \\(10^{{-5}}\\). The genuinely shared relative-abundance rule gives "
             f"**{shared_zero:.2f}** and **{shared_1e5:.2f}**, respectively. The incompatible denominator rules "
             "do not support a biological enrichment or reversal claim.\n")
    L.append(f"**(2) Abundance agreement is function-dependent.** Across {len(marq_rhos)} prevalent EQ signature "
             "functions (atmospheric CO `coxLMS`, hydrogenase `hyaB`, osmolyte `otsAB`/`ect`, and UV-repair "
             f"`phr`), PICRUSt2-vs-MAG across-sample Spearman correlations range from "
             f"**{marq_rho_min:.2f} to {marq_rho_max:.2f}** (median **{marq_rho_med:.2f}**). This is "
             "heterogeneous agreement, not a uniform failure.\n")
    L.append("**(3) The KO-space ceiling — the dominant failure.**")
    L.append(f"PICRUSt2 operates only in KEGG-KO space. Across {len(dark_rows):,} MAGs, "
             f"**{prot_frac_no_ko*100:.1f}% of all MAG-encoded proteins have no KO** "
             f"(median per-MAG {med_mag_no_ko*100:.1f}%); {prot_frac_no_egg*100:.1f}% have no eggNOG ortholog at all. "
             "This is a hard upper bound: even a perfect KO predictor represents none of that fraction, and the "
             f"{len(mag_only_global):,} KOs the MAGs encode that PICRUSt2 never predicts in any sample are an "
             "additional reference-database blind spot.\n")
    L.append("## Bottom line — the methodological pillar reinforces the dark-matter pillar")
    L.append("An amplicon-only survey of this biome recovers most of the KO core but cannot represent the no-KO "
             "majority of the proteome, including by construction the structural dark matter central to this paper. "
             "It also does not uniformly recover the sample abundance patterns of the atmospheric/osmolyte/UV "
             "signature functions. The KO-richness contrast must be reported with threshold sensitivity and is not "
             "retained as a biological enrichment.\n")
    L.append("_Outputs: per_sample_concordance.tsv, marquee_functions.tsv, dark_fraction.tsv, "
             "rebuild/threshold_sweep.tsv._")
    open(os.path.join(OUT, "SUMMARY.md"), "w").write("\n".join(L) + "\n")

    print(f"matched={len(matched)} recall_med={med('recall'):.3f} spearman_med={med('spearman'):.3f} "
          f"jaccard_med={med('jaccard'):.3f}")
    print(f"NSTI med={nsti_med:.3f}; proteome no-KO={prot_frac_no_ko*100:.1f}% no-eggNOG={prot_frac_no_egg*100:.1f}%; "
          f"MAG-only global KOs={len(mag_only_global)}")

if __name__ == "__main__":
    main()
