#!/usr/bin/env python3
"""Genome-level verification layer for the §2.3/§2.4 KO- and BGC-based censuses.

Three fully independent methods call the same two capacities on the same 5,438
representatives:
  KO census      protein KO annotation, panel lookup (osmolyte_clean, rep_ko_v2)
  antiSMASH      nucleotide rule-based cluster detection (bgc_clean, 8.0.4 --minimal)
                 -- ectoine only; trehalose is not a BGC product class
  gapsmith       MetaCyc pathway completeness from the DRAFT reconstruction
                 (06_darkmatter/gapsmith_models, 5,438 current reps, VERIFIED
                 2026-08-13). all-Pathways.tbl is written by `gapsmith find`
                 BEFORE gap-filling (log step [1/5] vs [4/5]) -- verified
                 2026-08-15 -- so Prediction/Completeness reflect the annotation,
                 not gap-filled (absent) reactions; filled-added.tsv remains the
                 audit trail. Gap-filled models are NOT used for any capacity claim.

Pathways:
  ectoine synthesis    P101-PWY
  trehalose synthesis  TRESYN-PWY, PWY-881, TREHALOSESYN-PWY, PWY-2622, PWY-2661,
                       PWY-5983, PWY-5985 (biosynthesis I-VII; any TRUE = capacity)

Out (results/):
  triangulation_per_rep.tsv  per rep: the three calls per capacity (+ gapsmith
                             max completeness)
  triangulation_summary.tsv  prevalences, pairwise agreement, Jaccard

Run: python3 scripts/check_1_triangulate.py   (reads 5,438 small .tbl files, ~1 GB)
"""
import csv, os
from multiprocessing import Pool

M   = "/ibex/scratch/projects/c2014/rob/metagenomics"
B   = f"{M}/06_metabolism/pathway_check_clean"
RES = f"{B}/results"
GAP = f"{M}/06_darkmatter/gapsmith_models"
OSMO = f"{M}/06_metabolism/osmolyte_clean/results/rep_traits.tsv"
BGC  = f"{M}/06_metabolism/bgc_clean/results/bgc_per_mag.tsv"

ECTO_PWY = {"P101-PWY"}
TRE_PWYS = {"TRESYN-PWY", "PWY-881", "TREHALOSESYN-PWY", "PWY-2622", "PWY-2661",
            "PWY-5983", "PWY-5985"}

# ---- KO-census calls (§2.3) --------------------------------------------------
ko_ect, ko_tre, annotated = {}, {}, {}
reps = []
for r in csv.DictReader(open(OSMO), delimiter="\t"):
    reps.append(r["rep"])
    annotated[r["rep"]] = r["annotated"] == "1"
    ko_ect[r["rep"]] = r["ectoine"] == "1"
    ko_tre[r["rep"]] = r["trehalose_any"] == "1"
print(f"[tri] reps from KO census: {len(reps)} (annotated {sum(annotated.values())})")

# ---- antiSMASH calls (§2.4) --------------------------------------------------
as_ect = {}
for r in csv.DictReader(open(BGC), delimiter="\t"):
    as_ect[r["mag"]] = "ectoine" in r["product_types"].split(";")
print(f"[tri] antiSMASH ectoine carriers: {sum(as_ect.values())}")

# ---- gapsmith calls ----------------------------------------------------------
def scan(rep):
    p = f"{GAP}/{rep}/{rep}-all-Pathways.tbl"
    if not os.path.exists(p):
        return rep, None, None, None, None
    e_pred = t_pred = False
    e_comp = t_comp = 0.0
    with open(p) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            pid = row["ID"].strip("|")
            if pid in ECTO_PWY:
                if row["Prediction"] == "TRUE":
                    e_pred = True
                e_comp = max(e_comp, float(row["Completeness"] or 0))
            elif pid in TRE_PWYS:
                if row["Prediction"] == "TRUE":
                    t_pred = True
                t_comp = max(t_comp, float(row["Completeness"] or 0))
    return rep, e_pred, e_comp, t_pred, t_comp

with Pool(8) as pool:
    scanned = pool.map(scan, reps, chunksize=50)
gp = {rep: (e, ec, t, tc) for rep, e, ec, t, tc in scanned}
missing = [r for r in reps if gp[r][0] is None]
print(f"[tri] gapsmith models scanned: {len(reps) - len(missing)}; missing: {len(missing)}")

# ---- per-rep table -----------------------------------------------------------
with open(f"{RES}/triangulation_per_rep.tsv", "w") as o:
    o.write("rep\tannotated\tko_ectoine\tas_ectoine\tgs_ectoine\tgs_ectoine_completeness\t"
            "ko_trehalose\tgs_trehalose\tgs_trehalose_completeness\n")
    for r in reps:
        e, ec, t, tc = gp[r]
        o.write(f"{r}\t{int(annotated[r])}\t{int(ko_ect[r])}\t{int(as_ect.get(r, False))}\t"
                f"{'' if e is None else int(e)}\t{'' if ec is None else f'{ec:.1f}'}\t"
                f"{int(ko_tre[r])}\t{'' if t is None else int(t)}\t"
                f"{'' if tc is None else f'{tc:.1f}'}\n")

# ---- summary -----------------------------------------------------------------
def pair(name, a, b, universe):
    A = {r for r in universe if a[r]}
    Bs = {r for r in universe if b[r]}
    both, either = len(A & Bs), len(A | Bs)
    agree = sum(1 for r in universe if a[r] == b[r])
    return [name, len(universe), len(A), len(Bs), both,
            f"{100*agree/len(universe):.1f}", f"{both/either:.3f}" if either else "NA"]

have_gs = [r for r in reps if gp[r][0] is not None]
gs_ect = {r: gp[r][0] for r in have_gs}
gs_tre = {r: gp[r][2] for r in have_gs}
uni_ann = [r for r in have_gs if annotated[r]]

rows = [pair("ectoine: KO vs antiSMASH", ko_ect, as_ect,
             [r for r in reps if annotated[r]]),
        pair("ectoine: KO vs gapsmith", ko_ect, gs_ect, uni_ann),
        pair("ectoine: antiSMASH vs gapsmith", as_ect, gs_ect, have_gs),
        pair("trehalose: KO vs gapsmith", ko_tre, gs_tre, uni_ann)]
with open(f"{RES}/triangulation_summary.tsv", "w") as o:
    o.write("comparison\tn_universe\tn_method1\tn_method2\tn_both\tpct_agree\tjaccard\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
        print("  " + "\t".join(str(x) for x in r))
print("[tri] done")
