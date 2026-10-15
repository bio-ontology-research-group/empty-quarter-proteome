#!/usr/bin/env python3
"""v2_thermo_relatives - do the 237 desert Thermoplasmatota species have any known relative, and how close.

Reads only results/ files (the downloads in tmp/thermo_relatives/ are deleted after the run):
  results/v2_thermo_ecology_by_rep.tsv          class of the 237 at r232
  results/v2_thermo_relatives_comparators.tsv   manifest of the 218 GTDB r232 reps + 83 MGnify
  results/v2_thermo_relatives_skani.tsv.gz      skani 0.3.1 triangle, sparse (ANI >= ~80 only)
  results/v2_thermo_relatives_aai.tsv.gz        RBH AAI for every genome pair with >= 1 RBH
  results/v2_thermo_relatives_gtdbtk_mgnify.tsv GTDB-Tk 2.7.2 / r232 on the MGnify genomes
  W/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv  closest r232 reference for the 237 (cross-check)

Thresholds: species ANI >= 95 % (reported with and without AF >= 50 %); AAI >= 65 % is a GENUS-LEVEL
HEURISTIC (Konstantinidis, Rosselló-Móra & Amann 2017 range 65-95 % for genus), not a GTDB rule.
An AAI resting on fewer than MIN_RBH reciprocal pairs is not used for a nearest-relative call.

Outputs: results/v2_thermo_relatives_pairs.tsv, _summary.tsv, _report.txt
"""
import os, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
GTDBTK_DESERT = f"{W}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
TYPE = "GCA_003009755.1"   # c__SW-10-69-26 class-name genome (Atacama halite)
MIN_RBH = 100
GENUS_AAI = 65.0
out_lines = []


def say(s=""):
    print(s)
    out_lines.append(s)


def gid_from_path(p):
    return re.sub(r"\.(fa|fna|fna\.gz)$", "", os.path.basename(p))


def rank(lin, letter):
    for part in str(lin).split(";"):
        if part.startswith(letter + "__"):
            return part[3:]
    return ""


# ------------------------------------------------------------------ genome table
eco = pd.read_csv(os.path.join(RES, "v2_thermo_ecology_by_rep.tsv"), sep="\t",
                  usecols=["genome", "phylum", "class"])
des = eco[eco.phylum.str.replace("p__", "") == "Thermoplasmatota"].rename(
    columns={"genome": "genome_id", "class": "cls"})
des = des.assign(source="desert", lineage="")[["genome_id", "source", "cls", "lineage"]]
assert len(des) == 237

man = pd.read_csv(os.path.join(RES, "v2_thermo_relatives_comparators.tsv"), sep="\t", dtype=str)
man = man[man.status == "ok"].copy()
gt = pd.read_csv(os.path.join(RES, "v2_thermo_relatives_gtdbtk_mgnify.tsv"), sep="\t", dtype=str)
gt = gt.rename(columns={"user_genome": "genome_id", "classification": "lineage_r232"})
man = man.merge(gt[["genome_id", "lineage_r232"]], on="genome_id", how="left")
man["lineage_r232"] = np.where(man.source == "GTDB_r232_rep", man.lineage, man.lineage_r232)
man["cls_r232"] = man.lineage_r232.map(lambda s: rank(s, "c"))
comp = man[["genome_id", "source", "cls", "lineage", "lineage_r232", "cls_r232", "biome",
            "catalogue", "same_assembly_as_gtdb_rep"]]
info = pd.concat([des, comp], ignore_index=True).set_index("genome_id")
DES = set(des.genome_id)
GTDB = set(man.genome_id[man.source == "GTDB_r232_rep"])
MGN = set(man.genome_id[man.source == "MGnify"])
say(f"genomes: desert {len(DES)}, GTDB r232 reps {len(GTDB)}, MGnify {len(MGN)}")

# ------------------------------------------------------------------ skani (symmetrise)
sk = pd.read_csv(os.path.join(RES, "v2_thermo_relatives_skani.tsv.gz"), sep="\t")
sk["a"] = sk.Ref_file.map(gid_from_path)
sk["b"] = sk.Query_file.map(gid_from_path)
assert set(sk.a) | set(sk.b) <= set(info.index), "unmapped skani genome"
ski = pd.concat([
    sk.rename(columns={"a": "x", "b": "y", "Align_fraction_ref": "af_x",
                       "Align_fraction_query": "af_y"})[["x", "y", "ANI", "af_x", "af_y"]],
    sk.rename(columns={"b": "x", "a": "y", "Align_fraction_query": "af_x",
                       "Align_fraction_ref": "af_y"})[["x", "y", "ANI", "af_x", "af_y"]]])
say(f"skani pairs reported (ANI >= ~80 %, AF >= 15 %): {len(sk)}")

# ------------------------------------------------------------------ AAI (symmetrise)
aa = pd.read_csv(os.path.join(RES, "v2_thermo_relatives_aai.tsv.gz"), sep="\t")
aai = pd.concat([
    aa.rename(columns={"genome_a": "x", "genome_b": "y"}),
    aa.rename(columns={"genome_b": "x", "genome_a": "y", "n_prot_a": "n_prot_b",
                       "n_prot_b": "n_prot_a"})])[["x", "y", "aai", "n_rbh", "frac_rbh_of_smaller"]]
say(f"AAI genome pairs with >= 1 RBH: {len(aa)}; with >= {MIN_RBH} RBH: {(aa.n_rbh >= MIN_RBH).sum()}")
aai_ok = aai[aai.n_rbh >= MIN_RBH]


def nearest(tab, x_set, y_set, col):
    t = tab[tab.x.isin(x_set) & tab.y.isin(y_set)]
    return t.sort_values(col, ascending=False).drop_duplicates("x").set_index("x")


# MGnify genomes that duplicate a GTDB rep: same NCBI assembly (manifest flag) or skani ANI >= 99.9
# with AF >= 90 both ways. They are dropped from every nearest-relative call (the GTDB rep stands for
# them) and lend their MGnify biome to the GTDB rep they duplicate.
dup = ski[ski.x.isin(MGN) & ski.y.isin(GTDB) & (ski.ANI >= 99.9) & (ski.af_x >= 90) & (ski.af_y >= 90)]
same_asm = set(man.genome_id[(man.source == "MGnify") & (man.same_assembly_as_gtdb_rep == "yes")])
MGN_DUP = same_asm | set(dup.x)
twin = dict(zip(dup.y, dup.x))
MGN_ONLY = MGN - MGN_DUP
say(f"MGnify duplicates of a GTDB rep: {len(MGN_DUP)} (same assembly {len(same_asm)}, skani >= 99.9 % "
    f"{dup.x.nunique()}, skani-only {len(set(dup.x) - same_asm)}); distinct MGnify comparators {len(MGN_ONLY)}")
for gg, mm in twin.items():
    info.at[gg, "biome"] = f"{info.at[mm, 'biome']} (via {mm})"

COMP = GTDB | MGN_ONLY
nani = nearest(ski, DES, COMP, "ANI")
naai = nearest(aai_ok, DES, COMP, "aai")
naai_g = nearest(aai_ok, DES, GTDB, "aai")
naai_m = nearest(aai_ok, DES, MGN_ONLY, "aai")
ndes = nearest(aai_ok[aai_ok.x != aai_ok.y], DES, DES, "aai")
to_type = aai[(aai.y == TYPE) & aai.x.isin(DES)].set_index("x")

gd = pd.read_csv(GTDBTK_DESERT, sep="\t", dtype=str).set_index("user_genome")

rows = []
for g in sorted(DES):
    r = dict(genome=g, cls=info.at[g, "cls"],
             lineage_r232=gd.at[g, "classification"] if g in gd.index else "")
    if g in nani.index:
        y = nani.loc[g]
        r.update(ani_nearest=y.y, ani=y.ANI, af_desert=y.af_x, af_comparator=y.af_y,
                 ani_source=info.at[y.y, "source"], ani_lineage=info.at[y.y, "lineage_r232"])
    if g in naai.index:
        y = naai.loc[g]
        r.update(aai_nearest=y.y, aai=y.aai, aai_n_rbh=y.n_rbh, aai_frac_rbh=y.frac_rbh_of_smaller,
                 aai_source=info.at[y.y, "source"], aai_lineage=info.at[y.y, "lineage_r232"],
                 aai_biome=info.at[y.y, "biome"])
    if g in naai_g.index:
        r.update(aai_nearest_gtdb=naai_g.loc[g].y, aai_gtdb=naai_g.loc[g].aai)
    if g in naai_m.index:
        r.update(aai_nearest_mgnify=naai_m.loc[g].y, aai_mgnify=naai_m.loc[g].aai)
    if g in to_type.index:
        r.update(aai_to_type_GCA_003009755=to_type.loc[g].aai, n_rbh_to_type=to_type.loc[g].n_rbh)
    if g in ndes.index:
        r.update(aai_nearest_desert=ndes.loc[g].y, aai_desert=ndes.loc[g].aai)
    if g in gd.index:
        r.update(gtdbtk_closest_ref=gd.at[g, "closest_genome_reference"],
                 gtdbtk_closest_ani=gd.at[g, "closest_genome_ani"],
                 gtdbtk_closest_af=gd.at[g, "closest_genome_af"],
                 gtdbtk_closest_taxonomy=gd.at[g, "closest_genome_taxonomy"])
    rows.append(r)
COLS = ["genome", "cls", "lineage_r232",
        "ani_nearest", "ani", "af_desert", "af_comparator", "ani_source", "ani_lineage",
        "aai_nearest", "aai", "aai_n_rbh", "aai_frac_rbh", "aai_source", "aai_lineage", "aai_biome",
        "aai_nearest_gtdb", "aai_gtdb", "aai_nearest_mgnify", "aai_mgnify",
        "aai_to_type_GCA_003009755", "n_rbh_to_type", "aai_nearest_desert", "aai_desert",
        "gtdbtk_closest_ref", "gtdbtk_closest_ani", "gtdbtk_closest_af", "gtdbtk_closest_taxonomy"]
pairs = pd.DataFrame(rows).reindex(columns=COLS)
for c in ["ani", "af_desert", "af_comparator", "aai", "aai_gtdb", "aai_mgnify",
          "aai_to_type_GCA_003009755", "aai_desert", "aai_frac_rbh"]:
    pairs[c] = pd.to_numeric(pairs[c], errors="coerce")
pairs["order_r232"] = pairs.lineage_r232.map(lambda s: rank(s, "o"))
pairs["family_r232"] = pairs.lineage_r232.map(lambda s: rank(s, "f"))
pairs["genus_r232"] = pairs.lineage_r232.map(lambda s: rank(s, "g"))
pairs["nearest_same_family"] = [bool(f) and rank(l, "f") == f
                                for f, l in zip(pairs.family_r232, pairs.aai_lineage)]
pairs.to_csv(os.path.join(RES, "v2_thermo_relatives_pairs.tsv"), sep="\t", index=False,
             float_format="%.3f")

# ------------------------------------------------------------------ summary
BINS = [0, 50, 55, 60, 65, 70, 80, 95, 101]
LAB = ["<50", "50-55", "55-60", "60-65", "65-70", "70-80", "80-95", ">=95"]
summ = []


def add(group, metric, num, den, note=""):
    summ.append(dict(group=group, metric=metric, value=num, denominator=den, note=note))


groups = [("all_237", pairs), ("SW-10-69-26", pairs[pairs.cls == "SW-10-69-26"]),
          ("E2", pairs[pairs.cls == "E2"]), ("Thermoplasmata", pairs[pairs.cls == "Thermoplasmata"])]
for name, t in groups:
    n = len(t)
    add(name, "n_reps", n, n)
    add(name, "any_skani_hit_to_comparator (ANI>=~80)", int(t.ani.notna().sum()), n)
    add(name, "ANI>=95_to_comparator", int((t.ani >= 95).sum()), n)
    add(name, "ANI>=95_and_AF>=50_to_comparator", int(((t.ani >= 95) & (t.af_desert >= 50)
                                                       & (t.af_comparator >= 50)).sum()), n)
    add(name, f"AAI>={GENUS_AAI:g}_to_comparator (genus-level heuristic)",
        int((t.aai >= GENUS_AAI).sum()), n)
    add(name, f"AAI>={GENUS_AAI:g}_to_GTDB_r232_rep", int((t.aai_gtdb >= GENUS_AAI).sum()), n)
    add(name, f"AAI>={GENUS_AAI:g}_to_MGnify", int((t.aai_mgnify >= GENUS_AAI).sum()), n)
    add(name, "nearest_by_AAI_is_MGnify", int((t.aai_source == "MGnify").sum()), n)
    add(name, "no_comparator_with_enough_RBH", int(t.aai.isna().sum()), n, f"n_rbh<{MIN_RBH}")
    q = t.aai.dropna()
    if len(q):
        for lab, v in zip(["min", "q25", "median", "q75", "max"],
                          np.percentile(q, [0, 25, 50, 75, 100])):
            add(name, f"best_AAI_{lab}", round(float(v), 3), len(q))
        cut = pd.cut(q, BINS, labels=LAB, right=False).value_counts().reindex(LAB)
        for lab, v in cut.items():
            add(name, f"best_AAI_bin_{lab}", int(v), len(q))
    d = t.aai_desert.dropna()
    if len(d):
        add(name, "median_AAI_to_nearest_other_desert_rep", round(float(d.median()), 3), len(d))

sw = pairs[pairs.cls == "SW-10-69-26"]
add("SW-10-69-26", f"nearest_by_AAI_is_{TYPE}", int((sw.aai_nearest == TYPE).sum()), len(sw))
tt = sw.aai_to_type_GCA_003009755.dropna()
for lab, v in zip(["min", "median", "max"], np.percentile(tt, [0, 50, 100])):
    add("SW-10-69-26", f"AAI_to_{TYPE}_{lab}", round(float(v), 3), len(tt))

# GTDB reps: nearest species-representative relative (desert or other GTDB rep)
POOL = DES | GTDB
ng = nearest(aai_ok[aai_ok.x != aai_ok.y], GTDB, POOL, "aai")
ng["y_is_desert"] = ng.y.isin(DES)
ng["genus_level"] = ng.aai >= GENUS_AAI
for c in ["SW-10-69-26", "E2", "all"]:
    sub = ng if c == "all" else ng[[info.at[x, "cls"] == c for x in ng.index]]
    den = len(GTDB) if c == "all" else int((man.cls[man.source == "GTDB_r232_rep"] == c).sum())
    add(f"GTDB_reps_{c}", "nearest_species_rep_is_desert", int(sub.y_is_desert.sum()), den)
    add(f"GTDB_reps_{c}", f"nearest_is_desert_and_AAI>={GENUS_AAI:g}",
        int((sub.y_is_desert & sub.genus_level).sum()), den)
    add(f"GTDB_reps_{c}", f"any_desert_rep_at_AAI>={GENUS_AAI:g}",
        int(aai_ok[aai_ok.x.isin(sub.index) & aai_ok.y.isin(DES)
                   & (aai_ok.aai >= GENUS_AAI)].x.nunique()), den)

# AAI calibration on GTDB r232 rep pairs, by the lowest rank they share
gl = man[man.source == "GTDB_r232_rep"].set_index("genome_id").lineage


def shared(a, b):
    la, lb = gl[a].split(";"), gl[b].split(";")
    for lvl, name in [(5, "same_genus"), (4, "same_family"), (3, "same_order"), (2, "same_class")]:
        if la[lvl] == lb[lvl] and len(la[lvl]) > 3:
            return name
    return "different_class"


cal = aa[aa.genome_a.isin(GTDB) & aa.genome_b.isin(GTDB) & (aa.n_rbh >= MIN_RBH)].copy()
cal["shared"] = [shared(a, b) for a, b in zip(cal.genome_a, cal.genome_b)]
CAL = ["same_genus", "same_family", "same_order", "same_class", "different_class"]
for c in CAL:
    x = cal.aai[cal.shared == c]
    if len(x):
        q05, q50, q95 = np.percentile(x, [5, 50, 95])
        add("calibration_GTDB_pairs", f"{c}_AAI_median", round(float(q50), 3), len(x),
            f"q05 {q05:.2f} q95 {q95:.2f}; frac>=65 {(x >= GENUS_AAI).mean():.3f}")
        if c == "same_family":
            STRICT = round(float(q95), 2)   # excludes 95 % of same-family, different-genus pairs
for name, t in groups:
    add(name, f"AAI>={STRICT:g}_to_comparator (same-family q95 cut)", int((t.aai >= STRICT).sum()), len(t))
top = pairs.sort_values("aai", ascending=False).head(3)
for _, r in top.iterrows():
    add("top_pairs", f"{r.genome} -> {r.aai_nearest}", r.aai, int(r.aai_n_rbh),
        f"{r.lineage_r232} | {r.aai_lineage} | {r.aai_source} {r.aai_biome}")

# desert placement depth and nearest-family agreement
add("all_237", "r232_genus_assigned", int((pairs.genus_r232 != "").sum()), 237)
add("all_237", "r232_family_assigned", int((pairs.family_r232 != "").sum()), 237)
add("all_237", "nearest_by_AAI_in_same_r232_family", int(pairs.nearest_same_family.sum()),
    int((pairs.family_r232 != "").sum()))
add("all_237", "median_RBH_fraction_to_nearest", round(float(pairs.aai_frac_rbh.median()), 4), 237)
osw = sw[sw.lineage_r232.str.contains("o__SW-10-69-26")]
add("SW-10-69-26_in_o__SW-10-69-26", "n_reps", len(osw), len(osw))
add("SW-10-69-26_in_o__SW-10-69-26", f"nearest_by_AAI_is_{TYPE}", int((osw.aai_nearest == TYPE).sum()),
    len(osw))
if len(osw):
    add("SW-10-69-26_in_o__SW-10-69-26", f"AAI_to_{TYPE}_median",
        round(float(osw.aai_to_type_GCA_003009755.median()), 3), len(osw),
        f"min {osw.aai_to_type_GCA_003009755.min():.2f} max {osw.aai_to_type_GCA_003009755.max():.2f}")

# MGnify at r232
mg = man[man.source == "MGnify"]
for c, v in mg.cls_r232.fillna("unclassified").value_counts().items():
    add("MGnify", f"class_r232_{c}", int(v), len(mg))
add("MGnify", "class_changed_vs_catalogue", int((mg.cls != mg.cls_r232).sum()), len(mg))
add("MGnify", "same_assembly_as_GTDB_r232_rep", int((mg.same_assembly_as_gtdb_rep == "yes").sum()),
    len(mg))
for b, v in mg.biome.value_counts().items():
    add("MGnify", f"biome_{b}", int(v), len(mg))
S = pd.DataFrame(summ)
S.to_csv(os.path.join(RES, "v2_thermo_relatives_summary.tsv"), sep="\t", index=False)


# ------------------------------------------------------------------ report
def v(group, metric):
    r = S[(S.group == group) & (S.metric == metric)]
    if not len(r):
        return ("NA", "NA")
    x = r.value.iloc[0]
    x = int(x) if float(x).is_integer() else x
    return (x, int(r.denominator.iloc[0]))


say("")
say("Known relatives of the 237 desert Thermoplasmatota species representatives")
say(f"Comparators: {len(GTDB)} GTDB r232 species reps (c__SW-10-69-26 + c__E2) and {len(MGN)} MGnify "
    "genomes in the same two classes. SMAG not used (Zenodo 8223844 is one split archive).")
say(f"AAI is reciprocal-best-hit AAI (Pyrodigal + DIAMOND); nearest calls need >= {MIN_RBH} RBH pairs. "
    f"AAI >= {GENUS_AAI:g} % is a genus-level heuristic, not a GTDB criterion.")
say("")
say(f"{'metric':58s} {'all':>9s} {'SW':>9s} {'E2':>9s} {'Tpl':>7s}")
for m in ["any_skani_hit_to_comparator (ANI>=~80)", "ANI>=95_to_comparator",
          "ANI>=95_and_AF>=50_to_comparator",
          f"AAI>={GENUS_AAI:g}_to_comparator (genus-level heuristic)",
          f"AAI>={GENUS_AAI:g}_to_GTDB_r232_rep", f"AAI>={GENUS_AAI:g}_to_MGnify",
          "nearest_by_AAI_is_MGnify", "best_AAI_q25", "best_AAI_median", "best_AAI_q75",
          "best_AAI_max", "median_AAI_to_nearest_other_desert_rep"] + [f"best_AAI_bin_{l}" for l in LAB]:
    vals = [v(g, m)[0] for g in ["all_237", "SW-10-69-26", "E2", "Thermoplasmata"]]
    say(f"{m:58s} " + " ".join(f"{str(x):>9s}" if i < 3 else f"{str(x):>7s}"
                              for i, x in enumerate(vals)))
say("")
a, d = v("SW-10-69-26", f"nearest_by_AAI_is_{TYPE}")
say(f"SW-10-69-26 reps whose nearest comparator by AAI is {TYPE} (Atacama halite): {a} / {d}; "
    f"AAI to it min {v('SW-10-69-26', f'AAI_to_{TYPE}_min')[0]}, median "
    f"{v('SW-10-69-26', f'AAI_to_{TYPE}_median')[0]}, max {v('SW-10-69-26', f'AAI_to_{TYPE}_max')[0]}.")
a, d = v("SW-10-69-26_in_o__SW-10-69-26", f"nearest_by_AAI_is_{TYPE}")
r = S[(S.group == "SW-10-69-26_in_o__SW-10-69-26") & S.metric.str.startswith("AAI_to")]
say(f"  of the {d} SW-10-69-26 reps placed in o__SW-10-69-26 (the type genome's order), {a} have it as "
    f"nearest; AAI to it median {r.value.iloc[0]} ({r.note.iloc[0]}). Other orders: " +
    "; ".join(f"{k or 'none'} {n}" for k, n in sw[~sw.index.isin(osw.index)].order_r232.value_counts().items()))
say(f"  desert reps with an r232 genus: {v('all_237', 'r232_genus_assigned')[0]} / 237; with a family: "
    f"{v('all_237', 'r232_family_assigned')[0]}; nearest-by-AAI comparator in the same r232 family: "
    f"{v('all_237', 'nearest_by_AAI_in_same_r232_family')[0]} / "
    f"{v('all_237', 'nearest_by_AAI_in_same_r232_family')[1]}; median RBH fraction to nearest "
    f"{v('all_237', 'median_RBH_fraction_to_nearest')[0]}")
nt = sw.aai_nearest.value_counts().head(5)
say("  most frequent nearest comparators of the SW-10-69-26 reps: " +
    "; ".join(f"{k} ({info.at[k, 'lineage_r232'].split(';g__')[-1] if isinstance(info.at[k, 'lineage_r232'], str) else ''}) {n}"
              for k, n in nt.items()))
ne = pairs[pairs.cls == "E2"].aai_nearest.value_counts().head(5)
say("  most frequent nearest comparators of the E2 reps: " +
    "; ".join(f"{k} ({str(info.at[k, 'lineage_r232']).split(';g__')[-1]}) {n}" for k, n in ne.items()))
say(f"  AAI >= {STRICT:g} (same-family q95 on this pipeline, a stricter genus-level cut): " +
    ", ".join(f"{g} {v(g, f'AAI>={STRICT:g}_to_comparator (same-family q95 cut)')[0]}"
              for g in ["all_237", "SW-10-69-26", "E2", "Thermoplasmata"]))
say("  highest AAI pairs: ")
for _, r in S[S.group == "top_pairs"].iterrows():
    say(f"    {r.metric} AAI {r.value:.2f} ({r.denominator} RBH): {r.note}")
say("  The single Thermoplasmata rep (o__RBG-16-68-12) has no same-class comparator in this set; its "
    "best AAI is a between-class value, not evidence of no relative.")
say("")
for c in ["all", "SW-10-69-26", "E2"]:
    a1, den = v(f"GTDB_reps_{c}", "nearest_species_rep_is_desert")
    a2, _ = v(f"GTDB_reps_{c}", f"nearest_is_desert_and_AAI>={GENUS_AAI:g}")
    a3, _ = v(f"GTDB_reps_{c}", f"any_desert_rep_at_AAI>={GENUS_AAI:g}")
    say(f"GTDB r232 reps ({c}, n={den}): nearest species rep (desert or GTDB, >= {MIN_RBH} RBH) is a "
        f"desert rep for {a1}; desert AND AAI >= {GENUS_AAI:g}: {a2}; any desert rep at AAI >= "
        f"{GENUS_AAI:g}: {a3}")
say("")
say("AAI calibration on this pipeline, GTDB r232 rep pairs by lowest shared rank (median; q05-q95; "
    "share >= 65):")
for _, r in S[S.group == "calibration_GTDB_pairs"].iterrows():
    say(f"  {r.metric.replace('_AAI_median', ''):16s} n={r.denominator:<6d} {r.value}  {r.note}")
say("")
say("MGnify genomes at GTDB r232 (GTDB-Tk 2.7.2): " +
    "; ".join(f"{m.replace('class_r232_', '')} {int(x)}" for m, x in
              S[(S.group == 'MGnify') & S.metric.str.startswith('class_r232_')][['metric', 'value']].values))
say(f"  class changed vs catalogue lineage: {v('MGnify', 'class_changed_vs_catalogue')[0]}; same NCBI "
    f"assembly as a GTDB r232 rep: {v('MGnify', 'same_assembly_as_GTDB_r232_rep')[0]} / {len(mg)}")
say("  biomes: " + "; ".join(f"{m.replace('biome_root:', '')} {int(x)}" for m, x in
                             S[(S.group == 'MGnify') & S.metric.str.startswith('biome_')][['metric', 'value']].values))
nb = pairs.aai_biome.fillna("").map(lambda s: re.sub(r" \(via .*\)", "", s) or "no biome (GTDB only)")
say(f"  MGnify duplicates of GTDB reps excluded from nearest calls: {len(MGN_DUP)}; distinct: {len(MGN_ONLY)}")
say("  biome of the nearest comparator of the 237 (GTDB reps take the biome of their MGnify twin): " +
    "; ".join(f"{k.replace('root:', '')} {n}" for k, n in nb.value_counts().items()))
say("  catalogues: " +"; ".join(f"{k} {n}" for k, n in mg.catalogue.value_counts().items()))
say("")
say("Cross-check against the GTDB-Tk r232 classification of the 237 (closest reference over ALL r232):")
has = pairs.gtdbtk_closest_ref.notna() & (pairs.gtdbtk_closest_ref != "N/A")
say(f"  reps with a GTDB-Tk closest reference (ANI screen): {int(has.sum())} / 237")
if has.any():
    out = pairs[has & ~pairs.gtdbtk_closest_taxonomy.fillna("").str.contains("c__SW-10-69-26|c__E2")]
    say(f"  of which outside c__SW-10-69-26/c__E2: {len(out)}")
    for _, r in pairs[has].iterrows():
        say(f"    {r.genome} -> {r.gtdbtk_closest_ref} ANI {r.gtdbtk_closest_ani} AF "
            f"{r.gtdbtk_closest_af} | skani nearest {r.get('ani_nearest')} {r.get('ani')}")
hi = pairs[pairs.ani.notna()]
if len(hi):
    say("  desert reps with any skani hit (ANI >= ~80 %) to a comparator:")
    for _, r in hi.iterrows():
        say(f"    {r.genome} -> {r.ani_nearest} ({r.ani_source}) ANI {r.ani:.2f} AF "
            f"{r.af_desert:.1f}/{r.af_comparator:.1f} {r.ani_lineage}")

with open(os.path.join(RES, "v2_thermo_relatives_report.txt"), "w") as o:
    o.write("\n".join(out_lines) + "\n")
