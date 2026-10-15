"""Thermoplasmatota ecology (occurrence, compartment, depth pairing) and clean
functional absences on the KO table, against Halobacteriota / all archaea / all bacteria.

Reads only; writes three tables under results/:
  v2_thermo_ecology_summary.tsv    long format, one row per reported number
  v2_thermo_ecology_by_sample.tsv  one row per profiled sample column
  v2_thermo_ecology_by_rep.tsv     one row per representative, marker presence 0/1

Denominators.
  5,438 species representatives = 591 archaea + 4,847 bacteria.
  237 Thermoplasmatota = 185 c__SW-10-69-26 + 51 c__E2 + 1 c__Thermoplasmata.
  269 Halobacteriota (all c__Halobacteria).
  The abundance matrix carries 550 sample columns of the 551 canonical samples
  (F13PRr1 is absent from it). F5Dr1 is entirely literal `nan`, so 549 are usable.
  All `nan` cells are coerced to absent (0), as in the other consumers of this matrix.

KO identities were checked against the KEGG REST API on 2026-09-20; see MARKERS below.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")

F_REPS = os.path.join(RES, "v2_reps.tsv")
F_AB = W + "/10_thermal/abundance/rep_relative_abundance.tsv"
F_META = W + "/metadata/sample_metadata.tsv"
F_KO = W + "/04_functional/repko_clean/results/rep_ko_v2.tsv"
F_BAC = W + "/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
F_AR = W + "/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"

OUT_SUM = os.path.join(RES, "v2_thermo_ecology_summary.tsv")
OUT_SAMP = os.path.join(RES, "v2_thermo_ecology_by_sample.tsv")
OUT_REP = os.path.join(RES, "v2_thermo_ecology_by_rep.tsv")

# KEGG names verified 2026-09-20 via https://rest.kegg.jp/list/<ko>.
# Rhodopsins: the KEGG set was taken from `find/ko/rhodopsin`; of the microbial
# entries only K04641 bop, K04642 hop and K04643 sop occur in rep_ko_v2.tsv.
# K23306 (brz, bacteriorhodopsin-regulating zinc finger) is absent from the table
# and is a regulator, not a rhodopsin. K04646 was rejected: KEGG calls it CLTC,
# clathrin heavy chain, not sensory rhodopsin II.
MARKERS = [
    ("mcr", "K00399", "mcrA methyl-coenzyme M reductase alpha"),
    ("mcr", "K00401", "mcrB methyl-coenzyme M reductase beta"),
    ("mcr", "K00402", "mcrG methyl-coenzyme M reductase gamma"),
    ("amo", "K10944", "pmoA-amoA ammonia/methane monooxygenase A"),
    ("amo", "K10945", "pmoB-amoB ammonia/methane monooxygenase B"),
    ("amo", "K10946", "pmoC-amoC ammonia/methane monooxygenase C"),
    ("rhodopsin", "K04641", "bop bacteriorhodopsin"),
    ("rhodopsin", "K04642", "hop halorhodopsin"),
    ("rhodopsin", "K04643", "sop sensory rhodopsin"),
    ("hydrogenase", "K06281", "hyaB/hybC NiFe hydrogenase large subunit"),
    ("coxL", "K03520", "coxL form-I CO dehydrogenase large subunit"),
    ("nif", "K02588", "nifH nitrogenase iron protein"),
    ("nif", "K02586", "nifD nitrogenase MoFe alpha"),
    ("nif", "K02591", "nifK nitrogenase MoFe beta"),
    ("photolyase", "K01669", "phr deoxyribodipyrimidine photolyase"),
    ("trehalose", "K00697", "otsA trehalose-6-phosphate synthase"),
    ("trehalose", "K01087", "otsB trehalose-6-phosphate phosphatase"),
    ("ectoine", "K06720", "ectC L-ectoine synthase"),
    ("aerobic", "K02274", "coxA cytochrome c oxidase subunit I"),
    ("aerobic", "K00425", "cydA cytochrome bd subunit I"),
    ("aerobic", "K00426", "cydB cytochrome bd subunit II"),
]
SETS = [
    ("mcr", ["K00399", "K00401", "K00402"], "methanogenesis / anaerobic methanotrophy"),
    ("amo", ["K10944", "K10945", "K10946"], "ammonia (methane) monooxygenase"),
    ("rhodopsin", ["K04641", "K04642", "K04643"], "microbial rhodopsins"),
    ("hydrogenase", ["K06281"], "NiFe hydrogenase large subunit"),
    ("coxL", ["K03520"], "form-I CO dehydrogenase"),
    ("nif", ["K02588", "K02586", "K02591"], "nitrogenase"),
    ("photolyase", ["K01669"], "photolyase"),
    ("trehalose", ["K00697", "K01087"], "trehalose synthesis"),
    ("ectoine", ["K06720"], "ectoine synthase"),
    ("aerobic_coxA", ["K02274"], "cytochrome c oxidase subunit I"),
    ("aerobic_cyd", ["K00425", "K00426"], "cytochrome bd"),
]

ROWS = []


def add(section, group, metric, value, numerator="", denominator="", note=""):
    ROWS.append(dict(section=section, group=group, metric=metric, value=value,
                     numerator=numerator, denominator=denominator, note=note))


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return format(round(float(x), nd), "." + str(nd) + "f")


def rank(classification, letter):
    if not isinstance(classification, str):
        return "unclassified"
    for part in classification.split(";"):
        part = part.strip()
        if part.startswith(letter + "__"):
            return part[3:] or "unclassified"
    return "unclassified"


# ---------------------------------------------------------------- inputs
def load_reps():
    reps = pd.read_csv(F_REPS, sep="\t")
    tax = pd.concat([
        pd.read_csv(F_BAC, sep="\t", usecols=["user_genome", "classification"]),
        pd.read_csv(F_AR, sep="\t", usecols=["user_genome", "classification"]),
    ]).rename(columns={"user_genome": "genome"})
    tax["class"] = tax["classification"].map(lambda s: rank(s, "c"))
    reps = reps.merge(tax[["genome", "class"]], on="genome", how="left")
    reps["class"] = reps["class"].fillna("unclassified")
    return reps


def load_abundance():
    ab = pd.read_csv(F_AB, sep="\t", index_col=0)
    ab = ab.apply(pd.to_numeric, errors="coerce")
    dead = [c for c in ab.columns if ab[c].isna().all() or (ab[c].fillna(0) == 0).all()]
    ab = ab.fillna(0.0)
    return ab, dead


def main():
    reps = load_reps()
    n_all = len(reps)
    ab, dead_cols = load_abundance()
    meta = pd.read_csv(F_META, sep="\t")

    n_cols = ab.shape[1]
    usable = [c for c in ab.columns if c not in dead_cols]
    n_use = len(usable)
    assert n_all == 5438 and ab.shape[0] == 5438

    add("inputs", "all", "species_representatives", str(n_all), str(n_all), str(n_all))
    add("inputs", "all", "abundance_sample_columns", str(n_cols), str(n_cols), "551",
        "551 canonical samples; F13PRr1 has no column")
    add("inputs", "all", "usable_sample_columns", str(n_use), str(n_use), str(n_cols),
        "dropped all-nan/all-zero: " + ",".join(dead_cols))

    groups = {
        "Thermoplasmatota": reps["phylum"] == "Thermoplasmatota",
        "Halobacteriota": reps["phylum"] == "Halobacteriota",
        "Archaea": reps["domain"] == "Archaea",
        "Bacteria": reps["domain"] == "Bacteria",
    }
    # v2_reps.tsv and the abundance matrix are not in the same row order, so every
    # group is selected by genome name, never by position.
    members = {n: set(reps.loc[m, "genome"]) for n, m in groups.items()}
    assert set(ab.index) == set(reps["genome"]), "abundance index and v2_reps disagree"
    for name, mask in groups.items():
        add("groups", name, "n_representatives", str(int(mask.sum())),
            str(int(mask.sum())), str(n_all))
        add("groups", name, "median_completeness", fmt(reps.loc[mask, "completeness"].median(), 2),
            "", str(int(mask.sum())), "IQR " + fmt(reps.loc[mask, "completeness"].quantile(.25), 2)
            + " to " + fmt(reps.loc[mask, "completeness"].quantile(.75), 2))
    for cls, sub in reps[groups["Thermoplasmatota"]].groupby("class"):
        add("groups", "Thermoplasmatota", "n_class_" + cls, str(len(sub)), str(len(sub)), "237")

    # ------------------------------------------------------- per sample
    meta_i = meta.set_index("sample")
    samp = pd.DataFrame(index=pd.Index(ab.columns, name="sample"))
    samp["usable"] = [c not in dead_cols for c in ab.columns]
    for col in ["site", "compartment", "trip"]:
        samp[col] = [meta_i[col].get(c, "unknown") for c in ab.columns]
    for name in groups:
        block = ab.loc[ab.index.isin(members[name])]
        assert len(block) == len(members[name])
        samp[name.lower() + "_sum_pct"] = block.sum(axis=0).values
        samp[name.lower() + "_n_reps_detected"] = (block > 0).sum(axis=0).values
    samp["total_mapped_pct"] = ab.sum(axis=0).values
    # Share of the mapped fraction, not of all reads: total mappability itself moves
    # between compartments, so the raw summed share confounds the two.
    for name in groups:
        samp[name.lower() + "_frac_of_mapped"] = np.where(
            samp["total_mapped_pct"] > 0,
            samp[name.lower() + "_sum_pct"] / samp["total_mapped_pct"], np.nan)
    samp = samp.reset_index()
    samp.to_csv(OUT_SAMP, sep="\t", index=False, float_format="%.6f")

    use = samp[samp["usable"]].copy()

    def quant(series, section, group, metric, denom, note=""):
        s = series.dropna()
        add(section, group, metric + "_median", fmt(s.median(), 4), "", str(denom), note)
        add(section, group, metric + "_q25", fmt(s.quantile(.25), 4), "", str(denom), note)
        add(section, group, metric + "_q75", fmt(s.quantile(.75), 4), "", str(denom), note)

    for name in groups:
        key = name.lower() + "_sum_pct"
        det = int((use[key] > 0).sum())
        add("occurrence", name, "samples_with_any_rep_detected", str(det), str(det), str(n_use))
        quant(use[key], "occurrence", name, "summed_share_pct", n_use)
        quant(use[name.lower() + "_frac_of_mapped"], "occurrence", name,
              "frac_of_mapped", n_use, "summed share divided by total mapped per cent")
        nrk = name.lower() + "_n_reps_detected"
        quant(use[nrk], "occurrence", name, "n_reps_detected_per_sample", n_use)

    for comp, sub in use.groupby("compartment"):
        for name in groups:
            key = name.lower() + "_sum_pct"
            det = int((sub[key] > 0).sum())
            add("by_compartment", name, "samples_with_any_rep_detected_" + comp,
                str(det), str(det), str(len(sub)))
            quant(sub[key], "by_compartment", name, "summed_share_pct_" + comp, len(sub))
    for trip, sub in use.groupby("trip"):
        for name in groups:
            key = name.lower() + "_sum_pct"
            det = int((sub[key] > 0).sum())
            add("by_trip", name, "samples_with_any_rep_detected_" + trip,
                str(det), str(det), str(len(sub)))
            quant(sub[key], "by_trip", name, "summed_share_pct_" + trip, len(sub))

    # ------------------------------------------------- site-paired tests
    def paired(group_name, label, mask_a, label_b, mask_b, suffix="_sum_pct", tag=""):
        key = group_name.lower() + suffix
        a = use[mask_a].groupby("site")[key].median()
        b = use[mask_b].groupby("site")[key].median()
        a = a.dropna()
        b = b.dropna()
        sites = sorted(set(a.index) & set(b.index))
        n_sites = len(sites)
        add("paired", group_name, "sites_with_both_" + label + "_vs_" + label_b + tag,
            str(n_sites), str(n_sites), str(use["site"].nunique()))
        if n_sites < 3:
            return
        va = a.loc[sites].values
        vb = b.loc[sites].values
        diff = va - vb
        higher_a = int((diff > 0).sum())
        higher_b = int((diff < 0).sum())
        tie = int((diff == 0).sum())
        try:
            st, p = stats.wilcoxon(va, vb)
        except ValueError:
            st, p = float("nan"), float("nan")
        direction = label + " higher" if np.median(diff) > 0 else (
            label_b + " higher" if np.median(diff) < 0 else "tie")
        base = label + "_vs_" + label_b + tag
        add("paired", group_name, base + "_median_" + label, fmt(np.median(va)), "", str(n_sites))
        add("paired", group_name, base + "_median_" + label_b, fmt(np.median(vb)), "", str(n_sites))
        add("paired", group_name, base + "_median_paired_difference", fmt(np.median(diff)),
            "", str(n_sites))
        add("paired", group_name, base + "_sites_" + label + "_higher",
            str(higher_a), str(higher_a), str(n_sites), str(tie) + " ties")
        add("paired", group_name, base + "_sites_" + label_b + "_higher",
            str(higher_b), str(higher_b), str(n_sites), str(tie) + " ties")
        add("paired", group_name, base + "_wilcoxon_nonzero_pairs",
            str(n_sites - tie), str(n_sites - tie), str(n_sites),
            "scipy drops zero-difference pairs; this is the effective n")
        add("paired", group_name, base + "_wilcoxon_statistic", fmt(st, 1), "", str(n_sites))
        add("paired", group_name, base + "_wilcoxon_p", "%.3e" % p, "", str(n_sites))
        add("paired", group_name, base + "_direction", direction, "", str(n_sites),
            "signed-rank on per-site medians of the summed share")

    is_s = use["compartment"] == "Surface"
    is_d = use["compartment"] == "Deep"
    is_pr = use["compartment"] == "PlantRoot"
    is_bulk = is_s | is_d
    for name in groups:
        paired(name, "Surface", is_s, "Deep", is_d)
        paired(name, "PlantRoot", is_pr, "bulk", is_bulk)
        paired(name, "Surface", is_s, "Deep", is_d, "_frac_of_mapped", "_fracmapped")
        paired(name, "PlantRoot", is_pr, "bulk", is_bulk, "_frac_of_mapped", "_fracmapped")

    # -------------------------------------------------------- KO markers
    ko = pd.read_csv(F_KO, sep="\t", header=None, names=["genome", "ko_list"])
    ko_sets = {g: set(str(s).split(",")) for g, s in zip(ko["genome"], ko["ko_list"])}
    rep = reps[["genome", "domain", "phylum", "class", "completeness", "contamination",
                "dark_fraction"]].copy()
    rep["has_ko_row"] = rep["genome"].map(lambda g: g in ko_sets)
    for _, kid, _ in MARKERS:
        rep[kid] = rep["genome"].map(lambda g, k=kid: int(k in ko_sets.get(g, set())))
    for sname, kids, _ in SETS:
        rep["any_" + sname] = rep[kids].max(axis=1)
    rep["group"] = np.where(rep["phylum"] == "Thermoplasmatota", "Thermoplasmatota",
                            np.where(rep["phylum"] == "Halobacteriota", "Halobacteriota",
                                     rep["domain"]))
    rep.to_csv(OUT_REP, sep="\t", index=False, float_format="%.6f")

    add("inputs", "all", "reps_with_a_ko_row", str(int(rep["has_ko_row"].sum())),
        str(int(rep["has_ko_row"].sum())), str(n_all))

    comp_groups = {
        "Thermoplasmatota": rep["phylum"] == "Thermoplasmatota",
        "Halobacteriota": rep["phylum"] == "Halobacteriota",
        "Archaea": rep["domain"] == "Archaea",
        "Bacteria": rep["domain"] == "Bacteria",
    }
    for gname, gmask in comp_groups.items():
        sub = rep[gmask]
        n = len(sub)
        for sect, kid, desc in MARKERS:
            c = int(sub[kid].sum())
            add("ko_marker", gname, kid, str(c), str(c), str(n), desc)
        for sname, kids, desc in SETS:
            c = int(sub["any_" + sname].sum())
            add("ko_set_any", gname, sname, str(c), str(c), str(n), desc + " (any of " +
                "+".join(kids) + ")")
            call = int((sub[kids].min(axis=1) == 1).sum()) if len(kids) > 1 else c
            add("ko_set_all", gname, sname, str(call), str(call), str(n),
                desc + " (all of " + "+".join(kids) + ")")
        add("ko_marker", gname, "median_completeness", fmt(sub["completeness"].median(), 2),
            "", str(n), "read absences against this")

    # Thermoplasmatota aerobic markers, per class, for the record.
    thermo = rep[rep["phylum"] == "Thermoplasmatota"]
    for cls, sub in thermo.groupby("class"):
        for kid in ["K02274", "K00425", "K00426"]:
            c = int(sub[kid].sum())
            add("aerobic_by_class", "Thermoplasmatota " + cls, kid, str(c), str(c), str(len(sub)))
        c = int(sub["any_aerobic_cyd"].sum())
        add("aerobic_by_class", "Thermoplasmatota " + cls, "any_cyd", str(c), str(c), str(len(sub)))

    out = pd.DataFrame(ROWS)[["section", "group", "metric", "value", "numerator",
                              "denominator", "note"]]
    out.to_csv(OUT_SUM, sep="\t", index=False)
    sys.stderr.write("wrote %s (%d rows), %s (%d rows), %s (%d rows)\n" % (
        OUT_SUM, len(out), OUT_SAMP, len(samp), OUT_REP, len(rep)))


if __name__ == "__main__":
    main()
