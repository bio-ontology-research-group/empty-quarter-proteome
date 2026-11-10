#!/usr/bin/env python3
"""Harmonised dark-versus-non-dark export comparison across Thermoplasmatota cohorts.

Inputs are the COMPLETE_VERIFIED producer tables (common Pyrodigal gene calls,
direct UniRef50 2026_01 searches, SignalP 6.0h fast mode, organism group other).
Every estimate uses the same gene caller, search rule and signal-peptide predictor
for the desert and the external genomes.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

BINS = [0, 100, 150, 200, 250, 350, 500, 10**9]
EXCLUDED = {
    "External_Comparators_CrossResource_Duplicate_54",
    "SPIRE_ProGenomes3_SpecI_Excluded_148",
}
SHORT = {
    "EQ_Thermoplasmatota_237": "eq",
    "External_Comparators_Deduplicated_247": "external247",
    "Maeke_2025_Penumbrarchaeia_35": "maeke35",
    "SMAG_Thermoplasmatota_354": "smag354",
    "SPIRE_Thermoplasmatota_1084_DeNovo": "spire1084",
}


def load(root: Path) -> pd.DataFrame:
    genes = pd.read_csv(
        root / "genecall/harmonized_gene_manifest.tsv.gz", sep="\t",
        usecols=["gene_id", "genome_id", "is_complete_orf", "contig_edge_status",
                 "aa_length", "n_internal_stops"])
    uni = pd.read_csv(
        root / "mmseqs/harmonized_uniref50_per_gene.tsv.gz", sep="\t",
        usecols=["gene_id", "uniref50_search_status", "primary_direct_uniref50_dark"])
    sp = pd.read_csv(
        root / "signalp/harmonized_signalp6_predictions.tsv.gz", sep="\t",
        usecols=["gene_id", "signalp_status", "signalp_class", "sp_any"])
    n = len(genes)
    if not (genes.gene_id.is_unique and uni.gene_id.is_unique and sp.gene_id.is_unique):
        raise SystemExit("duplicate gene identifiers")
    df = genes.merge(uni, on="gene_id", validate="1:1").merge(sp, on="gene_id", validate="1:1")
    if len(df) != n or len(uni) != n or len(sp) != n:
        raise SystemExit(f"identifier closure failed: {n} {len(uni)} {len(sp)} {len(df)}")
    if not df.uniref50_search_status.str.startswith("completed").all():
        raise SystemExit("incomplete UniRef50 searches present")
    if set(df.signalp_class.unique()) - {"OTHER", "SP", "LIPO", "TAT", "TATLIPO", "PILIN"}:
        raise SystemExit("unexpected SignalP class")
    df["dark"] = df.primary_direct_uniref50_dark.astype(bool)
    df["sp"] = df.sp_any.astype(bool)
    if not (df.sp == (df.signalp_class != "OTHER")).all():
        raise SystemExit("sp_any disagrees with class")
    df["lenbin"] = np.digitize(df.aa_length, BINS[1:-1])
    df["interior"] = (df.is_complete_orf.astype(bool) & (df.contig_edge_status == "interior")
                      & (df.n_internal_stops == 0))
    return df


def strata_table(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(["genome_id", "lenbin", "dark", "sp"], observed=True).size().unstack(
        ["dark", "sp"], fill_value=0)
    for key in [(True, True), (True, False), (False, True), (False, False)]:
        if key not in g.columns:
            g[key] = 0
    out = pd.DataFrame({
        "a": g[(True, True)], "b": g[(True, False)],
        "c": g[(False, True)], "d": g[(False, False)]}).reset_index()
    nn = out[["a", "b", "c", "d"]].sum(axis=1)
    out["R"] = out.a * out.d / nn
    out["S"] = out.b * out.c / nn
    return out


def per_genome(df: pd.DataFrame) -> pd.DataFrame:
    df = df.assign(dark_sp=df.dark & df.sp)
    g = df.groupby("genome_id").agg(
        n_genes=("dark", "size"), n_dark=("dark", "sum"), n_sp=("sp", "sum"),
        n_dark_sp=("dark_sp", "sum"))
    g["dark_frac"] = g.n_dark / g.n_genes
    g["sp_frac"] = g.n_sp / g.n_genes
    g["sp_in_dark"] = g.n_dark_sp / g.n_dark.where(g.n_dark > 0)
    g["sp_in_nondark"] = (g.n_sp - g.n_dark_sp) / (g.n_genes - g.n_dark)
    g["sp_diff"] = g.sp_in_dark - g.sp_in_nondark
    st = strata_table(df).groupby("genome_id")[["R", "S"]].sum()
    return g.join(st)


def boot(pg: pd.DataFrame, block: str, rng, draws: int):
    """Block bootstrap over `block`; returns statistics per draw."""
    codes, _ = pd.factorize(pg[block])
    k = codes.max() + 1
    cols = ["R", "S", "n_genes", "n_dark", "n_sp", "n_dark_sp"]
    sums = np.zeros((k, len(cols)))
    np.add.at(sums, codes, pg[cols].to_numpy(float))
    res = np.empty((draws, 4))
    for i in range(draws):
        s = sums[rng.integers(0, k, k)].sum(axis=0)
        R, S, n, nd, nsp, ndsp = s
        res[i] = [np.log(R / S), nd / n, ndsp / nd, (nsp - ndsp) / (n - nd)]
    return res


def summarise(pg: pd.DataFrame, rng, draws: int) -> dict:
    b = boot(pg, "block", rng, draws)
    q = lambda x: [float(v) for v in np.percentile(x, [2.5, 97.5])]
    n, nd, nsp, ndsp = (pg[c].sum() for c in ["n_genes", "n_dark", "n_sp", "n_dark_sp"])
    return {
        "n_genomes": int(len(pg)), "n_blocks": int(pg.block.nunique()), "n_genes": int(n),
        "dark_frac_pooled": nd / n, "dark_frac_ci": q(b[:, 1]),
        "dark_frac_genome_median": float(pg.dark_frac.median()),
        "sp_in_dark_pooled": ndsp / nd, "sp_in_dark_ci": q(b[:, 2]),
        "sp_in_nondark_pooled": (nsp - ndsp) / (n - nd), "sp_in_nondark_ci": q(b[:, 3]),
        "sp_diff_genome_median": float(pg.sp_diff.median()),
        "genomes_sp_diff_positive": int((pg.sp_diff > 0).sum()),
        "genomes_sp_diff_defined": int(pg.sp_diff.notna().sum()),
        "mh_or": float(pg.R.sum() / pg.S.sum()), "mh_or_ci": [float(np.exp(v)) for v in q(b[:, 0])],
        "_boot_logor": b[:, 0],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261110)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    gm = pd.read_csv(a.root / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
    gm = gm[~gm.cohort.isin(EXCLUDED)].copy()
    ani = pd.read_csv(a.root / "source-metadata/ani_groups/ani_95_groups_1957.tsv", sep="\t", dtype=str)
    gm = gm.merge(ani[["genome_id", "ani95_group_id"]], on="genome_id", how="left", validate="1:1")
    if gm.ani95_group_id.isna().any():
        raise SystemExit("genomes without ANI group")
    gm["cohort_short"] = gm.cohort.map(SHORT)
    # Resampling block: named genus, or the genome's ANI95 group when genus is unassigned.
    gm["block"] = np.where(gm.tax_genus.fillna("g__").isin(["g__", "", "missing", "none"]),
                           "ani::" + gm.ani95_group_id, gm.tax_genus)
    gm["comp"] = gm.completeness_primary.astype(float)
    gm["cont"] = gm.contamination_primary.astype(float)

    df = load(a.root)
    df = df[df.genome_id.isin(gm.genome_id)]
    missing = set(gm.genome_id) - set(df.genome_id)
    if missing:
        raise SystemExit(f"{len(missing)} genomes without genes")

    eq_orders = set(gm.loc[gm.cohort_short == "eq", "tax_order"]) - {"o__", "missing", "none"}
    rng = np.random.default_rng(a.seed)
    rows, boots = [], {}
    for gene_set, d in [("all_genes", df), ("complete_interior", df[df.interior])]:
        pg = per_genome(d).join(gm.set_index("genome_id"), how="inner")
        pg.reset_index().drop(columns=["fasta_path"], errors="ignore").to_csv(
            a.out / f"per_genome_{gene_set}.tsv", sep="\t", index=False)
        pg["pooled_external"] = pg.cohort_short != "eq"
        subsets = {"all": pg.index == pg.index,
                   "quality_70_5": (pg.comp >= 70) & (pg.cont <= 5),
                   "eq_named_orders": pg.tax_order.isin(eq_orders)}
        for sname, mask in subsets.items():
            sub = pg[mask]
            groups = {c: sub[sub.cohort_short == c] for c in SHORT.values()}
            groups["pooled_external"] = sub[sub.pooled_external]
            for cname, part in groups.items():
                if len(part) < 5 or part.S.sum() == 0:
                    continue
                r = summarise(part, rng, a.draws)
                boots[(gene_set, sname, cname)] = r.pop("_boot_logor")
                rows.append({"gene_set": gene_set, "subset": sname, "cohort": cname, **r})
    res = pd.DataFrame(rows)
    # Ratio of odds ratios, desert versus each external group (independent bootstraps).
    ratios = []
    for (gs, sn, cn), b in boots.items():
        if cn == "eq" or (gs, sn, "eq") not in boots:
            continue
        diff = boots[(gs, sn, "eq")] - b
        pt = res[(res.gene_set == gs) & (res.subset == sn)].set_index("cohort").mh_or
        ratios.append({"gene_set": gs, "subset": sn, "external": cn,
                       "ratio_of_or_eq_over_external": float(pt["eq"] / pt[cn]),
                       "ci_low": float(np.exp(np.percentile(diff, 2.5))),
                       "ci_high": float(np.exp(np.percentile(diff, 97.5)))})
    for c in [c for c in res.columns if c.endswith("_ci")]:
        res[c + "_low"] = res[c].str[0]
        res[c + "_high"] = res[c].str[1]
        res = res.drop(columns=c)
    res.to_csv(a.out / "export_by_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(ratios).to_csv(a.out / "export_or_ratio.tsv", sep="\t", index=False)
    (a.out / "receipt.json").write_text(json.dumps({
        "seed": a.seed, "draws": a.draws, "n_genes": int(len(df)), "n_genomes": int(len(gm)),
        "eq_named_orders": sorted(eq_orders), "length_bins": BINS[:-1],
        "block": "named genus, else ANI95 group", "status": "COMPLETE"}, indent=1))
    print(res.to_string())
    print(pd.DataFrame(ratios).to_string())


if __name__ == "__main__":
    main()
