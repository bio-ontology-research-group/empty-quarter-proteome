#!/usr/bin/env python3
"""Hostile-referee diagnostics for the compartment axis.

The scan in 10_differential.py is the analysis; this script is the evidence for the
method review in the README. It changes nothing in results/ except its own report, and
every number it prints is either a census of the design or a sensitivity rerun of the
SAME pairing code (imported from common.py, not copied) under one perturbation at a time.

Sections
  1  design census: samples, sites, replicates, paired site-visits per contrast
  2  denominator: is the mapped fraction itself different between compartments?
  3  detection bias: are the abundant genomes of one compartment more complete, or
     simply richer in KOs, than the other's? Either shifts every carriage difference.
  4  the effect distribution's own offset, and how much of the passing list sits within
     half a point of the 5 pp floor
  5  presence-rule sensitivity: the whole scan restricted to reps >= 90 % complete
  6  pairing sensitivity: PlantRoot vs bulk restricted to visits where BOTH Surface and
     Deep were sampled, so "bulk" means the same thing at every visit
  7  multiple testing: BH pooled across the two compartment contrasts instead of within
  8  what the module presence rule actually demanded, and module/KO test duplication
  9  within-phylum magnitudes behind the sign-only control
 10  the chemistry axis against trip (reported for the record; the axis is not used)

Writes results/review_diagnostics.txt.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (W, RES, MODULE_FRACTION, KO_MIN_PREV, KO_MAX_PREV, EFFECT_FLOOR_PP,
                    QCUT, load_abundance, load_ko_sets, load_design, load_env,
                    build_presence, carriage, bh_fdr, wilcoxon_signed_rank_matrix,
                    site_compartment_trip, compartment_diffs, run_compartment)

GENOMEINFO = W + "05_refine/_drep/out/data_tables/genomeInformation.csv"
CONTRASTS = ["plantroot_bulk", "surface_deep"]
HIGH_COMPLETENESS = 90.0

lines = []


def say(s=""):
    print(s, flush=True)
    lines.append(s)


def paired_scalar(series, design, label):
    """Run one per-sample scalar through the analysis' own pairing and test."""
    v = pd.DataFrame([series.values], index=[label], columns=series.index)
    agg = site_compartment_trip(v, design)
    out = {}
    for kind in CONTRASTS:
        diff, a, b = compartment_diffs(agg, kind)
        _, p, n = wilcoxon_signed_rank_matrix(diff.values)
        out[kind] = (float(diff.median().iloc[0]), float(a.median().iloc[0]),
                     float(b.median().iloc[0]), int(n[0]), float(p[0]))
    return out


def main():
    ab, dropped = load_abundance()
    design = load_design()
    ko_sets = load_ko_sets()
    env = load_env()
    spec = json.load(open(os.path.join(RES, "feature_definitions.json")))
    tests = pd.read_csv(os.path.join(RES, "compartment_tests.tsv"), sep="\t")
    d = design.reindex(ab.columns)

    # ---------------------------------------------------------------- 1 design
    say("=" * 78)
    say("1  DESIGN CENSUS")
    say("=" * 78)
    say(f"samples in the abundance matrix: {ab.shape[1] + len(dropped)}, "
        f"used {ab.shape[1]} (dropped for zero mapped abundance: {dropped})")
    say(f"species representatives: {ab.shape[0]}")
    say(f"sites: {d['site'].nunique()}   trips: {sorted(d['trip'].unique())}")
    say("samples per compartment: " +
        ", ".join(f"{k} {v}" for k, v in d["compartment"].value_counts().items()))
    units = d.groupby(["site", "compartment", "trip"]).size()
    say(f"site x compartment x trip units: {len(units)}; replicates per unit: " +
        ", ".join(f"{k} rep(s) x {v}" for k, v in units.value_counts().sort_index().items()))
    say("replicates inside a unit are collapsed by median before any test, so the unit of "
        "analysis is the site, never the sequencing replicate.")
    visits = d.reset_index().groupby(["site", "trip"])["compartment"].apply(set)
    sd = [k for k, v in visits.items() if {"Surface", "Deep"} <= v]
    pb = [k for k, v in visits.items() if "PlantRoot" in v and v & {"Surface", "Deep"}]
    say(f"site-visits pairing Surface with Deep: {len(sd)} over {len({k[0] for k in sd})} sites")
    say(f"site-visits pairing PlantRoot with bulk: {len(pb)} over {len({k[0] for k in pb})} sites")
    both = sum(1 for k in pb if {"Surface", "Deep"} <= visits[k])
    only_s = sum(1 for k in pb if visits[k] & {"Surface", "Deep"} == {"Surface"})
    only_d = sum(1 for k in pb if visits[k] & {"Surface", "Deep"} == {"Deep"})
    say(f"  of those, bulk = Surface and Deep at {both} visits, Surface only at {only_s}, "
        f"Deep only at {only_d}")
    say("sites per contrast in the results table: " +
        ", ".join(f"{c} {int(tests[tests.contrast == c]['n_sites'].max())}" for c in CONTRASTS))

    # ---------------------------------------------------------------- 2 denominator
    mapped = ab.sum(axis=0)
    say()
    say("=" * 78)
    say("2  THE MAPPED-COMMUNITY DENOMINATOR")
    say("=" * 78)
    say("carriage is a percentage of the mapped community, so the denominator is itself a "
        "per-sample quantity. If it differs by compartment, the invisible half of the "
        "community differs too.")
    say(mapped.groupby(d["compartment"]).agg(["count", "mean", "median", "min", "max"])
        .round(2).to_string())
    r = paired_scalar(mapped, design, "mapped_pct")
    for k, (diff, a, b, n, p) in r.items():
        say(f"paired {k:15s} A {a:6.2f} % - B {b:6.2f} % = {diff:+6.2f} pp "
            f"over {n} sites, Wilcoxon p = {p:.3g}")

    # ---------------------------------------------------------------- 3 detection bias
    gi = pd.read_csv(GENOMEINFO)
    gi["g"] = gi["genome"].str.replace(r"\.fa$", "", regex=True)
    comp_g = pd.Series(dict(zip(gi["g"], gi["completeness"]))).reindex(ab.index)
    nko = pd.Series({g: len(ko_sets.get(g, ())) for g in ab.index})
    say()
    say("=" * 78)
    say("3  DETECTION BIAS IN THE ABUNDANT FRACTION")
    say("=" * 78)
    say(f"rep completeness: median {comp_g.median():.2f}, "
        f"{int((comp_g >= HIGH_COMPLETENESS).sum())} of {len(comp_g)} at or above "
        f"{HIGH_COMPLETENESS:.0f} %. Presence is read off incomplete genomes, so carriage "
        "is biased low; the bias only distorts a CONTRAST if it differs between arms.")
    for nm, ser in (("abundance-weighted completeness", comp_g),
                    ("abundance-weighted KOs per genome", nko)):
        w = ab.mul(ser, axis=0).sum(axis=0) / ab.sum(axis=0)
        say(f"{nm}: " + ", ".join(
            f"{k} {v:.2f}" for k, v in w.groupby(d["compartment"]).mean().items()))
        for k, (diff, a, b, n, p) in paired_scalar(w, design, nm).items():
            say(f"   paired {k:15s} {a:8.2f} - {b:8.2f} = {diff:+7.2f} over {n} sites, "
                f"p = {p:.3g}")

    # ---------------------------------------------------------------- 4 offset and floor
    say()
    say("=" * 78)
    say("4  WHERE THE EFFECT DISTRIBUTION SITS, AND HOW CLOSE THE PASSING LIST IS TO THE FLOOR")
    say("=" * 78)
    for basis in ("module", "ko"):
        for c in CONTRASTS:
            s = tests[(tests.basis == basis) & (tests.contrast == c)]
            med = s["effect"].median()
            npass = s["passes"].sum()
            near = int(((s["passes"]) & (s["effect"].abs() < EFFECT_FLOOR_PP + 0.5)).sum())
            cen = s[(s["q"] < QCUT) & ((s["effect"] - med).abs() >= EFFECT_FLOOR_PP)]
            lost = len(set(s[s.passes]["feature"]) - set(cen["feature"]))
            gained = len(set(cen["feature"]) - set(s[s.passes]["feature"]))
            say(f"{basis:6s} {c:15s} n={len(s):5d} median effect {med:+.3f} pp, "
                f"{(s['effect'] > 0).mean() * 100:4.1f} % positive; passing {npass:4d}, "
                f"of which {near} within 0.5 pp of the floor; "
                f"floor applied to median-centred effects: {len(cen)} pass "
                f"({lost} lost, {gained} gained)")

    # ---------------------------------------------------------------- 5 and 6 reruns
    say()
    say("=" * 78)
    say(f"5  SENSITIVITY: REPS >= {HIGH_COMPLETENESS:.0f} % COMPLETE ONLY")
    say("=" * 78)
    hc_mask = (comp_g >= HIGH_COMPLETENESS).fillna(False)
    say(f"{int(hc_mask.sum())} reps retained; they hold "
        f"{(ab[hc_mask.values].sum(axis=0) / ab.sum(axis=0) * 100).mean():.1f} % of mapped "
        "abundance on average (min "
        f"{(ab[hc_mask.values].sum(axis=0) / ab.sum(axis=0) * 100).min():.1f} %).")
    say()
    say("=" * 78)
    say("6  SENSITIVITY: PlantRoot vs bulk ON VISITS WITH BOTH Surface AND Deep")
    say("=" * 78)
    keep_visits = set(sd)

    for basis in ("module", "ko"):
        if basis == "module":
            fdefs = {m: ([set(s) for s in st], MODULE_FRACTION)
                     for m, st in spec["module"]["step_groups"].items()}
        else:
            fdefs = {k: ([set(v)], 1.0) for k, v in spec["ko"]["definitions"].items()}
        pres = build_presence(list(ab.index), fdefs, ko_sets)
        prev = pres.mean(axis=0)
        keep = ((prev >= KO_MIN_PREV) & (prev <= KO_MAX_PREV)) if basis == "ko" \
            else ((prev > 0) & (prev < 1))
        pres = pres.loc[:, keep]

        val_hc, _ = carriage(ab, pres, genome_mask=hc_mask)
        hc = run_compartment(site_compartment_trip(val_hc, design))
        val, _ = carriage(ab, pres)
        agg = site_compartment_trip(val, design)
        idx = [i for i in agg.index if (i[0], i[2]) in keep_visits]
        sub_agg = agg.loc[idx]
        bulk_both = run_compartment(sub_agg)

        for c in CONTRASTS:
            base = tests[(tests.basis == basis) & (tests.contrast == c)].set_index("feature")
            pset = base[base["passes"]]
            for label, alt in (("5 high-completeness", hc.get(c)),
                               ("6 both-bulk visits", bulk_both.get(c) if c == "plantroot_bulk" else None)):
                if alt is None or pset.empty:
                    continue
                e = alt["median_diff_pp"].reindex(pset.index)
                q = pd.Series(bh_fdr(alt["p"].values), index=alt.index).reindex(pset.index)
                same = int((np.sign(e) == np.sign(pset["effect"])).sum())
                floor = int((e.abs() >= EFFECT_FLOOR_PP).sum())
                still = int(((e.abs() >= EFFECT_FLOOR_PP) & (q < QCUT) &
                             (np.sign(e) == np.sign(pset["effect"]))).sum())
                ratio = (e.abs() / pset["effect"].abs()).median()
                say(f"[{label}] {basis:6s} {c:15s} {len(pset):4d} passing features: "
                    f"sign held {same}, |effect| still >= {EFFECT_FLOOR_PP:.0f} pp {floor}, "
                    f"still passing floor+FDR {still}, median |effect| ratio {ratio:.2f}, "
                    f"n_sites {int(alt['n_sites'].max())}")
                if basis == "module":
                    drop = pset.index[~((e.abs() >= EFFECT_FLOOR_PP) & (q < QCUT))]
                    say("          drops below floor or FDR: " +
                        ", ".join(f"{f} ({pset.loc[f, 'effect']:+.1f} -> {e[f]:+.1f} pp"
                                  f"{', phylum ' + pset.loc[f, 'within_phylum']})"
                                  for f in drop))

    # ---------------------------------------------------------------- 7 pooled BH
    say()
    say("=" * 78)
    say("7  BH POOLED ACROSS THE TWO COMPARTMENT CONTRASTS")
    say("=" * 78)
    say("the run corrects within (basis x contrast). Pooling the two contrasts of one "
        "basis into a single family is the stricter alternative:")
    for basis in ("module", "ko"):
        s = tests[tests.basis == basis].copy()
        s["q_pooled"] = bh_fdr(s["p"].values)
        now = s["passes"].sum()
        pooled = int(((s["q_pooled"] < QCUT) & (s["effect"].abs() >= EFFECT_FLOOR_PP)).sum())
        lost = s[s["passes"] & ~((s["q_pooled"] < QCUT) &
                                 (s["effect"].abs() >= EFFECT_FLOOR_PP))]
        say(f"{basis:6s} rows passing now {now}, pooled {pooled}, lost {len(lost)} "
            f"({sorted(set(lost['feature']))[:8]}{' ...' if len(set(lost['feature'])) > 8 else ''})")
        say(f"       distinct features now {s[s.passes]['feature'].nunique()}, "
            f"pooled {s[(s['q_pooled'] < QCUT) & (s['effect'].abs() >= EFFECT_FLOOR_PP)]['feature'].nunique()}")

    # ---------------------------------------------------------------- 8 presence rule
    say()
    say("=" * 78)
    say("8  WHAT THE MODULE PRESENCE RULE DEMANDED")
    say("=" * 78)
    mod = tests[tests.basis == "module"].drop_duplicates("feature")
    say(f"tested modules {len(mod)}: steps per module median {mod['n_steps'].median():.0f}, "
        f"one-step {int((mod['n_steps'] == 1).sum())}, two or three steps "
        f"{int(mod['n_steps'].between(2, 3).sum())}, four or more "
        f"{int((mod['n_steps'] >= 4).sum())}")
    ko_universe = set(tests[tests.basis == "ko"]["feature"])
    one = mod[mod["n_steps"] == 1]
    dup = [f for f in one["feature"]
           if len(spec["module"]["definitions"][f]) == 1
           and spec["module"]["definitions"][f][0] in ko_universe]
    say(f"one-step modules are satisfied by ANY one of their alternative KOs "
        f"(up to {int(one['n_kos'].max())} alternatives), so they are KO-level tests "
        f"carrying a pathway name; {len(dup)} of them are a single KO that is also tested "
        f"in the KO family, i.e. the identical test appears in both BH families: {dup}")
    kot = tests[tests.basis == "ko"].set_index(["feature", "contrast"])
    for f in dup:
        k = spec["module"]["definitions"][f][0]
        for c in CONTRASTS:
            m = tests[(tests.feature == f) & (tests.contrast == c)].iloc[0]
            o = kot.loc[(k, c)]
            say(f"   {f} == {k} on {c:15s}: module effect {m['effect']:+.3f} pp q {m['q']:.3g} "
                f"passes {bool(m['passes'])} | KO effect {o['effect']:+.3f} pp q {o['q']:.3g} "
                f"passes {bool(o['passes'])}")
    passing = tests[tests["passes"] & (tests.basis == "module")]
    for c in CONTRASTS:
        s = passing[passing.contrast == c]
        say(f"passing modules on {c}: {len(s)}, of which one-step "
            f"{int((s['n_steps'] == 1).sum())} "
            f"({sorted(s[s['n_steps'] == 1]['feature'])})")

    # ---------------------------------------------------------------- 9 phylum magnitudes
    say()
    say("=" * 78)
    say("9  WITHIN-PHYLUM MAGNITUDES BEHIND THE SIGN-ONLY CONTROL")
    say("=" * 78)
    for basis in ("module", "ko"):
        s = tests[(tests.basis == basis) & tests["passes"] &
                  (tests["within_phylum"] == "pass")].copy()
        if s.empty:
            continue
        sm = s[["eff_Actinomycetota", "eff_Pseudomonadota"]].abs().min(axis=1)
        lg = s[["eff_Actinomycetota", "eff_Pseudomonadota"]].abs().max(axis=1)
        say(f"{basis:6s} {len(s)} phylum-passing rows: smaller within-phylum |effect| is a "
            f"median {(sm / s['effect'].abs()).median():.2f} of the whole-community "
            f"|effect| (larger one {(lg / s['effect'].abs()).median():.2f}); "
            f"{int(((sm / s['effect'].abs()) < 0.2).sum())} rows below 0.2, "
            f"{int((sm >= EFFECT_FLOOR_PP).sum())} rows where BOTH phyla clear the same "
            f"{EFFECT_FLOOR_PP:.0f} pp floor on their own")
        say(f"       sites behind the within-phylum effects: Actinomycetota median "
            f"{s['n_sites_Actinomycetota'].median():.0f}, Pseudomonadota median "
            f"{s['n_sites_Pseudomonadota'].median():.0f}, against "
            f"{s['n_sites'].median():.0f} for the whole community")

    # ---------------------------------------------------------------- 10 chemistry vs trip
    say()
    say("=" * 78)
    say("10  THE CHEMISTRY AXIS IS A TRIP AXIS (not reported; recorded for the file)")
    say("=" * 78)
    e = env.set_index("sample").reindex(ab.columns)
    x = pd.crosstab(e["trip"], e["S"].notna())
    x["coverage %"] = (100 * x[True] / (x[True] + x[False])).round(1)
    say(x.to_string())
    say("XRF coverage runs from 32 % of trip1 samples to 100 % of trip4, and a trip is an "
        "extraction and sequencing batch, so an element-carriage correlation cannot be "
        "separated from batch on this design.")

    path = os.path.join(RES, "review_diagnostics.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
