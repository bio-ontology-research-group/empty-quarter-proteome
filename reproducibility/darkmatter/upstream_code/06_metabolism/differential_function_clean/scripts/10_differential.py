#!/usr/bin/env python3
"""Round two: which functions differ by compartment, trip and chemistry?

Runs one pipeline over two feature bases:
  modules  -- KEGG modules, present in a genome when >= 2/3 of the module's top-level
              STEPS are satisfied, a step being satisfied by any one of its alternative
              KOs in rep_ko_v2.tsv (see 00_build_features.py; the rule is over steps, not
              over the raw KO list, and for a one-step module it reduces to "any one KO")
  ko       -- bare KOs, restricted to those present in 5-95 % of the 5,438 species reps

Three axes, all on abundance-weighted carriage (percent of the mapped community that
encodes the feature):
  compartment  site-paired Surface vs Deep, and PlantRoot vs bulk, Wilcoxon signed-rank
  trip         site x compartment units paired across visits, every trip pair with >= 15
  chemistry    partial Spearman vs S, Ca, Cl, Na, Fe, aridity_index, precip, controlling
               ranked lat and lon (the paper's section 2.5 method, reused verbatim)

Everything that passes an axis is re-run inside Actinomycetota only and inside
Pseudomonadota only; the phylum control passes when the sign of the effect holds in both.

Writes results/{compartment,trip,chemistry}_tests.tsv and results/passing_modules.tsv.

Every test row also carries the raw p, what the presence rule demanded (n_steps,
steps_required, n_kos) and the feature's prevalence across species; compartment rows
additionally carry the two carriage levels being differenced (median_a, median_b) and the
number of sites each within-phylum effect rests on. None of these change a test result;
they exist so a reported effect can be audited without rerunning the scan.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (W, R, RES, MODULE_FRACTION, KO_MIN_PREV, KO_MAX_PREV,
                    EFFECT_FLOOR_PP, RHO_FLOOR, QCUT, PHYLUM_MIN_FRAC, CONTROL_PHYLA,
                    load_abundance, load_ko_sets, load_taxonomy, load_genome_size,
                    load_design, load_env, build_presence, carriage, bh_fdr,
                    partial_spearman, wilcoxon_signed_rank_matrix,
                    site_compartment_trip, compartment_diffs, run_compartment)

FEATURE_JSON = os.path.join(RES, "feature_definitions.json")
CHEM_AXES = ["S", "Ca", "Cl", "Na", "Fe", "aridity_index", "precip"]
TRIP_MIN_UNITS = 15


# ------------------------------------------------------------------ aggregation
# site_compartment_trip, compartment_diffs and run_compartment now live in common.py so
# that the review diagnostics exercise the same pairing code; behaviour is unchanged.
META_COLS = ["n_steps", "steps_required", "n_kos", "prev_species"]


def feature_meta(fdefs, prev):
    """Per-feature description of what the presence test actually demanded.

    n_steps / steps_required make the 2/3 rule auditable feature by feature: a one-step
    module is satisfied by any single alternative KO, which is a KO test wearing a module
    name, and that has to be visible in the results table rather than inferred.
    """
    m = pd.DataFrame({
        "n_steps": {f: len(g) for f, (g, _) in fdefs.items()},
        "steps_required": {f: max(1, int(np.ceil(fr * len(g))))
                           for f, (g, fr) in fdefs.items()},
        "n_kos": {f: len(set().union(*g)) for f, (g, _) in fdefs.items()},
    })
    m["prev_species"] = prev.reindex(m.index)
    return m


def add_meta(df, meta, key=None):
    idx = df.index if key is None else df[key]
    for col in META_COLS:
        df[col] = meta[col].reindex(idx).values
    return df


def run_trip(agg, label):
    """Pair site x compartment units across visits."""
    wide = agg.unstack("trip")                   # index site x compartment
    trips = sorted(wide.columns.get_level_values(1).unique())
    feats = agg.columns
    per_trip = {t: wide.xs(t, axis=1, level=1).reindex(columns=feats) for t in trips}
    out = []
    for i, ta in enumerate(trips):
        for tb in trips[i + 1:]:
            A, B = per_trip[ta], per_trip[tb]
            common = A.dropna(how="all").index.intersection(B.dropna(how="all").index)
            if len(common) < TRIP_MIN_UNITS:
                continue
            diff = A.loc[common] - B.loc[common]
            med = diff.median(axis=0, skipna=True)
            _, p, n = wilcoxon_signed_rank_matrix(diff.values)
            out.append(pd.DataFrame({"feature": feats, "contrast": f"{ta}_vs_{tb}",
                                     "n_units": n, "median_diff_pp": med.values,
                                     "p": p}))
    return out


def run_chemistry(vals, env):
    """Partial Spearman vs each chemistry axis, controlling ranked lat and lon."""
    out = []
    e = env.set_index("sample")
    for axis in CHEM_AXES:
        if axis not in e.columns:
            continue
        sub = e.reindex(vals.columns)
        ok = sub[axis].notna() & sub["lat"].notna() & sub["lon"].notna()
        cols = list(vals.columns[ok.values])
        # phylum-restricted carriage blanks whole samples where the phylum is too rare;
        # those columns must go before the rank transform, not row by row
        usable = vals[cols].notna().all(axis=0)
        cols = [c for c in cols if usable[c]]
        if len(cols) < 30:
            continue
        M = vals[cols].values
        good = np.isfinite(M).all(axis=1)
        y = sub.loc[cols, axis].astype(float).values
        lat = sub.loc[cols, "lat"].astype(float).values
        lon = sub.loc[cols, "lon"].astype(float).values
        rho = np.full(M.shape[0], np.nan)
        pv = np.ones(M.shape[0])
        r_, p_ = partial_spearman(M[good], y, [lat, lon])
        rho[good], pv[good] = r_, p_
        out.append(pd.DataFrame({"feature": vals.index, "contrast": axis,
                                 "n_samples": len(cols), "rho": rho, "p": pv}))
    return out


# ------------------------------------------------------------------ phylum control
def phylum_carriage(ab, pres, tax, phylum, size_w=None):
    mask = pd.Series([tax.get(g, "") == phylum for g in ab.index], index=ab.index)
    if mask.sum() < 20:
        return None
    val, denom = carriage(ab, pres, weights=size_w, genome_mask=mask)
    frac = denom / ab.sum(axis=0)
    val.loc[:, frac.reindex(val.columns).values < PHYLUM_MIN_FRAC] = np.nan
    return val


def main():
    with open(FEATURE_JSON) as fh:
        spec = json.load(fh)

    ab, dropped = load_abundance()
    ko_sets = load_ko_sets()
    tax = load_taxonomy()
    sizes = load_genome_size()
    design = load_design()
    env = load_env()

    genomes = list(ab.index)
    size_w = pd.Series({g: 1.0 / sizes[g] for g in genomes if g in sizes})
    log = {"samples_used": ab.shape[1], "samples_dropped_zero_map": dropped,
           "species": len(genomes),
           "mapped_pct_mean": float(ab.sum(axis=0).mean())}

    comp_rows, trip_rows, chem_rows = [], [], []
    carriage_store = {}

    for basis in ("module", "ko"):
        names = spec[basis].get("names", {})
        if basis == "module":
            fdefs = {m: ([set(s) for s in steps], MODULE_FRACTION)
                     for m, steps in spec["module"]["step_groups"].items()}
        else:
            fdefs = {k: ([set(v)], 1.0) for k, v in spec["ko"]["definitions"].items()}
        pres = build_presence(genomes, fdefs, ko_sets)

        prev = pres.mean(axis=0)
        if basis == "ko":
            keep = (prev >= KO_MIN_PREV) & (prev <= KO_MAX_PREV)
        else:
            keep = (prev > 0) & (prev < 1)       # drop features no genome / every genome has
        pres = pres.loc[:, keep]
        meta = feature_meta(fdefs, prev).reindex(pres.columns)
        log[f"{basis}_defined"] = int(len(fdefs))
        log[f"{basis}_tested"] = int(pres.shape[1])
        print(f"[{basis}] {len(fdefs)} defined -> {pres.shape[1]} tested", flush=True)

        val, _ = carriage(ab, pres)                       # read-weighted
        val_sz, _ = carriage(ab, pres, weights=size_w)    # cell-weighted robustness
        agg = site_compartment_trip(val, design)
        agg_sz = site_compartment_trip(val_sz, design)

        carriage_store[basis] = val
        ph_val = {ph: phylum_carriage(ab, pres, tax, ph) for ph in CONTROL_PHYLA}
        ph_agg = {ph: (site_compartment_trip(v, design) if v is not None else None)
                  for ph, v in ph_val.items()}

        # ---------------- compartment
        base = run_compartment(agg, arms=True)
        sz = run_compartment(agg_sz)
        ph_base = {ph: (run_compartment(a) if a is not None else {})
                   for ph, a in ph_agg.items()}
        for kind, df in base.items():
            df = df.copy()
            df["q"] = bh_fdr(df["p"].values)
            df["basis"] = basis
            df["axis"] = "compartment"
            df["sizenorm_effect"] = sz[kind]["median_diff_pp"].reindex(df.index) if kind in sz else np.nan
            for ph in CONTROL_PHYLA:
                d = ph_base[ph].get(kind)
                df[f"eff_{ph}"] = d["median_diff_pp"].reindex(df.index) if d is not None else np.nan
                # the phylum series is blanked wherever the phylum holds < 1 % of the
                # mapped community, so its effect rests on fewer sites than the whole-
                # community one and cannot be read without its own n
                df[f"n_sites_{ph}"] = d["n_sites"].reindex(df.index) if d is not None else np.nan
            df["name"] = [names.get(f, "") for f in df.index]
            comp_rows.append(add_meta(df, meta).reset_index())

        # ---------------- trip
        base_t = run_trip(agg, basis)
        sz_t = {d["contrast"].iloc[0]: d.set_index("feature")["median_diff_pp"]
                for d in run_trip(agg_sz, basis)}
        ph_t = {ph: {d["contrast"].iloc[0]: d.set_index("feature")["median_diff_pp"]
                     for d in run_trip(a, basis)} if a is not None else {}
                for ph, a in ph_agg.items()}
        for df in base_t:
            df = df.copy()
            df["q"] = bh_fdr(df["p"].values)
            df["basis"] = basis
            df["axis"] = "trip"
            c = df["contrast"].iloc[0]
            df["sizenorm_effect"] = sz_t.get(c, pd.Series(dtype=float)).reindex(df["feature"]).values
            for ph in CONTROL_PHYLA:
                df[f"eff_{ph}"] = ph_t[ph].get(c, pd.Series(dtype=float)).reindex(df["feature"]).values
            df["name"] = [names.get(f, "") for f in df["feature"]]
            trip_rows.append(add_meta(df, meta, "feature"))

        # ---------------- chemistry
        base_c = run_chemistry(val, env)
        sz_c = {d["contrast"].iloc[0]: d.set_index("feature")["rho"]
                for d in run_chemistry(val_sz, env)}
        ph_c = {ph: {d["contrast"].iloc[0]: d.set_index("feature")["rho"]
                     for d in run_chemistry(v, env)} if v is not None else {}
                for ph, v in ph_val.items()}
        for df in base_c:
            df = df.copy()
            df["q"] = bh_fdr(df["p"].values)
            df["basis"] = basis
            df["axis"] = "chemistry"
            c = df["contrast"].iloc[0]
            df["sizenorm_effect"] = sz_c.get(c, pd.Series(dtype=float)).reindex(df["feature"]).values
            for ph in CONTROL_PHYLA:
                df[f"eff_{ph}"] = ph_c[ph].get(c, pd.Series(dtype=float)).reindex(df["feature"]).values
            df["name"] = [names.get(f, "") for f in df["feature"]]
            chem_rows.append(add_meta(df, meta, "feature"))

    comp = pd.concat(comp_rows, ignore_index=True)
    trip = pd.concat(trip_rows, ignore_index=True)
    chem = pd.concat(chem_rows, ignore_index=True)

    # ------------------------------------------------------------- passing table
    comp["effect"] = comp["median_diff_pp"]
    trip["effect"] = trip["median_diff_pp"]
    chem["effect"] = chem["rho"]
    comp["passes"] = (comp["q"] < QCUT) & (comp["effect"].abs() >= EFFECT_FLOOR_PP)
    trip["passes"] = (trip["q"] < QCUT) & (trip["effect"].abs() >= EFFECT_FLOOR_PP)
    chem["passes"] = (chem["q"] < QCUT) & (chem["effect"].abs() >= RHO_FLOOR)

    # a trip feature must additionally hold the same direction in every trip pair tested
    sign_ok = (trip.assign(s=np.sign(trip["effect"]))
                   .groupby(["basis", "feature"])["s"]
                   .agg(lambda s: s.dropna().nunique() <= 1))
    trip["dir_consistent"] = [
        bool(sign_ok.get((b, f), False)) for b, f in zip(trip["basis"], trip["feature"])]
    trip["passes"] = trip["passes"] & trip["dir_consistent"]

    for df in (comp, trip, chem):
        s0 = np.sign(df["effect"]).values
        n_known = np.zeros(len(df), dtype=int)
        n_agree = np.zeros(len(df), dtype=int)
        for ph in CONTROL_PHYLA:
            e = df[f"eff_{ph}"].values
            # an effect of exactly zero means the phylum does not carry the feature at
            # all, so it cannot speak to the direction: that is untestable, not a failure
            known = np.isfinite(e) & (e != 0)
            n_known += known
            n_agree += known & (np.sign(e) == s0)
        df["within_phylum"] = np.where(
            n_known == 0, "untestable",
            np.where(n_agree == n_known,
                     np.where(n_known == len(CONTROL_PHYLA), "pass", "partial"),
                     "fail"))
        df["n_phyla_testable"] = n_known
        df["sizenorm_sign_ok"] = np.sign(df["sizenorm_effect"]).values == s0

    cols_out = ["basis", "feature", "name", "axis", "contrast", "effect", "p", "q",
                "within_phylum", "n_phyla_testable",
                f"eff_{CONTROL_PHYLA[0]}", f"eff_{CONTROL_PHYLA[1]}",
                "sizenorm_effect", "sizenorm_sign_ok", "passes"] + META_COLS
    comp_extra = ["n_sites", "median_a", "median_b"] + [f"n_sites_{ph}" for ph in CONTROL_PHYLA]
    comp.to_csv(os.path.join(RES, "compartment_tests.tsv"), sep="\t", index=False,
                columns=cols_out + comp_extra, float_format="%.6g")
    trip.to_csv(os.path.join(RES, "trip_tests.tsv"), sep="\t", index=False,
                columns=cols_out + ["n_units", "dir_consistent"], float_format="%.6g")
    chem.to_csv(os.path.join(RES, "chemistry_tests.tsv"), sep="\t", index=False,
                columns=cols_out + ["n_samples"], float_format="%.6g")

    keep = pd.concat([comp[comp["passes"]], trip[trip["passes"]], chem[chem["passes"]]],
                     ignore_index=True)
    # one row per feature x axis: keep the strongest contrast within an axis
    keep["abs_eff"] = keep["effect"].abs()
    keep = (keep.sort_values("abs_eff", ascending=False)
                .drop_duplicates(["basis", "feature", "axis"])
                .sort_values(["axis", "abs_eff"], ascending=[True, False]))
    keep.to_csv(os.path.join(RES, "passing_modules.tsv"), sep="\t", index=False,
                columns=cols_out + ["abs_eff"], float_format="%.6g")

    # per-sample carriage of the passing features, so the figures do not rebuild presence
    parts = [carriage_store[b].loc[carriage_store[b].index.intersection(
        keep[keep["basis"] == b]["feature"].unique())] for b in carriage_store]
    parts = [p for p in parts if not p.empty]
    if parts:
        pd.concat(parts).to_csv(os.path.join(RES, "carriage_passing.tsv"), sep="\t",
                                float_format="%.5g")

    for ax in ("compartment", "trip", "chemistry"):
        for b in ("module", "ko"):
            s = keep[(keep["axis"] == ax) & (keep["basis"] == b)]
            log[f"pass_{b}_{ax}_before_phylum"] = int(s["feature"].nunique())
            log[f"pass_{b}_{ax}_after_phylum"] = int(
                s[s["within_phylum"] == "pass"]["feature"].nunique())
            log[f"pass_{b}_{ax}_phylum_partial"] = int(
                s[s["within_phylum"] == "partial"]["feature"].nunique())
    # how much of the trip axis is lost to the "same direction in every pair" rule
    for b in ("module", "ko"):
        s = trip[(trip["basis"] == b) & (trip["q"] < QCUT)
                 & (trip["effect"].abs() >= EFFECT_FLOOR_PP)]
        log[f"trip_{b}_floor_and_q_before_direction"] = int(s["feature"].nunique())
    log["trip_contrasts_tested"] = sorted(trip["contrast"].unique())
    log["compartment_contrasts_tested"] = sorted(comp["contrast"].unique())
    log["chem_n_by_axis"] = {c: int(chem[chem["contrast"] == c]["n_samples"].iloc[0])
                             for c in sorted(chem["contrast"].unique())}
    log["compartment_n_sites"] = {c: int(comp[comp["contrast"] == c]["n_sites"].max())
                                  for c in sorted(comp["contrast"].unique())}
    with open(os.path.join(RES, "run_log.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
