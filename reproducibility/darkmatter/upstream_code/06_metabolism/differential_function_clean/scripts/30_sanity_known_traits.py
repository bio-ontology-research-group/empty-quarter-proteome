#!/usr/bin/env python3
"""Sanity check: do the traits the paper already tested reappear in this unbiased scan?

The paper tested UV repair, osmolytes, trace gases, BGC and AMR. BGC and AMR have no
KEGG representation, so they cannot appear on this basis at all and are reported as such.
The other three are matched by keyword against the kofam KO definitions, and for each
matched KO the script reports whether it entered the tested universe (KOs carried by
5-95 % of species representatives) and, if so, its strongest result on any axis.

A KO can miss the universe in two different ways and they mean opposite things:
  core     carried by > 95 % of species, so it cannot differentiate anything here
  rare     carried by < 5 % of species
Neither is the same as "tested and null", so they are reported separately.

Writes results/known_traits_sanity.tsv.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (RES, KO_MIN_PREV, KO_MAX_PREV, load_abundance, load_ko_sets,
                    build_presence)

TRAITS = {
    "UV repair": r"photolyase|excinuclease|photo-lyase|spore photoproduct",
    "osmolyte": r"ectoine|trehalose|betaine|glucosylglycerol|mannosylglycerate|"
                r"proline/glycine betaine|choline dehydrogenase",
    # \b matters: without it "hydrogenase" matches inside "dehydrogenase" and drags in
    # every oxidoreductase in the catalogue
    "trace gas": r"carbon.monoxide dehydrogenase|\bhydrogenase|methane monooxygenase|"
                 r"\bCO dehydrogenase",
}


def main():
    spec = json.load(open(os.path.join(RES, "feature_definitions.json")))
    names = spec["ko"]["names"]

    tests = pd.concat([pd.read_csv(os.path.join(RES, f"{a}_tests.tsv"), sep="\t")
                       for a in ("compartment", "trip", "chemistry")], ignore_index=True)
    tests = tests[tests["basis"] == "ko"]

    hits = {}
    for trait, pat in TRAITS.items():
        for k, n in names.items():
            if n and pd.Series([n]).str.contains(pat, case=False, regex=True).iloc[0]:
                hits.setdefault(k, []).append(trait)

    ab, _ = load_abundance()
    ko_sets = load_ko_sets()
    pres = build_presence(list(ab.index), {k: ([{k}], 1.0) for k in hits}, ko_sets)
    prev = pres.mean(axis=0)

    rows = []
    for k, traits in sorted(hits.items()):
        p = float(prev[k])
        if p < KO_MIN_PREV:
            status = "not tested (rare, < 5 % of species)"
        elif p > KO_MAX_PREV:
            status = "not tested (core, > 95 % of species)"
        else:
            status = "tested"
        sub = tests[tests["feature"] == k]
        best = sub.loc[sub["effect"].abs().idxmax()] if len(sub) else None
        passing = sub[sub["passes"]] if len(sub) else sub
        # `best` is the largest effect on ANY axis, which for some KOs is a trip or
        # chemistry row that does not pass. The paper reports the compartment axis only,
        # so the compartment result is reported separately: the strongest PASSING
        # compartment row if there is one, otherwise the strongest compartment row.
        comp = sub[sub["axis"] == "compartment"] if len(sub) else sub
        cpass = comp[comp["passes"]] if len(comp) else comp
        pool = cpass if len(cpass) else comp
        cb = pool.loc[pool["effect"].abs().idxmax()] if len(pool) else None
        rows.append({
            "trait": "/".join(traits), "ko": k, "name": names.get(k, ""),
            "species_prevalence": round(p, 4), "status": status,
            "n_axes_passed": int(passing["axis"].nunique()) if len(sub) else 0,
            "best_axis": (best["axis"] if best is not None else ""),
            "best_contrast": (best["contrast"] if best is not None else ""),
            "best_effect": (round(float(best["effect"]), 3) if best is not None else np.nan),
            "best_q": (float(best["q"]) if best is not None else np.nan),
            "within_phylum": (best["within_phylum"] if best is not None else ""),
            "best_passes": (bool(best["passes"]) if best is not None else False),
            "comp_contrast": (cb["contrast"] if cb is not None else ""),
            "comp_effect": (round(float(cb["effect"]), 3) if cb is not None else np.nan),
            "comp_q": (float(cb["q"]) if cb is not None else np.nan),
            "comp_within_phylum": (cb["within_phylum"] if cb is not None else ""),
            "comp_eff_Actinomycetota": (round(float(cb["eff_Actinomycetota"]), 3)
                                        if cb is not None else np.nan),
            "comp_eff_Pseudomonadota": (round(float(cb["eff_Pseudomonadota"]), 3)
                                        if cb is not None else np.nan),
            "comp_passes": (bool(cb["passes"]) if cb is not None else False),
        })

    out = pd.DataFrame(rows).sort_values(["trait", "species_prevalence"],
                                         ascending=[True, False])
    path = os.path.join(RES, "known_traits_sanity.tsv")
    out.to_csv(path, sep="\t", index=False, float_format="%.6g")
    print(out.groupby(["trait", "status"]).size())
    print("\npassing at least one axis:")
    print(out[out["n_axes_passed"] > 0][
        ["trait", "ko", "name", "best_axis", "best_effect", "best_q", "within_phylum"]]
        .to_string(index=False, max_colwidth=45))
    print("\ncompartment axis only (what the paper reports):")
    print(out[out["comp_passes"]][
        ["trait", "ko", "name", "comp_contrast", "comp_effect", "comp_q",
         "comp_within_phylum", "comp_eff_Actinomycetota", "comp_eff_Pseudomonadota"]]
        .to_string(index=False, max_colwidth=45))
    print("\nBGC and AMR have no KEGG module or KO representation; not testable here.")
    print("wrote", path)


if __name__ == "__main__":
    main()
