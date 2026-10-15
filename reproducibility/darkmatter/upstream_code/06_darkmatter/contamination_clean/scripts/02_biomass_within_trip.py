#!/usr/bin/env python3
"""Within-trip biomass relations, so the biomass test cannot be a trip confound.

The global Spearman of E. coli against mapped fraction is near zero only because the trips differ
in both prevalence and depth (Simpson's paradox). Stratifying by trip is the honest test.

Two independent biomass proxies are used and they are uncorrelated with each other
(Spearman rho 0.009, n = 187), so agreement between them is real evidence:
  mapped_pct  = 100 - unmapped_percent, available for all 550 samples
  dna_ng_ul   = input DNA concentration, available for 188 samples (trips 1 to 3 only)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean"
RES = os.path.join(P, "results")

TARGETS = {
    "V42PRr2__single_SemiBin_33": "Escherichia coli",
    "F4Sr2__single_SemiBin_30_sub": "Acinetobacter baumannii",
    "F1PRr2__self_SemiBin_738": "Enterococcus_B lactis",
    "30PRr3__single_SemiBin_4_sub": "Bacillus_A anthracis",
}
COMPARATORS = {
    "F13PRr2__self_SemiBin_135": "Brevundimonas sp. (F13)",
    "F27PRr2__self_SemiBin_29": "Brevundimonas sp. (F27)",
    "S34PRr1__self_SemiBin_99": "Brevundimonas sp. (S34)",
    "30PRr3__self_SemiBin_293": "Metabacillus niabensis",
}
TRIPS = ["trip1", "trip2", "trip3", "trip4", "trip5"]

ab = pd.read_csv(os.path.join(W, "10_thermal/abundance/rep_relative_abundance.tsv"),
                 sep="\t", index_col=0).fillna(0.0)
meta = pd.read_csv(os.path.join(W, "metadata/sample_metadata.tsv"),
                   sep="\t").set_index("sample").reindex(ab.columns)
des = pd.read_csv(os.path.join(W, "10_thermal/abundance/sample_design.tsv"),
                  sep="\t", index_col="sample_id").reindex(ab.columns)
dna = pd.read_csv(os.path.join(W, "metadata/dna_conc.tsv"),
                  sep="\t", index_col="sample")["dna_ng_ul"].reindex(ab.columns)
mapped = 100.0 - des["unmapped_percent"]

rows = []
for g, sp in {**TARGETS, **COMPARATORS}.items():
    v = ab.loc[g]
    rec = {"genome": g, "species": sp,
           "role": "target" if g in TARGETS else "comparator"}
    rhos_m = []
    for t in TRIPS:
        m = (meta["trip"] == t) & (mapped > 0)
        # A trip in which the genome is never detected gives a constant vector, for which
        # Spearman is undefined. Such trips are NOT evaluable and must not be counted as
        # agreeing with any direction; n_trips_evaluable records how many really contributed.
        detected_here = int((v[m] > 0).sum())
        r = (stats.spearmanr(v[m], mapped[m])[0]
             if m.sum() >= 12 and detected_here >= 1 else np.nan)
        rec[f"rho_mapped_{t}"] = r
        rec[f"n_det_{t}"] = detected_here
        if not np.isnan(r):
            rhos_m.append(r)
        d = (meta["trip"] == t) & dna.notna() & (mapped > 0)
        rec[f"rho_dna_{t}"] = (stats.spearmanr(v[d], dna[d])[0]
                               if d.sum() >= 12 and (v[d] > 0).any() else np.nan)
    rec["n_trips_evaluable"] = len(rhos_m)
    rec["n_trips_rho_mapped_negative"] = int(sum(1 for r in rhos_m if r < 0))
    rec["median_rho_mapped_within_trip"] = float(np.median(rhos_m)) if rhos_m else np.nan
    dd = dna.notna() & (mapped > 0)
    rec["rho_dna_global"] = stats.spearmanr(v[dd], dna[dd])[0]
    rows.append(rec)

out = pd.DataFrame(rows)
out.to_csv(os.path.join(RES, "biomass_within_trip.tsv"), sep="\t",
           index=False, float_format="%.4f")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
print(out[["species", "role", "n_trips_evaluable", "n_trips_rho_mapped_negative",
           "median_rho_mapped_within_trip", "rho_dna_global",
           "rho_dna_trip1", "rho_dna_trip3"]].to_string(index=False))
print("\nper-trip detections (a trip with 0 is not evaluable; 1 to 4 is a rank on almost no data):")
print(out[["species"] + [f"n_det_{t}" for t in TRIPS]].to_string(index=False))
m = dna.notna() & (mapped > 0)
print("\nThe two biomass proxies are independent: Spearman(dna, mapped) = "
      f"{stats.spearmanr(dna[m], mapped[m])[0]:.3f} (n={int(m.sum())})")
