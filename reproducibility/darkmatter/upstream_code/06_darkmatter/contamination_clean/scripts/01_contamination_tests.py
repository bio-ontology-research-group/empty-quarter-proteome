#!/usr/bin/env python3
"""Test whether the four MAG species shared with the cultured catalogue are kit/lab contamination.

Read-only on the workbench. Writes results/*.tsv and figures/abundance_by_trip.png.

Tests
  1  prevalence + abundance overall, by trip, by compartment
  2  biomass relation: Spearman rho of abundance vs mapped-read fraction and vs input DNA ng/ul
  3  spatial structure: replicate concordance and site recurrence across trips
  4  prevalence-matched random controls give the null band each target is judged against
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean"
RES = os.path.join(P, "results")
FIG = os.path.join(P, "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
SEED = 20260910

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

# ---------------------------------------------------------------- load
ab = pd.read_csv(os.path.join(W, "10_thermal/abundance/rep_relative_abundance.tsv"),
                 sep="\t", index_col=0)
ab = ab.fillna(0.0)
meta = pd.read_csv(os.path.join(W, "metadata/sample_metadata.tsv"), sep="\t")
meta = meta.set_index("sample").reindex(ab.columns)
design = pd.read_csv(os.path.join(W, "10_thermal/abundance/sample_design.tsv"), sep="\t",
                     index_col="sample_id").reindex(ab.columns)
dna = pd.read_csv(os.path.join(W, "metadata/dna_conc.tsv"), sep="\t",
                  index_col="sample")["dna_ng_ul"].reindex(ab.columns)

mapped_pct = 100.0 - design["unmapped_percent"]          # = column sum of the matrix
usable = mapped_pct > 0                                   # F5Dr1 has zero mapped reads
TRIPS = ["trip1", "trip2", "trip3", "trip4", "trip5"]
COMPS = ["Surface", "Deep", "PlantRoot"]

print(f"matrix {ab.shape}; samples with mapped reads {usable.sum()}; "
      f"samples with DNA conc {dna.notna().sum()}")

# ---------------------------------------------------------------- matched randoms
prev_all = (ab > 0).sum(axis=1)
rng = np.random.default_rng(SEED)
matched = {}
pool = [g for g in ab.index if g not in TARGETS and g not in COMPARATORS]
for g in TARGETS:
    n = prev_all[g]
    lo, hi = 0.7 * n, 1.3 * n
    band = [x for x in pool if lo <= prev_all[x] <= hi]
    take = min(20, len(band))
    matched[g] = list(rng.choice(band, size=take, replace=False))
    print(f"matched randoms for {g}: n_det={n} band=[{lo:.0f},{hi:.0f}] pool={len(band)} used={take}")


def genome_stats(g):
    """Prevalence, abundance and structure statistics for one representative genome."""
    v = ab.loc[g]
    det = v > 0
    row = {
        "genome": g,
        "n_detected": int(det.sum()),
        "prevalence_pct": 100.0 * det.sum() / len(v),
        "median_rel_abund_all": float(v.median()),
        "median_rel_abund_detected": float(v[det].median()) if det.any() else np.nan,
        "max_rel_abund": float(v.max()),
        "max_sample": str(v.idxmax()),
    }
    row["max_over_median_detected"] = (row["max_rel_abund"] / row["median_rel_abund_detected"]
                                       if det.any() and row["median_rel_abund_detected"] > 0 else np.nan)
    for t in TRIPS:
        m = meta["trip"] == t
        row[f"prev_{t}"] = 100.0 * (det & m).sum() / max(m.sum(), 1)
        row[f"medabund_{t}"] = float(v[m & det].median()) if (m & det).any() else 0.0
    for c in COMPS:
        m = meta["compartment"] == c
        row[f"prev_{c}"] = 100.0 * (det & m).sum() / max(m.sum(), 1)

    # test 2: biomass relations (Spearman, over samples with a mapped fraction)
    u = usable & v.notna()
    row["rho_mapped_all"], row["p_mapped_all"] = stats.spearmanr(v[u], mapped_pct[u])
    d = u & dna.notna()
    if d.sum() > 10:
        row["rho_dna_all"], row["p_dna_all"] = stats.spearmanr(v[d], dna[d])
    else:
        row["rho_dna_all"], row["p_dna_all"] = np.nan, np.nan
    dd = u & det
    row["rho_mapped_detected"] = (stats.spearmanr(v[dd], mapped_pct[dd])[0] if dd.sum() > 10 else np.nan)
    ddd = dd & dna.notna()
    row["rho_dna_detected"] = (stats.spearmanr(v[ddd], dna[ddd])[0] if ddd.sum() > 10 else np.nan)

    # test 3a: replicate concordance inside site x compartment x trip
    grp = meta.groupby(["site", "compartment", "trip"], dropna=True).groups
    n_grp = n_all = n_none = 0
    for _, idx in grp.items():
        idx = [s for s in idx if s in det.index]
        if len(idx) < 2:
            continue
        n_grp += 1
        k = det.loc[idx].sum()
        if k == len(idx):
            n_all += 1
        elif k == 0:
            n_none += 1
    p = det.mean()
    row["n_repgroups"] = n_grp
    row["frac_groups_concordant"] = (n_all + n_none) / n_grp if n_grp else np.nan
    # expectation under random placement at the same overall prevalence
    exp = np.mean([p ** len(v_) + (1 - p) ** len(v_)
                   for v_ in [[s for s in i if s in det.index] for i in grp.values()]
                   if len(v_) >= 2]) if n_grp else np.nan
    row["exp_frac_concordant"] = exp
    row["concordance_excess"] = row["frac_groups_concordant"] - exp if n_grp else np.nan

    # test 3b: site recurrence across trips. P(present at site in trip B | present in trip A)
    st = pd.DataFrame({"site": meta["site"], "trip": meta["trip"], "det": det})
    sxt = st.groupby(["site", "trip"])["det"].max().unstack()
    pres = []
    for site in sxt.index:
        r = sxt.loc[site].dropna()
        if len(r) >= 2 and r.sum() >= 1:
            pres.append((r.sum(), len(r)))
    if pres:
        # among site-trip cells of sites that are positive somewhere, what fraction are positive
        row["site_recurrence"] = sum(a for a, _ in pres) / sum(b for _, b in pres)
        row["n_sites_positive"] = len(pres)
    else:
        row["site_recurrence"], row["n_sites_positive"] = np.nan, 0
    row["site_baseline"] = float(sxt.stack().mean())
    return row


rows = []
for g, sp in {**TARGETS, **COMPARATORS}.items():
    r = genome_stats(g)
    r["species"] = sp
    r["role"] = "target" if g in TARGETS else "comparator"
    rows.append(r)
for tg, gl in matched.items():
    for g in gl:
        r = genome_stats(g)
        r["species"] = "(matched random)"
        r["role"] = f"random_for_{TARGETS[tg].split()[0]}"
        rows.append(r)
out = pd.DataFrame(rows)
cols = ["genome", "species", "role"] + [c for c in out.columns if c not in ("genome", "species", "role")]
out = out[cols]
out.to_csv(os.path.join(RES, "per_genome_stats.tsv"), sep="\t", index=False, float_format="%.5g")

# ---------------------------------------------------------------- target vs its null band
summ = []
for g, sp in TARGETS.items():
    tr = out.loc[out.genome == g].iloc[0]
    nul = out.loc[out.role == f"random_for_{sp.split()[0]}"]
    rec = {"genome": g, "species": sp, "n_detected": int(tr.n_detected),
           "prevalence_pct": tr.prevalence_pct,
           "median_rel_abund_detected": tr.median_rel_abund_detected,
           "max_rel_abund": tr.max_rel_abund, "max_sample": tr.max_sample,
           "max_over_median_detected": tr.max_over_median_detected}
    for met in ["rho_mapped_all", "rho_dna_all", "concordance_excess", "site_recurrence",
                "median_rel_abund_detected"]:
        vals = nul[met].astype(float).dropna()
        rec[f"{met}__target"] = float(tr[met]) if pd.notna(tr[met]) else np.nan
        rec[f"{met}__null_median"] = float(vals.median()) if len(vals) else np.nan
        rec[f"{met}__null_p5"] = float(vals.quantile(0.05)) if len(vals) else np.nan
        rec[f"{met}__null_p95"] = float(vals.quantile(0.95)) if len(vals) else np.nan
        rec[f"{met}__pctile"] = (100.0 * (vals < tr[met]).mean()
                                 if len(vals) and pd.notna(tr[met]) else np.nan)
    summ.append(rec)
pd.DataFrame(summ).to_csv(os.path.join(RES, "target_vs_null.tsv"), sep="\t",
                          index=False, float_format="%.5g")

# ---------------------------------------------------------------- trip table
trip_rows = []
for g, sp in {**TARGETS, **COMPARATORS}.items():
    v = ab.loc[g]
    det = v > 0
    for t in TRIPS:
        m = meta["trip"] == t
        trip_rows.append({"genome": g, "species": sp,
                          "role": "target" if g in TARGETS else "comparator",
                          "trip": t, "n_samples": int(m.sum()),
                          "n_detected": int((det & m).sum()),
                          "prevalence_pct": 100.0 * (det & m).sum() / max(m.sum(), 1),
                          "median_detected": float(v[m & det].median()) if (m & det).any() else 0.0,
                          "mean_all": float(v[m].mean())})
pd.DataFrame(trip_rows).to_csv(os.path.join(RES, "by_trip.tsv"), sep="\t",
                               index=False, float_format="%.5g")

# ---------------------------------------------------------------- figure
# sharex must stay False: shared axes collapse the per-panel tick labels onto one set, which
# would print the last panel's prevalences under every panel.
fig, axes = plt.subplots(2, 4, figsize=(17, 8), sharex=False)
order = list(TARGETS) + list(COMPARATORS)
labels = {**TARGETS, **COMPARATORS}
for ax, g in zip(axes.ravel(), order):
    v = ab.loc[g]
    data, prev = [], []
    for t in TRIPS:
        m = (meta["trip"] == t).values
        x = v[m]
        data.append(x[x > 0].values)
        prev.append(100.0 * (x > 0).mean())
    pos = np.arange(1, 6)
    for i, d in enumerate(data):
        if len(d):
            ax.scatter(np.full(len(d), pos[i]) + rng.normal(0, .07, len(d)), d,
                       s=9, alpha=.55,
                       color="#c0392b" if g in TARGETS else "#2471a3")
    bx = [d if len(d) else np.array([np.nan]) for d in data]
    ax.boxplot(bx, positions=pos, widths=.55, showfliers=False,
               medianprops=dict(color="k"))
    ax.set_yscale("log")
    ax.set_title(f"{labels[g]}\n{g}   n={int((v > 0).sum())}/{len(v)}", fontsize=8)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"t{i}\n{p:.0f}%" for i, p in zip(range(1, 6), prev)], fontsize=7)
    ax.set_xlim(0.4, 5.6)
    ax.grid(alpha=.25, axis="y")
for ax in axes[:, 0]:
    ax.set_ylabel("relative abundance (%), detected samples")
fig.suptitle("Relative abundance by trip. Top row: shared-with-culture targets (red). "
             "Bottom row: comparators (blue). x label = trip and prevalence in that trip.",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "abundance_by_trip.png"), dpi=170)
print("wrote", RES, FIG)
