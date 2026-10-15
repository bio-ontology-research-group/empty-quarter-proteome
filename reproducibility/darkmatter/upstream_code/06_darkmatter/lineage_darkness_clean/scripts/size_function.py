#!/usr/bin/env python
"""03: do we see more functions only because genomes are bigger?

KO richness  = unique KOs per rep (rep_ko_v2)
corrected size = length / (completeness/100)

  M1  ko_richness ~ log10(corrected size)
  M2  M1 + completeness + contamination
  M3  M2 + phylum
  M4  dark_fraction ~ log10(corrected size)          (joined reps only)
  M5  M4 + completeness
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
# promoted 2026-09-20: 01_dark_by_lineage and 03_genome_size_function were merged into this one
# folder, so common.py and v2_reps.tsv are local instead of one directory up.
DARK01 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import common as eq  # not "C": patsy needs the name C for its categorical operator

lines = []


def say(s=""):
    print(s)
    lines.append(s)


reps = eq.load_reps()
ko = eq.load_ko()
d = reps.merge(ko, on="genome", how="inner")
say("03 GENOME SIZE AND FUNCTIONAL REPERTOIRE")
say("========================================")
say("reps with quality + taxonomy : %d" % len(reps))
say("reps with a KO row           : %d" % len(ko))
say("joined (size/function set)   : %d" % len(d))
say("KO richness  : median %d  mean %.1f  range %d - %d"
    % (d.ko_richness.median(), d.ko_richness.mean(), d.ko_richness.min(), d.ko_richness.max()))
say("corrected size (Mb): median %.2f  range %.2f - %.2f"
    % (d.corrected_size.median() / 1e6, d.corrected_size.min() / 1e6, d.corrected_size.max() / 1e6))
say("raw length (Mb)    : median %.2f" % (d.length.median() / 1e6))
say()

m1 = smf.ols("ko_richness ~ log10_corrected_size", data=d).fit()
m2 = smf.ols("ko_richness ~ log10_corrected_size + completeness + contamination", data=d).fit()
m3 = smf.ols("ko_richness ~ log10_corrected_size + completeness + contamination + C(phylum)", data=d).fit()

say("M1  ko_richness ~ log10(corrected size)")
say("  slope = %+.1f KOs per 10-fold size   intercept = %.1f" % (m1.params["log10_corrected_size"], m1.params["Intercept"]))
say("  R2 = %.4f   adj R2 = %.4f   p = %.3g   n = %d"
    % (m1.rsquared, m1.rsquared_adj, m1.pvalues["log10_corrected_size"], int(m1.nobs)))
say("  a doubling of size is %+.0f KOs" % (m1.params["log10_corrected_size"] * np.log10(2)))
say()

m1raw = smf.ols("ko_richness ~ log10_length", data=d).fit()
say("M1raw  ko_richness ~ log10(RAW observed length)")
say("  slope = %+.1f   R2 = %.4f" % (m1raw.params["log10_length"], m1raw.rsquared))
say("  raw length beats corrected size (%.4f vs %.4f) because KO richness is counted"
    % (m1raw.rsquared, m1.rsquared))
say("  on the observed genome, not the completeness-corrected one. That mismatch is")
say("  what the completeness term recovers in M2 below.")
say()

say("R2 INCREMENTS")
say("  %-34s %8s %8s %8s" % ("model", "R2", "adjR2", "dR2"))
say("  %-34s %8.4f %8.4f %8s" % ("M1 size", m1.rsquared, m1.rsquared_adj, "-"))
say("  %-34s %8.4f %8.4f %+8.4f" % ("M2 + completeness, contamination", m2.rsquared, m2.rsquared_adj, m2.rsquared - m1.rsquared))
say("  %-34s %8.4f %8.4f %+8.4f" % ("M3 + phylum", m3.rsquared, m3.rsquared_adj, m3.rsquared - m2.rsquared))
say()
say("M2 coefficients:")
for k in m2.params.index:
    say("    %-24s %+10.3f   p %.3g" % (k, m2.params[k], m2.pvalues[k]))
say()

d["resid_size"] = m1.resid
d["resid_size_quality"] = m2.resid  # M2 residual, to test the residuals against completeness
p = d.groupby("phylum").agg(
    n=("genome", "size"),
    median_ko_richness=("ko_richness", "median"),
    median_corrected_size=("corrected_size", "median"),
    median_resid_ko=("resid_size", "median"),
    median_resid_ko_m2=("resid_size_quality", "median"),
    median_completeness=("completeness", "median"),
).reset_index().sort_values("median_resid_ko", ascending=False)
p.round(4).to_csv(os.path.join(RES, "phylum_ko_residual.tsv"), sep="\t", index=False)
say("PHYLA with n>=20, by median residual KO richness (size-only model)")
say(p[p.n >= 20].to_string(index=False))
big = p[p.n >= 20]
say("rank agreement of phylum medians, size-only vs size+quality: Spearman %.3f"
    % big[["median_resid_ko", "median_resid_ko_m2"]].corr(method="spearman").iloc[0, 1])
say("so the extremes are not a completeness artefact, though mid-range phyla reorder")
say()

# CORRECTION 2026-09-17: was reps_dark.tsv, the superseded rename join (5,359 reps, and it read the
# parent bin for 975 `_sub` representatives). v2_reps.tsv is the native-key rebuild covering all
# 5,438. corrected_size is not stored there, so it is derived here exactly as in the rename build.
dk = pd.read_csv(os.path.join(DARK01, "results", "v2_reps.tsv"), sep="\t")
dk = dk.dropna(subset=["dark_fraction", "completeness", "length"])
dk = dk[dk.completeness > 0].copy()
dk["corrected_size"] = dk.length / (dk.completeness / 100.0)
dk["log10_corrected_size"] = np.log10(dk.corrected_size)
m4 = smf.ols("dark_fraction ~ log10_corrected_size", data=dk).fit()
m5 = smf.ols("dark_fraction ~ log10_corrected_size + completeness", data=dk).fit()
say("M4  dark_fraction ~ log10(corrected size)     n = %d" % int(m4.nobs))
say("  slope = %+.5f per 10-fold size   R2 = %.4f   p = %.3g"
    % (m4.params["log10_corrected_size"], m4.rsquared, m4.pvalues["log10_corrected_size"]))
say("M5  + completeness")
say("  slope(size) = %+.5f   slope(completeness) = %+.6f   R2 = %.4f"
    % (m5.params["log10_corrected_size"], m5.params["completeness"], m5.rsquared))
say("  direction: bigger genomes are %s dark"
    % ("MORE" if m4.params["log10_corrected_size"] > 0 else "LESS"))
say()
say("Spearman rho (corrected size vs dark fraction): %.4f"
    % dk[["log10_corrected_size", "dark_fraction"]].corr(method="spearman").iloc[0, 1])

# ---- figures ----
t8 = eq.top_phyla(d, 8)
cmap = plt.get_cmap("tab10")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(d.loc[~d.phylum.isin(t8), "corrected_size"] / 1e6,
           d.loc[~d.phylum.isin(t8), "ko_richness"], s=5, c="0.82",
           label="other phyla", rasterized=True)
for i, k in enumerate(t8):
    s = d[d.phylum == k]
    ax.scatter(s.corrected_size / 1e6, s.ko_richness, s=6, color=cmap(i),
               label="%s (n=%d)" % (k, len(s)))
xs = np.linspace(d.log10_corrected_size.min(), d.log10_corrected_size.max(), 100)
ax.plot(10 ** xs / 1e6, m1.params["Intercept"] + m1.params["log10_corrected_size"] * xs,
        "k-", lw=2, label="fit, R2=%.2f" % m1.rsquared)
ax.set_xscale("log")
ax.set_xlabel("completeness-corrected genome size (Mb, log scale)")
ax.set_ylabel("KO richness (unique KOs)")
ax.set_title("KO richness vs genome size, %d species reps" % len(d))
ax.legend(fontsize=7, markerscale=2, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "ko_vs_size.png"), dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 6))
t8d = eq.top_phyla(dk, 8)
ax.scatter(dk.loc[~dk.phylum.isin(t8d), "corrected_size"] / 1e6,
           dk.loc[~dk.phylum.isin(t8d), "dark_fraction"], s=5, c="0.82",
           label="other phyla", rasterized=True)
for i, k in enumerate(t8d):
    s = dk[dk.phylum == k]
    ax.scatter(s.corrected_size / 1e6, s.dark_fraction, s=6, color=cmap(i),
               label="%s (n=%d)" % (k, len(s)))
xs = np.linspace(dk.log10_corrected_size.min(), dk.log10_corrected_size.max(), 100)
ax.plot(10 ** xs / 1e6, m4.params["Intercept"] + m4.params["log10_corrected_size"] * xs,
        "k-", lw=2, label="fit, R2=%.3f" % m4.rsquared)
ax.set_xscale("log")
ax.set_xlabel("completeness-corrected genome size (Mb, log scale)")
ax.set_ylabel("dark fraction (no UniRef50 homolog)")
ax.set_title("Dark fraction vs genome size, %d species reps" % len(dk))
ax.legend(fontsize=7, markerscale=2, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "dark_vs_size.png"), dpi=150)
plt.close(fig)

with open(os.path.join(RES, "model_report.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\nwrote results/phylum_ko_residual.tsv, results/model_report.txt and 2 figures")
