#!/usr/bin/env python
"""ROUND TWO, GOAL A (models) + GOAL B (the tautology test).

A  native-key dark fraction from results/v2_rep_dark.tsv, compared with round one,
   then the round-one model ladder rerun on the native values.
B  taxonomic novelty (first unnamed GTDB rank) added to the ladder, and the
   named-genus versus placeholder-genus contrast inside single phyla.

Writes results/v2_model_report.txt, results/v2_reps.tsv, results/v2_phylum_dark.tsv,
results/v2_genus_dark.tsv, results/v2_unjoined_reps.tsv and
figures/v2_residual_dark_by_phylum.png.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
sys.path.insert(0, HERE)
import common as eq  # patsy needs the name C for its categorical operator

W = eq.W
lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ---------------------------------------------------------------- novelty
STEM = r"[A-Z][a-z]+(?:_[A-Z]{1,3})?"
LATIN = re.compile(r"^%s$" % STEM)
LATIN_SP = re.compile(r"^%s [a-z]+(?:_[A-Z]{1,3})?$" % STEM)   # s__ is "Genus species"
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
LETTER = {"phylum": "p", "class": "c", "order": "o", "family": "f",
          "genus": "g", "species": "s"}


def is_named(v, rank="genus"):
    """A GTDB name counts as Latin if it has no digit and looks like a Latin stem.

    Placeholders such as JAJPHT01, UBA1234, DASUME01, PMG-095, AWTP1-13,
    SW-10-69-26 and EX4484-52 all carry digits; blanks are unnamed.
    """
    if not isinstance(v, str) or v == "":
        return False
    if any(ch.isdigit() for ch in v):
        return False
    return bool((LATIN_SP if rank == "species" else LATIN).match(v))


def split_ranks(s):
    out = {}
    if isinstance(s, str):
        for part in s.split(";"):
            part = part.strip()
            if len(part) > 3 and part[1:3] == "__":
                out[part[0]] = part[3:]
    return out


def load_taxonomy():
    tax = pd.concat([
        pd.read_csv(W + "/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
                    sep="\t", usecols=["user_genome", "classification"]),
        pd.read_csv(W + "/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv",
                    sep="\t", usecols=["user_genome", "classification"]),
    ]).rename(columns={"user_genome": "genome"})
    rows = []
    odd = set()
    for g, cl in zip(tax["genome"], tax["classification"]):
        d = split_ranks(cl)
        rec = {"genome": g, "domain": d.get("d", "")}
        first = 7                       # 7 = every rank phylum..species is Latin
        for k, r in enumerate(RANKS, start=1):
            v = d.get(LETTER[r], "")
            rec[r] = v
            named = is_named(v, r)
            rec[r + "_named"] = named
            if v and not named and not any(ch.isdigit() for ch in v):
                odd.add(v)
            if not named and first == 7:
                first = k
        rec["first_unnamed_rank"] = first
        rows.append(rec)
    return pd.DataFrame(rows), sorted(odd)


# ---------------------------------------------------------------- load
q = eq.load_reps()[["genome", "completeness", "contamination", "length", "N50",
                    "log10_length", "log10_N50", "phylum", "genus"]]
v2 = pd.read_csv(os.path.join(RES, "v2_rep_dark.tsv"), sep="\t").rename(columns={"rep": "genome"})
tax, odd_names = load_taxonomy()

d = q.merge(v2, on="genome", how="left").merge(
    tax.drop(columns=["phylum", "genus"]), on="genome", how="left")
d["named_genus"] = d["genus_named"].fillna(False)
d["first_unnamed_rank"] = d["first_unnamed_rank"].fillna(1).astype(int)

n_genus = d["genus"].value_counts()
small = set(n_genus[n_genus < 5].index)
d["genus_grp"] = np.where(d["genus"].isin(small), "other", d["genus"])
d["genus_class"] = np.where(d["genus"] == "unclassified", "blank",
                            np.where(d["named_genus"], "named", "placeholder"))

say("ROUND TWO, 01 dark fraction by lineage")
say("======================================")
say()
say("GOAL A. NATIVE-KEY DARK FRACTION")
say("--------------------------------")
say("reps in Wdb                     : %d" % len(d))
say("reps with a native dark fraction: %d" % int((d["n_proteins"] > 0).sum()))
say("n_proteins  median %.0f   total %d" % (d.n_proteins.median(), int(d.n_proteins.sum())))
say("dark_fraction median %.4f  mean %.4f  IQR %.4f - %.4f"
    % (d.dark_fraction.median(), d.dark_fraction.mean(),
       d.dark_fraction.quantile(.25), d.dark_fraction.quantile(.75)))
say()

# ---------------------------------------------------------------- round one comparison
r1 = pd.read_csv(os.path.join(RES, "reps_dark.tsv"), sep="\t")[
    ["genome", "form", "total_proteins", "dark_proteins", "dark_fraction"]]
r1 = r1.rename(columns={"total_proteins": "r1_proteins", "dark_proteins": "r1_dark",
                        "dark_fraction": "r1_fraction"})
allforms = pd.read_csv(os.path.join(RES, "reps_all.tsv"), sep="\t")[["genome", "form"]]
d = d.merge(allforms, on="genome", how="left")
m = d.merge(r1.drop(columns=["form"]), on="genome", how="left")

j = m[m["r1_fraction"].notna()]
rho, p = spearmanr(j["dark_fraction"], j["r1_fraction"])
say("COMPARISON WITH ROUND ONE (the %d reps round one could join)" % len(j))
say("  Spearman rho %.4f (p %.3g)   median |v2 - v1| %.4f   mean |v2 - v1| %.4f"
    % (rho, p, (j.dark_fraction - j.r1_fraction).abs().median(),
       (j.dark_fraction - j.r1_fraction).abs().mean()))
say("  identical protein count: %d / %d (%.1f%%)"
    % (int((j.n_proteins == j.r1_proteins).sum()), len(j),
       100.0 * (j.n_proteins == j.r1_proteins).mean()))
say("  %-8s %6s %10s %14s %14s %12s" % ("form", "n", "rho", "median|diff|", "exact n_prot", "median v2-v1"))
for f, sub in j.groupby("form"):
    r_, _ = spearmanr(sub["dark_fraction"], sub["r1_fraction"])
    say("  %-8s %6d %10.4f %14.4f %13.1f%% %+12.4f"
        % (f, len(sub), r_, (sub.dark_fraction - sub.r1_fraction).abs().median(),
           100.0 * (sub.n_proteins == sub.r1_proteins).mean(),
           (sub.dark_fraction - sub.r1_fraction).median()))
say()

miss = m[m["r1_fraction"].isna()].copy()
miss["reason"] = np.where(miss["form"] == "no_assembly_token",
                          "no assembly token (48PRr1 family)", "unmatched _sub")
cols = ["genome", "reason", "phylum", "genus", "completeness", "n_proteins", "n_dark", "dark_fraction"]
miss[cols].sort_values(["reason", "genome"]).to_csv(
    os.path.join(RES, "v2_unjoined_reps.tsv"), sep="\t", index=False)
say("THE %d REPS ROUND ONE DROPPED, NOW COMPUTED NATIVELY" % len(miss))
for reason, sub in miss.groupby("reason"):
    say("  %s: n = %d   median dark_fraction %.4f (round-one modelled median %.4f)"
        % (reason, len(sub), sub.dark_fraction.median(), j.dark_fraction.median()))
    say("  %-34s %-20s %10s %8s %8s" % ("rep", "phylum", "n_proteins", "n_dark", "dark_fr"))
    for _, r in sub.sort_values("dark_fraction", ascending=False).iterrows():
        say("  %-34s %-20s %10d %8d %8.4f"
            % (r.genome, r.phylum[:20], r.n_proteins, r.n_dark, r.dark_fraction))
    say()

d = m

# ---------------------------------------------------------------- GOAL A models
d = d.dropna(subset=["dark_fraction", "completeness", "contamination", "log10_N50", "log10_length"])
Q = "completeness + contamination + log10_N50 + log10_length"
mA = smf.ols("dark_fraction ~ " + Q, data=d).fit()
mB = smf.ols("dark_fraction ~ " + Q + " + C(phylum)", data=d).fit()
mC = smf.ols("dark_fraction ~ " + Q + " + C(phylum) + C(genus_grp)", data=d).fit()
mL1 = smf.ols("dark_fraction ~ C(phylum)", data=d).fit()
mL2 = smf.ols("dark_fraction ~ C(phylum) + C(genus_grp)", data=d).fit()
mL3 = smf.ols("dark_fraction ~ C(phylum) + C(genus_grp) + " + Q, data=d).fit()

say("MODEL LADDER ON THE NATIVE VALUES (n = %d)" % int(mA.nobs))
say("  %-34s %8s %8s %9s" % ("model", "R2", "adjR2", "dR2"))
say("  %-34s %8.4f %8.4f %9s" % ("A quality", mA.rsquared, mA.rsquared_adj, "-"))
say("  %-34s %8.4f %8.4f %+9.4f" % ("B = A + phylum", mB.rsquared, mB.rsquared_adj, mB.rsquared - mA.rsquared))
say("  %-34s %8.4f %8.4f %+9.4f" % ("C = B + genus", mC.rsquared, mC.rsquared_adj, mC.rsquared - mB.rsquared))
say("  %-34s %8.4f %8.4f %9s" % ("lineage first: phylum", mL1.rsquared, mL1.rsquared_adj, "-"))
say("  %-34s %8.4f %8.4f %+9.4f" % ("phylum + genus", mL2.rsquared, mL2.rsquared_adj, mL2.rsquared - mL1.rsquared))
say("  %-34s %8.4f %8.4f %+9.4f" % ("phylum + genus + quality", mL3.rsquared, mL3.rsquared_adj, mL3.rsquared - mL2.rsquared))
say()
say("VARIANCE PARTITION (full model R2 %.4f)" % mL3.rsquared)
say("  unique to quality %.4f   unique to lineage %.4f   shared %.4f   unexplained %.4f"
    % (mL3.rsquared - mL2.rsquared, mL3.rsquared - mA.rsquared,
       mA.rsquared + mL2.rsquared - mL3.rsquared, 1 - mL3.rsquared))
say("MODEL A slopes: " + "  ".join("%s %+.5f" % (k, mA.params[k])
                                   for k in mA.params.index if k != "Intercept"))
say()

d["resid_A"] = mA.resid

# ---------------------------------------------------------------- GOAL B
say("GOAL B. THE TAUTOLOGY TEST")
say("--------------------------")
say("first_unnamed_rank: deepest Latin rank + 1 on the ladder phylum=1 class=2 order=3")
say("family=4 genus=5 species=6; 7 = every rank phylum..species carries a Latin name.")
say("Placeholder = contains a digit, or does not match ^[A-Z][a-z]+(_[A-Z]{1,3})?$; blank = unnamed.")
if odd_names:
    say("digit-free names classed as placeholders (eyeball list, %d): %s"
        % (len(odd_names), ", ".join(odd_names[:40])))
else:
    say("every placeholder call was driven by a digit in the name.")
say()
vc = d["first_unnamed_rank"].value_counts().sort_index()
say("  %-22s %7s %10s" % ("first unnamed rank", "n", "median dark"))
for k in vc.index:
    lab = (RANKS[k - 1] if k <= 6 else "none (all Latin)")
    say("  %-22s %7d %10.4f" % (lab, vc[k], d.loc[d.first_unnamed_rank == k, "dark_fraction"].median()))
say("  named genus: %d of %d (%.1f%%);  placeholder genus %d;  blank g__ %d"
    % (int((d.genus_class == "named").sum()), len(d), 100.0 * (d.genus_class == "named").mean(),
       int((d.genus_class == "placeholder").sum()), int((d.genus_class == "blank").sum())))
say()

nA = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank", data=d).fit()
nAg = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank + named_genus", data=d).fit()
nB = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank + C(phylum)", data=d).fit()
nC = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank + C(phylum) + C(genus_grp)", data=d).fit()
say("NOVELTY LADDER")
say("  %-40s %8s %8s %9s" % ("model", "R2", "adjR2", "dR2"))
say("  %-40s %8.4f %8.4f %9s" % ("A quality", mA.rsquared, mA.rsquared_adj, "-"))
say("  %-40s %8.4f %8.4f %+9.4f" % ("N = A + first_unnamed_rank", nA.rsquared, nA.rsquared_adj, nA.rsquared - mA.rsquared))
say("  %-40s %8.4f %8.4f %+9.4f" % ("N + named_genus flag", nAg.rsquared, nAg.rsquared_adj, nAg.rsquared - nA.rsquared))
say("  %-40s %8.4f %8.4f %+9.4f" % ("N + phylum", nB.rsquared, nB.rsquared_adj, nB.rsquared - nA.rsquared))
say("  %-40s %8.4f %8.4f %+9.4f" % ("N + phylum + genus", nC.rsquared, nC.rsquared_adj, nC.rsquared - nB.rsquared))
say("  first_unnamed_rank slope in N: %+.5f (p %.3g); named_genus slope %+.5f (p %.3g)"
    % (nA.params["first_unnamed_rank"], nA.pvalues["first_unnamed_rank"],
       nAg.params.get("named_genus[T.True]", np.nan), nAg.pvalues.get("named_genus[T.True]", np.nan)))
say("  phylum increment WITHOUT novelty %+.4f -> WITH novelty %+.4f (survives %.1f%%)"
    % (mB.rsquared - mA.rsquared, nB.rsquared - nA.rsquared,
       100.0 * (nB.rsquared - nA.rsquared) / max(mB.rsquared - mA.rsquared, 1e-9)))
say()

# within-phylum named vs placeholder
bac = d[d["domain"] == "Bacteria"]
big_bac = list(bac["phylum"].value_counts().head(4).index)
targets = ["Thermoplasmatota"] + big_bac
say("WITHIN-PHYLUM: named genus versus placeholder genus")
say("  %-18s %6s %20s %20s %20s %12s %12s"
    % ("phylum", "n", "named n / median", "placeholder n / med", "blank n / median",
       "MWU named-plc", "MWU named-rest"))
within = []
for ph in targets:
    s = d[d["phylum"] == ph]
    grp = {k: s.loc[s.genus_class == k, "dark_fraction"].values
           for k in ("named", "placeholder", "blank")}
    rest = np.concatenate([grp["placeholder"], grp["blank"]])
    def cell(v):
        return "%d / %s" % (len(v), ("%.4f" % np.median(v)) if len(v) else "-")
    def mwu(a, b):
        if len(a) < 5 or len(b) < 5:
            return "n<5"
        return "%.3g" % mannwhitneyu(a, b, alternative="two-sided").pvalue
    say("  %-18s %6d %20s %20s %20s %12s %12s"
        % (ph, len(s), cell(grp["named"]), cell(grp["placeholder"]), cell(grp["blank"]),
           mwu(grp["named"], grp["placeholder"]), mwu(grp["named"], rest)))
    within.append((ph, len(s), len(grp["named"]), len(grp["placeholder"]), len(grp["blank"]),
                   np.median(grp["named"]) if len(grp["named"]) else np.nan,
                   np.median(grp["placeholder"]) if len(grp["placeholder"]) else np.nan,
                   mwu(grp["named"], grp["placeholder"])))
say()

# novelty-matched control for Thermoplasmatota
th = d[d["phylum"] == "Thermoplasmatota"]
say("NOVELTY-MATCHED CONTROL for Thermoplasmatota (n = %d)" % len(th))
say("  %-24s %8s %12s %8s %12s %12s"
    % ("first_unnamed_rank", "Thermo n", "Thermo med", "other n", "other med", "MWU p"))
for k in sorted(th["first_unnamed_rank"].unique()):
    a = th.loc[th.first_unnamed_rank == k, "dark_fraction"].values
    b = d.loc[(d.first_unnamed_rank == k) & (d.phylum != "Thermoplasmatota"), "dark_fraction"].values
    pv = ("%.3g" % mannwhitneyu(a, b, alternative="two-sided").pvalue) if (len(a) >= 5 and len(b) >= 5) else "n<5"
    say("  %-24s %8d %12.4f %8d %12.4f %12s"
        % (RANKS[k - 1] if k <= 6 else "none", len(a), np.median(a), len(b), np.median(b), pv))
say()

# ---------------------------------------------------------------- per-lineage tables
g = d.groupby("genus").agg(
    phylum=("phylum", lambda s: s.mode().iat[0]),
    genus_class=("genus_class", lambda s: s.mode().iat[0]),
    n_reps=("genome", "size"),
    median_dark_fraction=("dark_fraction", "median"),
    median_resid_A=("resid_A", "median"),
    median_completeness=("completeness", "median"),
).reset_index().sort_values("median_resid_A", ascending=False)
g.round(6).to_csv(os.path.join(RES, "v2_genus_dark.tsv"), sep="\t", index=False)

p = d.groupby("phylum").agg(
    n_reps=("genome", "size"),
    median_dark_fraction=("dark_fraction", "median"),
    median_resid_A=("resid_A", "median"),
    median_completeness=("completeness", "median"),
    frac_named_genus=("named_genus", "mean"),
    median_first_unnamed=("first_unnamed_rank", "median"),
).reset_index()
p.round(6).to_csv(os.path.join(RES, "v2_phylum_dark.tsv"), sep="\t", index=False)
pk = p[p.n_reps >= 20].sort_values("median_resid_A", ascending=False)
say("PHYLA (n >= 20) by median Model A residual")
say(pk[["phylum", "n_reps", "median_dark_fraction", "median_resid_A",
        "frac_named_genus", "median_first_unnamed"]].to_string(index=False))
say()
g10 = g[(g.n_reps >= 10) & (g.genus != "unclassified")]
say("TOP 12 genera by median Model A residual (n >= 10)")
say(g10.head(12)[["genus", "phylum", "genus_class", "n_reps", "median_dark_fraction", "median_resid_A"]].to_string(index=False))
say()
say("BOTTOM 12 genera by median Model A residual (n >= 10)")
say(g10.tail(12)[["genus", "phylum", "genus_class", "n_reps", "median_dark_fraction", "median_resid_A"]].to_string(index=False))
say()

# ---------------------------------------------------------------- figure
keep = list(pk.phylum)
sub = d[d.phylum.isin(keep)]
order = sub.groupby("phylum").resid_A.median().sort_values().index
fig, ax = plt.subplots(figsize=(12, max(5, 0.28 * len(order))))
ax.boxplot([sub.loc[sub.phylum == k, "resid_A"].values for k in order],
           vert=False, labels=list(order), showfliers=False, widths=0.65)
ax.axvline(0, color="crimson", lw=1)
ax.set_xlabel("Model A residual dark fraction (observed minus quality prediction)")
ax.set_title("Native-key dark fraction residual by phylum (phyla with >=20 reps, n=%d)" % len(sub))
ax.tick_params(axis="y", labelsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "v2_residual_dark_by_phylum.png"), dpi=150)
plt.close(fig)

out = d[["genome", "form", "completeness", "contamination", "length", "N50",
         "log10_length", "log10_N50", "phylum", "genus", "genus_grp", "genus_class",
         "named_genus", "first_unnamed_rank", "domain",
         "n_proteins", "n_dark", "dark_fraction", "resid_A",
         "r1_proteins", "r1_dark", "r1_fraction"]]
out.to_csv(os.path.join(RES, "v2_reps.tsv"), sep="\t", index=False)

with open(os.path.join(RES, "v2_model_report.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\nwrote results/v2_reps.tsv, v2_model_report.txt, v2_phylum_dark.tsv, "
      "v2_genus_dark.tsv, v2_unjoined_reps.tsv and figures/v2_residual_dark_by_phylum.png")
