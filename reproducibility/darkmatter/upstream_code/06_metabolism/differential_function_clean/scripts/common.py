#!/usr/bin/env python3
"""Shared loaders and statistics for the round-two differential-function scan.

Everything here is read-only with respect to the workbench (W) and the paper repo (R).
Outputs are written only under this folder, which locates itself from the script path, so
the tree can be moved without editing a constant.

Units convention used throughout
--------------------------------
`rep_relative_abundance.tsv` columns are percentages of TOTAL READS in the sample that
map to each species representative. Column sums run ~40-57 % (mean 47.7 %), i.e. about
half the reads map to no genome in the catalogue. Every carriage number below is
therefore normalised by the column sum, making it

    carriage(feature, sample) = 100 * sum_i ra_i * pres_i(feature) / sum_i ra_i

= the percentage of the MAPPED community that encodes the feature. Differences are then
in percentage points of carriage, which is what the 5 pp effect floor refers to.
"""
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as tdist, norm

W = "/ibex/scratch/projects/c2014/rob/metagenomics/"
R = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags/"
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, "..")) + "/"
RES = os.path.join(P, "results")
FIG = os.path.join(P, "figures")
LOGS = os.path.join(P, "logs")

ABUND = W + "10_thermal/abundance/rep_relative_abundance.tsv"
REPKO = W + "04_functional/repko_clean/results/rep_ko_v2.tsv"
META = W + "metadata/sample_metadata.tsv"
ENVM = R + "results/darkmatter_ecology/env_matrix.tsv"
GENOMEINFO = W + "05_refine/_drep/out/data_tables/genomeInformation.csv"
GTDB = [W + "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
        W + "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"]

MODULE_FRACTION = 2.0 / 3.0   # a module is present when >= 2/3 of its KOs are present
KO_MIN_PREV = 0.05            # KO-level universe: present in >= 5 % of species
KO_MAX_PREV = 0.95            # ... and <= 95 %
EFFECT_FLOOR_PP = 5.0         # |median paired difference| floor, percentage points
RHO_FLOOR = 0.30              # |partial Spearman rho| floor
QCUT = 0.05                   # BH-FDR level
PHYLUM_MIN_FRAC = 0.01        # phylum must hold >= 1 % of mapped abundance in a sample
CONTROL_PHYLA = ["Actinomycetota", "Pseudomonadota"]


# ---------------------------------------------------------------- loaders
def load_abundance():
    """species x sample percentages. Drops samples with zero mapped abundance."""
    ab = pd.read_csv(ABUND, sep="\t", index_col=0)
    keep = ab.sum(axis=0) > 0
    dropped = list(ab.columns[~keep])
    ab = ab.loc[:, keep]
    return ab, dropped


def load_ko_sets():
    """genome -> set of KO ids. rep_ko_v2.tsv has no header: genome TAB comma-joined KOs."""
    out = {}
    with open(REPKO) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts[0]:
                continue
            kos = parts[1].split(",") if len(parts) > 1 and parts[1] else []
            out[parts[0]] = set(k for k in kos if k)
    return out


def load_taxonomy():
    """genome -> phylum (GTDB r232, bacteria + archaea)."""
    tax = {}
    for path in GTDB:
        d = pd.read_csv(path, sep="\t", usecols=["user_genome", "classification"])
        for g, c in zip(d["user_genome"], d["classification"]):
            fields = dict(f.split("__", 1) for f in c.split(";") if "__" in f)
            tax[g] = fields.get("p", "")
    return tax


def load_genome_size():
    """genome -> estimated complete genome size in Mb.

    dRep `length` is the assembled length of an incomplete bin, so it is divided by the
    CheckM2 completeness fraction to estimate the size of the underlying genome.
    """
    gi = pd.read_csv(GENOMEINFO)
    gi["g"] = gi["genome"].str.replace(r"\.fa$", "", regex=True)
    gi = gi.dropna(subset=["length", "completeness"])
    gi = gi[gi["completeness"] > 0]
    size = gi["length"] / (gi["completeness"] / 100.0) / 1e6
    return dict(zip(gi["g"], size))


def load_design():
    """sample -> site, compartment, trip (from the workbench metadata table)."""
    m = pd.read_csv(META, sep="\t", usecols=["sample", "site", "compartment", "trip"])
    m["site"] = m["site"].astype(str)
    return m.set_index("sample")


def load_env():
    return pd.read_csv(ENVM, sep="\t")


# ---------------------------------------------------------------- feature matrix
def build_presence(genomes, feature_defs, ko_sets):
    """genomes x features 0/1 presence.

    feature_defs: dict feature_id -> (list of STEP GROUPS, fraction threshold), where a
    group is a set of mutually alternative KOs. A step counts as satisfied when the genome
    carries any KO in the group, and the feature is present when the satisfied fraction
    reaches the threshold. A bare KO is the degenerate case: one group, threshold 1.0, so
    modules and KOs run through exactly the same downstream code.
    """
    ids = list(feature_defs)
    kos_used = sorted({k for groups, _ in feature_defs.values() for grp in groups for k in grp})
    kidx = {k: i for i, k in enumerate(kos_used)}
    G = np.zeros((len(genomes), len(kos_used)), dtype=bool)   # genome x KO
    for i, g in enumerate(genomes):
        cols = [kidx[k] for k in ko_sets.get(g, ()) if k in kidx]
        if cols:
            G[i, cols] = True

    M = np.zeros((len(genomes), len(ids)), dtype=np.float32)
    for j, fid in enumerate(ids):
        groups, frac = feature_defs[fid]
        need = max(1, int(np.ceil(frac * len(groups))))
        hits = np.zeros(len(genomes), dtype=np.int16)
        for grp in groups:
            cols = [kidx[k] for k in grp if k in kidx]
            if cols:
                hits += G[:, cols].any(axis=1)
        M[:, j] = (hits >= need).astype(np.float32)
    return pd.DataFrame(M, index=genomes, columns=ids)


def carriage(ab, pres, weights=None, genome_mask=None):
    """features x samples carriage in percent of the (weighted) mapped community.

    ab       : species x sample abundance
    pres     : species x feature 0/1
    weights  : optional per-species multiplier (e.g. 1/size_Mb for the cell-weighted
               variant). Applied to BOTH numerator and denominator so the result stays
               a 0-100 percentage.
    genome_mask : optional boolean over species (e.g. one phylum only).
    """
    A = ab
    Pm = pres
    if genome_mask is not None:
        A = A.loc[genome_mask]
        Pm = Pm.loc[genome_mask]
    if weights is not None:
        w = weights.reindex(A.index).values[:, None]
        A = A * w
    denom = A.sum(axis=0)
    num = Pm.T.values @ A.values          # features x samples
    out = pd.DataFrame(num, index=Pm.columns, columns=A.columns)
    out = 100.0 * out.div(denom.replace(0, np.nan), axis=1)
    return out, denom


# ---------------------------------------------------------------- statistics
def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        prev = min(prev, p[idx] * n / (i + 1))
        q[idx] = prev
    return q


def ranks(a, axis=None):
    return rankdata(a, axis=axis).astype(float)


def partial_spearman(Msub, y, covars):
    """Partial Spearman of each ROW of Msub vs y, controlling covars.

    Verbatim method from the paper's code/phase_e_partial_osmolyte.py (section 2.5):
    rank-transform everything, residualise on [1, ranked covars] with a hat matrix,
    correlate the residuals, t-test with df = n - 2 - len(covars).
    """
    n = len(y)
    X = np.column_stack([np.ones(n)] + [ranks(c) for c in covars])
    H = X @ np.linalg.pinv(X.T @ X) @ X.T
    ry = ranks(y)
    ey = ry - H @ ry
    rM = ranks(Msub, axis=1)
    eM = rM - rM @ H.T
    num = eM @ ey
    den = np.sqrt((eM ** 2).sum(axis=1) * (ey ** 2).sum())
    r = np.divide(num, den, out=np.zeros(Msub.shape[0]), where=den > 0)
    r = np.clip(r, -0.999999, 0.999999)
    df = n - 2 - len(covars)
    tstat = r * np.sqrt(df / (1 - r ** 2))
    p = 2 * tdist.sf(np.abs(tstat), df)
    return r, p


MIN_PAIRS = 6


def wilcoxon_signed_rank_matrix(D, min_n=MIN_PAIRS):
    """Wilcoxon signed-rank for EVERY column of D at once.

    D is pairs x features; NaN marks a pair that is missing for that feature. Rather than
    calling scipy once per feature (which is exact-mode and far too slow at 4,127 features
    x every contrast), this uses the sign-permutation form of the test:

        S   = sum_i sign(d_i) * rank(|d_i|)
        Var = sum_{i: d_i != 0} rank(|d_i|)^2        (conditional on the observed ranks)
        z   = S / sqrt(Var)

    Conditioning on the ranks makes the variance exact under the sign-symmetry null and
    handles ties and zeros without a separate correction, so this is Pratt's treatment of
    zeros (they are ranked but contribute no variance) with an exact tie correction. It
    agrees with scipy's zero_method='pratt', mode='approx' to within the continuity term.

    Returns z, two-sided p, and the number of usable pairs per feature.
    """
    D = np.asarray(D, float)
    if D.ndim == 1:
        D = D[:, None]
    valid = ~np.isnan(D)
    A = np.where(valid, np.abs(D), np.inf)     # NaNs rank last, so valid ranks are right
    Rk = rankdata(A, axis=0)
    sg = np.sign(np.where(valid, D, 0.0))
    S = (sg * Rk).sum(axis=0)
    var = (((sg != 0) * Rk) ** 2).sum(axis=0)
    n = valid.sum(axis=0)
    z = np.divide(S, np.sqrt(var), out=np.zeros_like(S), where=var > 0)
    p = 2.0 * norm.sf(np.abs(z))
    p = np.where((n >= min_n) & (var > 0), p, 1.0)
    return z, p, n


# ---------------------------------------------------------------- pairing
# Moved here from 10_differential.py when the analysis was promoted, so that
# 40_review_diagnostics.py tests the same pairing code the results are built from
# instead of a copy of it. Behaviour is unchanged.
def site_compartment_trip(vals, design):
    """features x sample -> long frame of medians per site x compartment x trip."""
    d = design.reindex(vals.columns)
    keys = pd.MultiIndex.from_arrays([d["site"], d["compartment"], d["trip"]],
                                     names=["site", "compartment", "trip"])
    t = vals.T
    t.index = keys
    return t.groupby(level=["site", "compartment", "trip"]).median()


def compartment_diffs(agg, kind):
    """Paired differences per site, already collapsed over trips.

    Pairing is done WITHIN a visit (site x trip) so a trip1 Surface is never compared
    with a trip3 Deep, then the per-site value is the median over the visits that pair.

    Returns (difference, arm A level, arm B level), all one row per site. A is Surface
    or PlantRoot, B is Deep or bulk (the median of the Surface and Deep medians at the
    same visit); the arms are masked to the visits that actually paired, so the reported
    effect is the difference of the two reported levels.
    """
    wide = agg.unstack("compartment")            # index site x trip, cols feature x comp
    feats = agg.columns
    have = wide.columns.get_level_values(1).unique()
    if kind == "surface_deep":
        if "Surface" not in have or "Deep" not in have:
            return None, None, None
        a = wide.xs("Surface", axis=1, level=1).reindex(columns=feats)
        b = wide.xs("Deep", axis=1, level=1).reindex(columns=feats)
        idx = a.index.union(b.index)
        a, b = a.reindex(idx), b.reindex(idx)
    else:                                        # plantroot vs bulk (median of S and D)
        if "PlantRoot" not in have:
            return None, None, None
        a = wide.xs("PlantRoot", axis=1, level=1).reindex(columns=feats)
        parts = [wide.xs(c, axis=1, level=1).reindex(columns=feats)
                 for c in ("Surface", "Deep") if c in have]
        b = pd.concat(parts).groupby(level=[0, 1]).median() if len(parts) > 1 else parts[0]
        common = a.index.intersection(b.index)
        a, b = a.loc[common], b.loc[common]
    diff = (a - b).dropna(how="all")
    a, b = a.reindex(diff.index), b.reindex(diff.index)
    ok = diff.notna()
    return (diff.groupby(level="site").median(),
            a.where(ok).groupby(level="site").median(),
            b.where(ok).groupby(level="site").median())


def run_compartment(agg, arms=False):
    """Wilcoxon signed-rank over sites, one frame per compartment contrast.

    `arms=True` also reports the two carriage levels being differenced, so an effect in
    percentage points can be read against its baseline.
    """
    rows = {}
    for kind in ("surface_deep", "plantroot_bulk"):
        per_site, arm_a, arm_b = compartment_diffs(agg, kind)
        if per_site is None or per_site.empty:
            continue
        med = per_site.median(axis=0, skipna=True)
        _, p, n = wilcoxon_signed_rank_matrix(per_site.values)
        df = pd.DataFrame(
            {"feature": per_site.columns, "contrast": kind, "n_sites": n,
             "median_diff_pp": med.values, "p": p}, ).set_index("feature")
        if arms:
            df["median_a"] = arm_a.median(axis=0, skipna=True).reindex(df.index).values
            df["median_b"] = arm_b.median(axis=0, skipna=True).reindex(df.index).values
        rows[kind] = df
    return rows
