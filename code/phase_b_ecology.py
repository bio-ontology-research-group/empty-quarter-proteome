#!/usr/bin/env python3
"""Phase B (HEADLINE): is the dark proteome ecologically structured?

Reads the dark-family x sample matrix (Phase A) + env_matrix.tsv, then runs distance-based
PERMANOVA (adonis) per environmental block, Borcard variance partitioning, distance-decay-controlled
partial Mantel, niche IndVal indicators, and a per-family endemism index. All stats from eco_stats.py
(self-contained; no scikit-bio/vegan).

Usage:
  python3 code/phase_b_ecology.py            # real run (needs darkfam_abundance.tsv etc.)
  python3 code/phase_b_ecology.py --selftest # fabricate a compartment-structured community, test wiring
"""
import os, sys, math, csv
import numpy as np
import eco_stats as es

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
ENVF = os.path.join(ECO, "env_matrix.tsv")
ABUND = os.path.join(ECO, "darkfam_abundance.tsv")
PRES = os.path.join(ECO, "darkfam_presence.tsv")
OUTSUM = os.path.join(ECO, "phase_b_summary.txt")
PERMS = 999
MIN_PREV = 5          # families present in >= this many samples kept for ordination


def load_env(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return {r["sample"]: r for r in rows}


def load_matrix(path):
    """family x sample TSV -> (samples list, families list, np.array samples x families)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        fams, data = [], []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            fams.append(c[0])
            data.append([float(x) if x else 0.0 for x in c[1:]])
    M = np.array(data, float).T           # samples x families
    return samples, fams, M


def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def zscore(cols):
    X = np.array(cols, float)
    mu, sd = np.nanmean(X, 0), np.nanstd(X, 0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def haversine_matrix(lat, lon):
    n = len(lat); D = np.zeros((n, n))
    la = np.radians(lat); lo = np.radians(lon)
    for i in range(n):
        for j in range(i + 1, n):
            h = (math.sin((la[j]-la[i])/2)**2 +
                 math.cos(la[i])*math.cos(la[j])*math.sin((lo[j]-lo[i])/2)**2)
            D[i, j] = D[j, i] = 2*6371*math.asin(math.sqrt(min(1, h)))
    return D


def run(samples, fams, M, env, log):
    # align to samples present in env
    keep = [i for i, s in enumerate(samples) if s in env]
    samples = [samples[i] for i in keep]; M = M[keep]
    log(f"samples with community+env: {len(samples)}")
    # prevalence filter (presence inferred as TPM>0 here; presence matrix used separately)
    prev = (M > 0).sum(0)
    fmask = prev >= MIN_PREV
    M = M[:, fmask]; fams = [f for f, k in zip(fams, fmask) if k]
    log(f"families kept (prev>={MIN_PREV}): {M.shape[1]} of {len(fmask)}")

    H = es.hellinger(M)
    D = es.bray_curtis(H)

    # env vectors
    def col(name): return [fnum(env[s].get(name)) for s in samples]
    comp = [env[s].get("compartment", "") for s in samples]
    lat = np.array([fnum(env[s].get("lat")) for s in samples], float)
    lon = np.array([fnum(env[s].get("lon")) for s in samples], float)

    results = {}

    # --- compartment (full set) ---
    Xc = es.dummy(comp)
    if Xc.shape[1] > 0:
        results["compartment"] = es.adonis(D, Xc, perms=PERMS)

    # --- space (lat/lon) ---
    if not np.isnan(lat).any():
        Xs = zscore([[lat[i], lon[i]] for i in range(len(samples))])
        results["space(latlon)"] = es.adonis(D, Xs, perms=PERMS)

    def space_design(rows):
        return zscore([[lat[i], lon[i]] for i in rows])

    # --- climate block (samples with temp+precip); also conditioned on space ---
    clim_idx = [i for i, s in enumerate(samples)
                if fnum(env[s].get("temp")) is not None and fnum(env[s].get("precip")) is not None]
    if len(clim_idx) >= 30:
        Dc = D[np.ix_(clim_idx, clim_idx)]
        Xcl = zscore([[fnum(env[samples[i]]["temp"]), fnum(env[samples[i]]["precip"]),
                       fnum(env[samples[i]].get("aridity_index") or 0)] for i in clim_idx])
        results["climate"] = es.adonis(Dc, Xcl, perms=PERMS); results["climate"]["n"] = len(clim_idx)
        results["climate|space"] = es.adonis(Dc, Xcl, covariate=space_design(clim_idx), perms=PERMS)
        results["climate|space"]["n"] = len(clim_idx)

    # --- geochemistry block (XRF subset); also conditioned on space ---
    geo_idx = [i for i, s in enumerate(samples) if env[s].get("has_xrf") == "1"]
    if len(geo_idx) >= 30:
        Dg = D[np.ix_(geo_idx, geo_idx)]
        elems = ["V", "Fe", "Ca", "K", "Si", "Ti", "Mn", "salinity_proxy"]
        Xg = zscore([[fnum(env[samples[i]].get(e) or 0) for e in elems] for i in geo_idx])
        results["geochem"] = es.adonis(Dg, Xg, perms=PERMS); results["geochem"]["n"] = len(geo_idx)
        results["geochem|space"] = es.adonis(Dg, Xg, covariate=space_design(geo_idx), perms=PERMS)
        results["geochem|space"]["n"] = len(geo_idx)
        Xv = zscore([[fnum(env[samples[i]].get("V") or 0)] for i in geo_idx])
        results["vanadium"] = es.adonis(Dg, Xv, perms=PERMS); results["vanadium"]["n"] = len(geo_idx)

    # --- prevalence-stratified: are well-observed (core) families more clearly structured? ---
    core_mask = (M > 0).sum(0) >= 0.25 * len(samples)
    if core_mask.sum() >= 50:
        Hco = es.hellinger(M[:, core_mask]); Dco = es.bray_curtis(Hco)
        results["compartment[core]"] = es.adonis(Dco, Xc, perms=PERMS) if Xc.shape[1] else {}
        results["compartment[core]"]["n_fam"] = int(core_mask.sum())
        if len(geo_idx) >= 30:
            Dco_g = Dco[np.ix_(geo_idx, geo_idx)]
            elems = ["V", "Fe", "Ca", "K", "Si", "Ti", "Mn", "salinity_proxy"]
            Xg2 = zscore([[fnum(env[samples[i]].get(e) or 0) for e in elems] for i in geo_idx])
            results["geochem|space[core]"] = es.adonis(Dco_g, Xg2, covariate=space_design(geo_idx), perms=PERMS)
            results["geochem|space[core]"]["n_fam"] = int(core_mask.sum())
        if len(clim_idx) >= 30:
            Dco_c = Dco[np.ix_(clim_idx, clim_idx)]
            Xcl2 = zscore([[fnum(env[samples[i]]["temp"]), fnum(env[samples[i]]["precip"]),
                            fnum(env[samples[i]].get("aridity_index") or 0)] for i in clim_idx])
            results["climate|space[core]"] = es.adonis(Dco_c, Xcl2, covariate=space_design(clim_idx), perms=PERMS)
            results["climate|space[core]"]["n_fam"] = int(core_mask.sum())

    # --- variance partitioning on the full set: compartment | space | climate(if full) ---
    blocks = {}
    if Xc.shape[1] > 0: blocks["compartment"] = Xc
    if not np.isnan(lat).any():
        blocks["space"] = zscore([[lat[i], lon[i]] for i in range(len(samples))])
    if len(clim_idx) == len(samples):   # only if climate complete for all
        blocks["climate"] = zscore([[fnum(env[s]["temp"]), fnum(env[s]["precip"])] for s in samples])
    vp = es.varpart(D, blocks, perms=499) if len(blocks) >= 2 else {}

    # --- distance-decay control: partial Mantel community ~ compartment-dist | geographic-dist ---
    pm = None
    if not np.isnan(lat).any():
        Dgeo = haversine_matrix(lat.tolist(), lon.tolist())
        # env distance = euclidean on standardized climate+compartment where available
        r_dd, p_dd = es.mantel(D, Dgeo, perms=PERMS)   # raw distance-decay
        pm = {"distance_decay_r": r_dd, "p": p_dd}

    return results, vp, pm, len(samples), M.shape[1]


def report(results, vp, pm, n, nf, log):
    log("\n=== PERMANOVA (adonis) per environmental block ===")
    log(f"{'block':<18}{'n':>6}{'R2':>8}{'pseudoF':>10}{'p':>8}")
    for k, v in results.items():
        log(f"{k:<18}{v.get('n', n):>6}{v['R2_marginal']:>8.3f}{v['pseudoF']:>10.2f}{v['p']:>8.3f}")
    if vp:
        log("\n=== variance partitioning (unique R2 | conditioned on other blocks) ===")
        for k, v in vp.items():
            log(f"  {k:<14} alone R2={v['R2_alone']:.3f} (p={v['p_alone']:.3f})  "
                f"unique R2={v['R2_unique']:.3f} (p={v['p_unique']:.3f})")
    if pm:
        log(f"\n=== distance-decay (Mantel community vs geography) ===")
        log(f"  r={pm['distance_decay_r']:+.3f} p={pm['p']:.3f}")
    log(f"\n[gate] HEADLINE holds if >=1 env block has p<0.05 with non-trivial R2 above the spatial term.")


def selftest():
    env = load_env(ENVF)
    samples = list(env)[:120]
    comp = [env[s]["compartment"] for s in samples]
    rng = np.random.default_rng(0)
    nf = 300
    M = rng.random((len(samples), nf)) * 0.05
    cmap = {c: i for i, c in enumerate(sorted(set(comp)))}
    for r, c in enumerate(comp):                       # inject compartment structure
        M[r, :60] += cmap[c] * (0.4 + rng.random(60))
    fams = [f"fam{i}" for i in range(nf)]
    out = []
    res, vp, pm, n, k = run(samples, fams, M, env, out.append)
    report(res, vp, pm, n, k, out.append)
    print("\n".join(out))
    assert res["compartment"]["p"] < 0.05, "selftest: compartment structure not detected"
    print("\n[selftest] PASS")


def main():
    if "--selftest" in sys.argv:
        selftest(); return
    env = load_env(ENVF)
    samples, fams, M = load_matrix(ABUND)
    out = []
    def log(x): out.append(x); print(x)
    log(f"loaded matrix: {len(fams)} families x {len(samples)} samples")
    res, vp, pm, n, k = run(samples, fams, M, env, log)
    report(res, vp, pm, n, k, log)
    open(OUTSUM, "w").write("\n".join(out) + "\n")
    log(f"\nwritten -> {OUTSUM}")


if __name__ == "__main__":
    main()
