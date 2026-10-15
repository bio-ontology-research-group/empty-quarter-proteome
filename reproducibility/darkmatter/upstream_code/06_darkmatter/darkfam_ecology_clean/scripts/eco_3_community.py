#!/usr/bin/env python3
"""eco_3 — community-level structure (§2.5 rerun on the 106,042 universe).

Copy of paper-repo code/phase_b_ecology.py (June 2026, 23,382 object) repointed to the
rebuilt basis. Method verbatim: Hellinger + Bray-Curtis on the family x sample matrix,
distance-based PERMANOVA (adonis) per environmental block (compartment / space / climate /
geochemistry, the latter two also conditioned on space), prevalence-stratified core check,
Borcard variance partitioning, distance-decay Mantel. Stats from eco_stats.py (verbatim copy).

Additions vs June (documented in ../README.md): dumps PCoA coordinates + explained variance
(results/pcoa_coords.tsv) and the PERMANOVA block table (results/permanova_blocks.tsv) so the
figure step can be regenerated without rerunning the permutations; the June script only wrote
a text summary.

Hardening 2026-09-20 (two referee objections, neither changes any R2):
  * nesting. The 547 samples sit in 59 sites, up to three positions and three replicates per
    site visit, so free permutation of sample labels treats nested samples as independent
    units. Every block is now tested twice, with the free permutation (column `p`, the
    registered value) and with sample labels permuted only within site (`p_within_site`).
  * sequencing depth. Deeper libraries detect more rare dark families, so a block could be a
    depth gradient. eco_0 writes per-sample depth from the fastp reports; the geochemistry
    beyond-space block and the aridity block are rerun additionally conditioned on log10 depth
    (new rows, the existing rows untouched).

Usage:
  python3 eco_3_community.py            # real run
  python3 eco_3_community.py --selftest # fabricate a compartment-structured community, test wiring
"""
import os, sys, math, csv, shutil, tempfile
import numpy as np
import eco_stats as es

B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter"
REPO = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
ENVF = os.path.join(REPO, "results", "darkmatter_ecology", "env_matrix.tsv")
ABUND = os.path.join(B, "darkfam_basis_clean", "results", "fam_abundance.tsv")
RES = os.path.join(B, "darkfam_ecology_clean", "results")
DEPTHF = os.path.join(RES, "sample_depth.tsv")     # written by eco_0_depth.py
OUTSUM = os.path.join(RES, "community_summary.txt")
PERMS = 999
MIN_PREV = 5          # families present in >= this many samples kept for ordination
# Corrected 2026-09-17: S, Cl, Na added (they carry the family-level signal). V1 = the block as
# first shipped, retained only to reproduce the registered 0.123 / 0.080.
GEO_ELEMS_V1 = ["V", "Fe", "Ca", "K", "Si", "Ti", "Mn", "salinity_proxy"]
GEO_ELEMS = GEO_ELEMS_V1 + ["S", "Cl", "Na"]


def load_env(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return {r["sample"]: r for r in rows}


def load_depth(path):
    """sample -> bases surviving fastp (eco_0_depth.py). Empty dict if the table is absent."""
    if not os.path.exists(path):
        return {}
    rows = csv.DictReader(open(path), delimiter="\t")
    return {r["sample"]: float(r["total_bases_after_filtering"]) for r in rows}


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


def run(samples, fams, M, env, log, outdir=RES):
    # align to samples present in env
    keep = [i for i, s in enumerate(samples) if s in env]
    samples = [samples[i] for i in keep]; M = M[keep]
    log(f"samples with community+env: {len(samples)}")
    # prevalence filter (presence inferred as TPM>0 here)
    prev = (M > 0).sum(0)
    fmask = prev >= MIN_PREV
    M = M[:, fmask]; fams = [f for f, k in zip(fams, fmask) if k]
    log(f"families kept (prev>={MIN_PREV}): {M.shape[1]} of {len(fmask)}")

    H = es.hellinger(M)
    D = es.bray_curtis(H)

    # env vectors
    comp = [env[s].get("compartment", "") for s in samples]
    site = [env[s].get("site", "") for s in samples]       # 59 sites; the nesting unit
    lat = np.array([fnum(env[s].get("lat")) for s in samples], float)
    lon = np.array([fnum(env[s].get("lon")) for s in samples], float)
    if any(x == "" for x in site):
        raise SystemExit("eco_3: env table has samples without a site; strata undefined")
    log(f"sites (permutation strata): {len(set(site))} over {len(samples)} samples")

    def blk(Dm, X, covariate=None, rows=None, perms=PERMS):
        """One block, tested twice: free permutation (the registered p) and restricted to
        within-site exchanges. Same seed and same design both times, so R2 and pseudo-F are the
        free-permutation values unchanged and only p differs. rows = the sample indices Dm was
        built from (None = all)."""
        strata = site if rows is None else [site[i] for i in rows]
        r = es.adonis(Dm, X, covariate=covariate, perms=perms)
        rw = es.adonis(Dm, X, covariate=covariate, perms=perms, strata=strata)
        r["p_within_site"] = rw["p"]
        return r

    results = {}

    # --- compartment (full set) ---
    Xc = es.dummy(comp)
    if Xc.shape[1] > 0:
        results["compartment"] = blk(D, Xc)

    # --- space (lat/lon) ---
    if not np.isnan(lat).any():
        Xs = zscore([[lat[i], lon[i]] for i in range(len(samples))])
        results["space(latlon)"] = blk(D, Xs)

    def space_design(rows):
        return zscore([[lat[i], lon[i]] for i in rows])

    # --- climate block (samples with temp+precip); also conditioned on space ---
    clim_idx = [i for i, s in enumerate(samples)
                if fnum(env[s].get("temp")) is not None and fnum(env[s].get("precip")) is not None]
    if len(clim_idx) >= 30:
        Dc = D[np.ix_(clim_idx, clim_idx)]
        Xcl = zscore([[fnum(env[samples[i]]["temp"]), fnum(env[samples[i]]["precip"]),
                       fnum(env[samples[i]].get("aridity_index") or 0)] for i in clim_idx])
        results["climate"] = blk(Dc, Xcl, rows=clim_idx); results["climate"]["n"] = len(clim_idx)
        results["climate|space"] = blk(Dc, Xcl, covariate=space_design(clim_idx), rows=clim_idx)
        results["climate|space"]["n"] = len(clim_idx)

    # --- geochemistry block (XRF subset); also conditioned on space ---
    geo_idx = [i for i, s in enumerate(samples) if env[s].get("has_xrf") == "1"]
    if len(geo_idx) >= 30:
        Dg = D[np.ix_(geo_idx, geo_idx)]
        # CORRECTION 2026-09-17: S, Cl and Na were missing from this block although they carry
        # the family-level signal; Si, Ti and V are associated with no family and V is the
        # declared null block. GEO_ELEMS is now the corrected 11-axis block. GEO_ELEMS_V1 is
        # kept only so the registered 0.123 / 0.080 stay reproducible.
        Xg = zscore([[fnum(env[samples[i]].get(e) or 0) for e in GEO_ELEMS] for i in geo_idx])
        results["geochem"] = blk(Dg, Xg, rows=geo_idx); results["geochem"]["n"] = len(geo_idx)
        results["geochem|space"] = blk(Dg, Xg, covariate=space_design(geo_idx), rows=geo_idx)
        results["geochem|space"]["n"] = len(geo_idx)
        # compartment control: the XRF cohort is ~51% root-adjacent against 20% surface
        cov_sc = np.hstack([space_design(geo_idx), Xc[geo_idx]]) if Xc.shape[1] else space_design(geo_idx)
        results["geochem|space+compartment"] = blk(Dg, Xg, covariate=cov_sc, rows=geo_idx)
        results["geochem|space+compartment"]["n"] = len(geo_idx)
        Xg1 = zscore([[fnum(env[samples[i]].get(e) or 0) for e in GEO_ELEMS_V1] for i in geo_idx])
        results["geochem[v1_no_salt]"] = blk(Dg, Xg1, rows=geo_idx)
        results["geochem[v1_no_salt]"]["n"] = len(geo_idx)
        results["geochem[v1_no_salt]|space"] = blk(Dg, Xg1, covariate=space_design(geo_idx), rows=geo_idx)
        results["geochem[v1_no_salt]|space"]["n"] = len(geo_idx)
        Xv = zscore([[fnum(env[samples[i]].get("V") or 0)] for i in geo_idx])
        results["vanadium"] = blk(Dg, Xv, rows=geo_idx); results["vanadium"]["n"] = len(geo_idx)

    # --- prevalence-stratified: are well-observed (core) families more clearly structured? ---
    core_mask = (M > 0).sum(0) >= 0.25 * len(samples)
    if core_mask.sum() >= 50:
        Hco = es.hellinger(M[:, core_mask]); Dco = es.bray_curtis(Hco)
        results["compartment[core]"] = blk(Dco, Xc) if Xc.shape[1] else {}
        results["compartment[core]"]["n_fam"] = int(core_mask.sum())
        if len(geo_idx) >= 30:
            Dco_g = Dco[np.ix_(geo_idx, geo_idx)]
            Xg2 = zscore([[fnum(env[samples[i]].get(e) or 0) for e in GEO_ELEMS] for i in geo_idx])
            results["geochem|space[core]"] = blk(Dco_g, Xg2, covariate=space_design(geo_idx), rows=geo_idx)
            results["geochem|space[core]"]["n_fam"] = int(core_mask.sum())
        if len(clim_idx) >= 30:
            Dco_c = Dco[np.ix_(clim_idx, clim_idx)]
            Xcl2 = zscore([[fnum(env[samples[i]]["temp"]), fnum(env[samples[i]]["precip"]),
                            fnum(env[samples[i]].get("aridity_index") or 0)] for i in clim_idx])
            results["climate|space[core]"] = blk(Dco_c, Xcl2, covariate=space_design(clim_idx), rows=clim_idx)
            results["climate|space[core]"]["n_fam"] = int(core_mask.sum())

    # --- sequencing-depth covariate (added 2026-09-20; new rows only, nothing above is touched) ---
    # TPM removes library size from each family's value but not from what the library detected:
    # a deeper library sees more of the rare dark tail, which moves a Bray-Curtis profile. The
    # depth-conditioned rows say whether the geochemistry and aridity signals survive that.
    depth = load_depth(DEPTHF)
    missing = [s for s in samples if s not in depth]
    if missing:
        log(f"[depth] SKIPPED: {len(missing)} samples absent from {DEPTHF} "
            f"(e.g. {missing[:3]}); run eco_0_depth.py first")
    else:
        logd = np.log10(np.array([depth[s] for s in samples], float))
        log(f"[depth] log10 bases after filtering: min {logd.min():.2f} median "
            f"{np.median(logd):.2f} max {logd.max():.2f}")

        def depth_design(rows):
            return zscore([[logd[i]] for i in rows])

        allrows = list(range(len(samples)))
        results["depth(log10)"] = blk(D, depth_design(allrows))
        if len(geo_idx) >= 30:
            cov_sd = np.hstack([space_design(geo_idx), depth_design(geo_idx)])
            results["geochem|space+depth"] = blk(Dg, Xg, covariate=cov_sd, rows=geo_idx)
            results["geochem|space+depth"]["n"] = len(geo_idx)
        if len(clim_idx) >= 30:
            # aridity on its own, and the whole climate block, each beyond space and beyond
            # space+depth. The aridity cohort is the climate cohort (456; coverage is identical).
            Xar = zscore([[fnum(env[samples[i]]["aridity_index"])] for i in clim_idx])
            cov_s = space_design(clim_idx)
            cov_sd = np.hstack([cov_s, depth_design(clim_idx)])
            results["aridity|space"] = blk(Dc, Xar, covariate=cov_s, rows=clim_idx)
            results["aridity|space"]["n"] = len(clim_idx)
            results["aridity|space+depth"] = blk(Dc, Xar, covariate=cov_sd, rows=clim_idx)
            results["aridity|space+depth"]["n"] = len(clim_idx)
            results["climate|space+depth"] = blk(Dc, Xcl, covariate=cov_sd, rows=clim_idx)
            results["climate|space+depth"]["n"] = len(clim_idx)

    # --- variance partitioning ---
    # CORRECTION 2026-09-17: the previous version admitted climate only when it covered every
    # sample (456 of 547, so it was always dropped) and never offered geochemistry, leaving a
    # two-block partition of compartment and space that just reproduced the PERMANOVA, although
    # Methods described a partition of the environmental blocks. It now runs on the XRF cohort,
    # where compartment, space and geochemistry are all complete.
    blocks = {}
    if Xc.shape[1] > 0: blocks["compartment"] = Xc[geo_idx]
    if not np.isnan(lat).any(): blocks["space"] = space_design(geo_idx)
    blocks["geochem"] = zscore([[fnum(env[samples[i]].get(e) or 0) for e in GEO_ELEMS]
                                for i in geo_idx])
    # Same two permutation schemes as the PERMANOVA blocks: R2_alone / R2_unique are the free
    # values unchanged, the within-site p is appended beside each free p.
    vp = {}
    if len(blocks) >= 2:
        Dvp = D[np.ix_(geo_idx, geo_idx)]
        vp = es.varpart(Dvp, blocks, perms=499)
        vp_ws = es.varpart(Dvp, blocks, perms=499, strata=[site[i] for i in geo_idx])
        for nm in vp:
            vp[nm]["p_alone_within_site"] = vp_ws[nm]["p_alone"]
            vp[nm]["p_unique_within_site"] = vp_ws[nm]["p_unique"]

    # --- distance-decay control: Mantel community vs geography ---
    pm = None
    if not np.isnan(lat).any():
        Dgeo = haversine_matrix(lat.tolist(), lon.tolist())
        r_dd, p_dd = es.mantel(D, Dgeo, perms=PERMS)   # raw distance-decay
        pm = {"distance_decay_r": r_dd, "p": p_dd}

    # --- dumps for the figure step (addition vs June) ---
    coords, explained = es.pcoa(D, k=3)
    with open(os.path.join(outdir, "pcoa_coords.tsv"), "w") as o:
        o.write("sample\tPC1\tPC2\tPC3\tcompartment\tsite\ttrip\n")
        for i, s in enumerate(samples):
            r = env[s]
            o.write(f"{s}\t{coords[i,0]:.6f}\t{coords[i,1]:.6f}\t{coords[i,2]:.6f}\t"
                    f"{r.get('compartment','')}\t{r.get('site','')}\t{r.get('trip_num','')}\n")
    with open(os.path.join(outdir, "pcoa_explained.tsv"), "w") as o:
        o.write("axis\texplained\n")
        for i, e in enumerate(explained, 1):
            o.write(f"PC{i}\t{e:.6f}\n")

    return results, vp, pm, len(samples), M.shape[1]


def report(results, vp, pm, n, nf, log, outdir=RES):
    log("\n=== PERMANOVA (adonis) per environmental block ===")
    log(f"{'block':<28}{'n':>6}{'R2':>8}{'pseudoF':>10}{'p':>8}{'p_within_site':>15}")
    with open(os.path.join(outdir, "permanova_blocks.tsv"), "w") as o:
        # columns through n_fam are the registered layout and are not reordered; p_within_site is
        # the same test with sample labels permuted only within site (59 sites).
        o.write("block\tn\tR2_marginal\tR2_partial\tpseudoF\tp\tn_fam\tp_within_site\n")
        for k, v in results.items():
            pw = v.get("p_within_site", float("nan"))
            log(f"{k:<28}{v.get('n', n):>6}{v['R2_marginal']:>8.3f}{v['pseudoF']:>10.2f}"
                f"{v['p']:>8.3f}{pw:>15.3f}")
            o.write(f"{k}\t{v.get('n', n)}\t{v['R2_marginal']:.4f}\t{v.get('R2_partial', float('nan')):.4f}\t"
                    f"{v['pseudoF']:.3f}\t{v['p']:.4f}\t{v.get('n_fam', '')}\t{pw:.4f}\n")
    if vp:
        log("\n=== variance partitioning (unique R2 | conditioned on other blocks) ===")
        with open(os.path.join(outdir, "varpart.tsv"), "w") as o:
            o.write("block\tR2_alone\tp_alone\tR2_unique\tp_unique\t"
                    "p_alone_within_site\tp_unique_within_site\n")
            for k, v in vp.items():
                log(f"  {k:<14} alone R2={v['R2_alone']:.3f} (p={v['p_alone']:.3f}, "
                    f"within-site p={v['p_alone_within_site']:.3f})  "
                    f"unique R2={v['R2_unique']:.3f} (p={v['p_unique']:.3f}, "
                    f"within-site p={v['p_unique_within_site']:.3f})")
                o.write(f"{k}\t{v['R2_alone']:.4f}\t{v['p_alone']:.4f}\t"
                        f"{v['R2_unique']:.4f}\t{v['p_unique']:.4f}\t"
                        f"{v['p_alone_within_site']:.4f}\t{v['p_unique_within_site']:.4f}\n")
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
    # fabricated tables go to a scratch directory: the selftest must never touch results/.
    tmp = tempfile.mkdtemp(prefix="eco3_selftest_")
    res, vp, pm, n, k = run(samples, fams, M, env, out.append, outdir=tmp)
    report(res, vp, pm, n, k, out.append, outdir=tmp)
    print("\n".join(out))
    print(f"\n[selftest] fabricated tables written to {tmp} (results/ untouched)")
    assert res["compartment"]["p"] < 0.05, "selftest: compartment structure not detected"
    # compartment varies inside a site, so the restricted test must still see it; lat/lon are
    # constant inside a site, so the restricted test of the spatial term must return p=1.
    assert res["compartment"]["p_within_site"] < 0.05, "selftest: compartment lost under within-site perms"
    assert res["space(latlon)"]["p_within_site"] == 1.0, \
        "selftest: a site-constant term must be untestable within site (p=1)"
    shutil.rmtree(tmp, ignore_errors=True)
    print("[selftest] PASS")


def main():
    if "--selftest" in sys.argv:
        selftest(); return
    env = load_env(ENVF)
    samples, fams, M = load_matrix(ABUND)
    out = []
    def log(x): out.append(x); print(x, flush=True)
    log(f"loaded matrix: {len(fams)} families x {len(samples)} samples")
    res, vp, pm, n, k = run(samples, fams, M, env, log)
    report(res, vp, pm, n, k, log)
    open(OUTSUM, "w").write("\n".join(out) + "\n")
    log(f"\nwritten -> {OUTSUM}")


if __name__ == "__main__":
    main()
