#!/usr/bin/env python3
"""Self-contained community-ecology statistics in numpy/scipy (no scikit-bio / vegan available).

Implements the McArdle & Anderson (2001) distance-based PERMANOVA / adonis on a Gower-centred
distance matrix, with covariate conditioning (partial adonis), Borcard variance partitioning,
classical PCoA, Mantel + partial Mantel, and IndVal indicator analysis. Used by phase_b_ecology.py.

Run `python3 code/eco_stats.py` to execute the built-in self-tests (synthetic data with known
structure: a real env effect must be detected; a null must not).
"""
import numpy as np
from scipy.spatial.distance import pdist, squareform


# ---------- distances ----------
def bray_curtis(M):
    """M: samples x features (abundances >=0). Returns square BC distance matrix."""
    return squareform(pdist(M, metric="braycurtis"))


def hellinger(M):
    """Hellinger transform (sqrt of relative abundance); good for compositional community data."""
    rs = M.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    return np.sqrt(M / rs)


# ---------- Gower centring ----------
def gower(D):
    """G = -1/2 J D^2 J  (double-centred). trace(G) = total sum of squares."""
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J


def _hat(X):
    """Hat (projection) matrix H = X (X'X)^-1 X' for design X (n x p, include intercept)."""
    XtX = X.T @ X
    return X @ np.linalg.pinv(XtX) @ X.T


# ---------- restricted permutation ----------
def strata_groups(strata):
    """Positions of each stratum level, in first-appearance order. strata: length-n labels."""
    groups, order = {}, []
    for i, s in enumerate(strata):
        if s not in groups:
            groups[s] = []
            order.append(s)
        groups[s].append(i)
    return [np.asarray(groups[s], dtype=int) for s in order]


def permute_within(rng, groups, n):
    """Permutation of 0..n-1 that only exchanges positions inside a stratum.

    Used for the restricted PERMANOVA test: samples are nested in sites, so free permutation of
    sample labels treats replicates and positions within one site visit as independent units
    (pseudoreplication). Permuting within site keeps every sample in its own site, so the null
    is exchangeability within site rather than across the whole design.
    """
    pr = np.empty(n, dtype=int)
    for g in groups:
        pr[g] = rng.permutation(g)
    return pr


# ---------- PERMANOVA / adonis (distance-based) ----------
def adonis(D, X, covariate=None, perms=999, seed=1, strata=None):
    """Distance-based PERMANOVA (McArdle-Anderson). D: n x n distance matrix.
    X: n x k design for the term(s) of interest (NO intercept column; added here).
    covariate: optional n x c design conditioned out (partial adonis). Returns dict with R2,
    pseudo-F, p (by permutation of the Gower matrix). R2 is partial when covariate given.

    strata: optional length-n labels (e.g. site). When given, sample labels are permuted only
    WITHIN a stratum; when None the permutation is free over all n samples, bit-identical to
    the unstratified version (same seed, same rng draws). Only p changes: the observed
    statistic, R2 and pseudo-F do not depend on the permutation scheme.
    """
    n = D.shape[0]
    G = gower(D)
    ones = np.ones((n, 1))
    if covariate is not None:
        Xc = np.hstack([ones, covariate])
    else:
        Xc = ones
    Xf = np.hstack([Xc, X])               # full model
    Hc, Hf = _hat(Xc), _hat(Xf)
    # SS explained by X beyond covariate; residual = beyond full model.
    # trace(H G H) = trace(G H) = sum(G*H) since H is symmetric & idempotent -> O(n^2) per perm.
    def ss(H):
        return np.sum(G * H)
    trG = np.trace(G)
    ss_total = trG
    ss_cov = ss(Hc)
    ss_full = ss(Hf)
    ss_term = ss_full - ss_cov            # added by X
    ss_res = ss_total - ss_full
    p_term = X.shape[1]
    df_res = n - Xf.shape[1]
    F = (ss_term / p_term) / (ss_res / df_res) if ss_res > 0 and df_res > 0 else np.nan
    r2 = ss_term / (ss_total - ss_cov) if (ss_total - ss_cov) > 0 else np.nan  # partial R2
    r2_marg = ss_term / ss_total if ss_total > 0 else np.nan
    # permutation: row/col-permute the Gower matrix (trace(Gp)=trace(G); use sum(Gp*H) for trace(Gp H))
    rng = np.random.default_rng(seed)
    groups = strata_groups(strata) if strata is not None else None
    # Count a permutation as >= the observed F within rounding error, as vegan's permutest does.
    # This matters for restricted permutation: a term that is constant inside every stratum is
    # exactly invariant under within-stratum exchange, so every Fp equals F and the honest answer
    # is p=1 (the term is not testable at this level). Without the tolerance, summation noise of
    # order 1e-10 would split those ties at random and manufacture a p-value.
    tol = 1e-9 * abs(F) if np.isfinite(F) else 0.0
    ge = 1
    for _ in range(perms):
        pr = rng.permutation(n) if groups is None else permute_within(rng, groups, n)
        Gp = G[np.ix_(pr, pr)]
        ss_full_p = np.sum(Gp * Hf)
        ss_cov_p = np.sum(Gp * Hc)
        ss_res_p = trG - ss_full_p
        Fp = ((ss_full_p - ss_cov_p) / p_term) / (ss_res_p / df_res) if ss_res_p > 0 else -np.inf
        if Fp >= F - tol:
            ge += 1
    return {"R2_partial": r2, "R2_marginal": r2_marg, "pseudoF": F,
            "p": ge / (perms + 1), "ss_term": ss_term, "ss_total": ss_total,
            "perm_scheme": "free" if groups is None else "within_strata"}


def dummy(labels):
    """One-hot design (drop first level) for a categorical factor; returns n x (k-1)."""
    levs = sorted(set(labels))
    if len(levs) < 2:
        return np.zeros((len(labels), 0))
    idx = {l: i for i, l in enumerate(levs)}
    X = np.zeros((len(labels), len(levs)))
    for r, l in enumerate(labels):
        X[r, idx[l]] = 1.0
    return X[:, 1:]


def varpart(D, blocks, perms=499, seed=1, strata=None):
    """Borcard variance partitioning for up to 3 named env blocks (dict name->design matrix).
    Reports each block's marginal R2 (alone) and unique/partial R2 (conditioned on the others).
    strata: passed straight to adonis; None keeps the free permutation (identical results)."""
    names = list(blocks)
    out = {}
    for nm in names:
        others = [blocks[o] for o in names if o != nm]
        cov = np.hstack(others) if others else None
        a_alone = adonis(D, blocks[nm], covariate=None, perms=perms, seed=seed, strata=strata)
        a_uniq = adonis(D, blocks[nm], covariate=cov, perms=perms, seed=seed, strata=strata)
        out[nm] = {"R2_alone": a_alone["R2_marginal"], "p_alone": a_alone["p"],
                   "R2_unique": a_uniq["R2_partial"], "p_unique": a_uniq["p"]}
    return out


# ---------- PCoA (classical MDS) ----------
def pcoa(D, k=3):
    G = gower(D)
    w, V = np.linalg.eigh(G)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    pos = w > 0
    coords = V[:, :k] * np.sqrt(np.abs(w[:k]))
    explained = w[pos] / w[pos].sum()
    return coords, explained[:k]


# ---------- Mantel + partial Mantel ----------
def mantel(A, B, perms=999, seed=1):
    iu = np.triu_indices_from(A, 1)
    a, b = A[iu], B[iu]
    if a.std() == 0 or b.std() == 0:
        return np.nan, np.nan
    r = np.corrcoef(a, b)[0, 1]
    rng = np.random.default_rng(seed); n = A.shape[0]; ge = 1
    for _ in range(perms):
        p = rng.permutation(n)
        if abs(np.corrcoef(A[np.ix_(p, p)][iu], b)[0, 1]) >= abs(r):
            ge += 1
    return r, ge / (perms + 1)


def partial_mantel(A, B, C, perms=999, seed=1):
    """Correlation of A,B controlling for C (all square distance matrices)."""
    iu = np.triu_indices_from(A, 1)
    a, b, c = A[iu], B[iu], C[iu]
    def resid(x, z):
        z1 = np.vstack([np.ones_like(z), z]).T
        beta = np.linalg.lstsq(z1, x, rcond=None)[0]
        return x - z1 @ beta
    ra, rb = resid(a, c), resid(b, c)
    r = np.corrcoef(ra, rb)[0, 1]
    rng = np.random.default_rng(seed); n = A.shape[0]; ge = 1
    for _ in range(perms):
        p = rng.permutation(n)
        Ap = A[np.ix_(p, p)][iu]
        rap = resid(Ap, c)
        if abs(np.corrcoef(rap, rb)[0, 1]) >= abs(r):
            ge += 1
    return r, ge / (perms + 1)


# ---------- IndVal (Dufrene-Legendre) ----------
def indval(M, labels, fams):
    """Indicator value per feature per group. M: samples x features (abundance). Returns list of
    (family, best_group, indval, A_specificity, B_fidelity)."""
    labels = np.asarray(labels)
    groups = sorted(set(labels))
    res = []
    colsum = M.sum(axis=0) + 1e-12
    for j in range(M.shape[1]):
        best = (None, 0.0, 0.0, 0.0)
        for g in groups:
            mask = labels == g
            A = M[mask, j].sum() / colsum[j]                  # specificity
            B = (M[mask, j] > 0).sum() / max(mask.sum(), 1)   # fidelity
            iv = A * B
            if iv > best[1]:
                best = (g, iv, A, B)
        res.append((fams[j], best[0], best[1], best[2], best[3]))
    return res


# ---------- self-tests ----------
def _selftest():
    rng = np.random.default_rng(0)
    n, p = 60, 200
    env = np.concatenate([np.zeros(30), np.ones(30)])         # 2 groups
    # structured community: 40 families respond to env, 160 noise
    M = rng.random((n, p)) * 0.1
    for j in range(40):
        M[:, j] += env * (0.5 + rng.random())
    M = hellinger(M)
    D = bray_curtis(M)
    X = dummy(env.tolist())
    a = adonis(D, X, perms=499)
    print(f"[test] structured env: R2={a['R2_marginal']:.3f} F={a['pseudoF']:.2f} p={a['p']:.3f}  (expect p<0.05, R2 high)")
    # null env
    enull = rng.permutation(env)
    Mn = hellinger(rng.random((n, p)))
    an = adonis(bray_curtis(Mn), dummy(enull.tolist()), perms=499)
    print(f"[test] null env:       R2={an['R2_marginal']:.3f} F={an['pseudoF']:.2f} p={an['p']:.3f}  (expect p NOT <<0.05)")
    # varpart sanity (env block real, random block null)
    vp = varpart(D, {"env": X, "rand": rng.random((n, 2))}, perms=299)
    print(f"[test] varpart env unique R2={vp['env']['R2_unique']:.3f} p={vp['env']['p_unique']:.3f} ;"
          f" rand unique R2={vp['rand']['R2_unique']:.3f} p={vp['rand']['p_unique']:.3f}")
    assert a["p"] < 0.05 and a["R2_marginal"] > 0.1, "structured env not detected"

    # --- restricted permutation (strata) ---
    # 1. the permutation itself: every sample stays inside its own stratum, so each stratum keeps
    #    exactly its own members and therefore its size, and the whole thing is still a bijection.
    strata = [f"site{i % 12}" for i in range(n)]           # 12 unequal-use strata over 60 samples
    strata[0] = strata[1] = "singleton_pair"               # sizes are not all equal
    groups = strata_groups(strata)
    sizes = {s: strata.count(s) for s in set(strata)}
    rngp = np.random.default_rng(7)
    moved = 0
    for _ in range(200):
        pr = permute_within(rngp, groups, n)
        assert sorted(pr.tolist()) == list(range(n)), "strata permutation is not a bijection"
        for g in groups:
            assert set(pr[g].tolist()) == set(g.tolist()), "sample left its stratum"
        permuted_sizes = {}
        for i in range(n):
            permuted_sizes[strata[pr[i]]] = permuted_sizes.get(strata[pr[i]], 0) + 1
        assert permuted_sizes == sizes, "strata permutation changed group sizes"
        moved += int((pr != np.arange(n)).sum())
    print(f"[test] strata permutation: {len(groups)} strata, sizes {sorted(sizes.values())} preserved "
          f"over 200 draws, mean {moved/200:.1f} of {n} positions moved  (expect sizes fixed, movement > 0)")
    assert moved > 0, "strata permutation never moved anything"

    # 2. adonis wiring: strata=None must stay bit-identical to the unstratified call, and a term
    #    that is constant within stratum must become untestable (p -> 1) under within-strata perms.
    a_free = adonis(D, X, perms=499, strata=None)
    assert (a_free["p"], a_free["pseudoF"]) == (a["p"], a["pseudoF"]), "strata=None changed the free test"
    blk = [f"s{i // 5}" for i in range(n)]                  # 12 blocks of 5, env constant in block
    a_str = adonis(D, X, perms=499, strata=blk)
    print(f"[test] strata adonis: R2={a_str['R2_marginal']:.3f} (free {a_free['R2_marginal']:.3f}, "
          f"identical) F={a_str['pseudoF']:.2f} p={a_str['p']:.3f}  (expect R2/F unchanged, p high: "
          f"env is a block-level term)")
    assert abs(a_str["R2_marginal"] - a_free["R2_marginal"]) < 1e-12, "strata changed R2"
    assert abs(a_str["pseudoF"] - a_free["pseudoF"]) < 1e-12, "strata changed pseudo-F"
    assert a_str["perm_scheme"] == "within_strata" and a_free["perm_scheme"] == "free"
    assert a_str["p"] > a_free["p"], "restricted test not more conservative for a block-level term"

    # 3. a term that is CONSTANT inside every stratum is exactly invariant under within-stratum
    #    exchange, so the restricted test must return p=1, not a tie broken by rounding noise.
    Xconst = np.array([[float(i // 5)] for i in range(n)])   # one value per block of 5
    a_const = adonis(D, Xconst, perms=199, strata=blk)
    print(f"[test] strata, term constant within stratum: p={a_const['p']:.3f}  (expect exactly 1.000)")
    assert a_const["p"] == 1.0, "constant-within-stratum term did not give p=1"
    print("[test] PASS")


if __name__ == "__main__":
    _selftest()
