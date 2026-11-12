#!/usr/bin/env python3
"""Are exported sequence-dark protein families shared between desert (EQ) Thermoplasmatota
genomes and external Thermoplasmatota genomes?

Inputs: MMseqs2 families over all proteins of the 1,957 genomes (stage 1), the sensitive search of
every EQ protein against every external protein (stage 2), and the harmonised producer tables
(loaded with the loader of analysis1_export/export_comparison.py).

Conventions
- Family class (primary): strict majority of the focal (EQ) members dark and strict majority with a
  predicted signal peptide. Protein class: the protein's own dark and signal-peptide status.
- "Shared" with a target genome set: (cluster) the family has at least one member in a target genome;
  (search) at least one focal member has a search hit in a target genome (e-value <= 1e-5, query
  coverage >= 0.5); (either); (search_gain) search hit among families without a clustered member.
- Uncertainty: block bootstrap over focal blocks (named genus, else ANI95 group). Target genomes are
  held fixed. Class and prevalence are recomputed in each draw; the family length bin is fixed.
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CLASSES = ["dark_exported", "dark_not_exported", "nondark_exported", "nondark_not_exported"]
METHODS = ["cluster", "search", "either", "search_gain_if_cluster_unshared"]
GB_EDGES = [2, 5, 10]          # focal genomes: 1, 2-4, 5-9, >=10
BB_EDGES = [2, 3, 5, 10]       # focal blocks: 1, 2, 3-4, 5-9, >=10
NG, NB = len(GB_EDGES) + 1, len(BB_EDGES) + 1


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(4 << 20), b""):
            h.update(b)
    return h.hexdigest()


def q95(x, axis=0):
    with np.errstate(all="ignore"):
        return np.nanpercentile(x, [2.5, 97.5], axis=axis)


def ratio(a, b):
    with np.errstate(all="ignore"):
        return np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)


class Focal:
    """Sharing statistics of the families of a focal genome set with fixed target genome sets."""

    def __init__(self, fam, genome, dark, sp, aa_len, genome_block, shared, hits, lenbins):
        # fam: global family code per focal gene; genome: focal genome index per gene
        self.nl = len(lenbins) - 1
        uniq, self.f = np.unique(fam, return_inverse=True)
        self.fam_global = uniq
        self.F = len(uniq)
        self.g = genome
        self.G = len(genome_block)
        self.gb = genome_block
        self.dark = dark.astype(float)
        self.sp = sp.astype(float)
        self.Sh = shared[uniq].astype(np.float64)          # F x T
        self.T = self.Sh.shape[1]
        self.hits = hits                                     # n_genes x T bool, or None
        if hits is not None:
            self.hg, self.ht = np.nonzero(hits)
        pg = np.unique(self.f.astype(np.int64) * self.G + genome)
        self.pg_f, self.pg_g = pg // self.G, pg % self.G
        nbk = genome_block.max() + 1
        pb = np.unique(self.f.astype(np.int64) * nbk + genome_block[genome])
        self.pb_f, self.pb_b = pb // nbk, pb % nbk
        self.nblocks = nbk
        med = pd.Series(aa_len).groupby(self.f).median().to_numpy()
        self.len_f = np.digitize(med, lenbins[1:-1])
        self.med_len = med
        self.len_gene = np.digitize(aa_len, lenbins[1:-1])
        self.K = self.nl * NG * NB * 4
        self.own_cls = np.where(dark, 0, 2) + np.where(sp, 0, 1)

    def fam_state(self, wg):
        wgene = wg[self.g]
        n = np.bincount(self.f, wgene, self.F)
        d = np.bincount(self.f, wgene * self.dark, self.F)
        p = np.bincount(self.f, wgene * self.sp, self.F)
        cls = np.where(2 * d > n, 0, 2) + np.where(2 * p > n, 0, 1)
        ng = np.bincount(self.pg_f, (wg > 0)[self.pg_g].astype(float), self.F)
        wbpos = (np.bincount(self.gb, wg, self.nblocks) > 0).astype(float)
        nb = np.bincount(self.pb_f, wbpos[self.pb_b], self.F)
        return n, d, p, cls, ng, nb, wgene

    def code(self, lenbin, ng, nb, cls):
        return ((lenbin * NG + np.digitize(ng, GB_EDGES)) * NB + np.digitize(nb, BB_EDGES)) * 4 + cls

    def search_shared(self, wgene):
        if self.hits is None:
            return np.zeros((self.F, self.T))
        h = np.bincount(self.f[self.hg].astype(np.int64) * self.T + self.ht, wgene[self.hg],
                        self.F * self.T)
        return (h.reshape(self.F, self.T) > 0).astype(np.float64)

    def _tab(self, code, mat, nrow):
        """Sum the columns of mat (n x T) within code -> (K x T)."""
        T = mat.shape[1]
        idx = (code[:, None].astype(np.int64) * T + np.arange(T)[None, :]).ravel()
        return np.bincount(idx, mat.ravel(), self.K * T).reshape(self.K, T)

    def family_counts(self, wg):
        """(K, M, T, 2[num, den]) family counts for genome weights wg."""
        n, d, p, cls, ng, nb, wgene = self.fam_state(wg)
        pres = n > 0
        code = self.code(self.len_f, ng, nb, cls)[pres]
        Sh = self.Sh[pres]
        Hs = self.search_shared(wgene)[pres]
        one = np.ones_like(Sh)
        X = np.empty((self.K, len(METHODS), self.T, 2))
        den = self._tab(code, one, None)
        X[:, 0, :, 0] = self._tab(code, Sh, None)
        X[:, 1, :, 0] = self._tab(code, Hs, None)
        X[:, 2, :, 0] = self._tab(code, np.maximum(Sh, Hs), None)
        X[:, 3, :, 0] = self._tab(code, Hs * (1 - Sh), None)
        X[:, 0, :, 1] = X[:, 1, :, 1] = X[:, 2, :, 1] = den
        X[:, 3, :, 1] = self._tab(code, 1 - Sh, None)
        return X

    def protein_tensor(self):
        """(G, K, M, T, 2) protein counts per focal genome; strata use full-sample family prevalence."""
        n, d, p, cls, ng, nb, _ = self.fam_state(np.ones(self.G))
        self.full = dict(n=n, d=d, p=p, cls=cls, ng=ng, nb=nb)
        code = self.code(self.len_gene, ng[self.f], nb[self.f], self.own_cls)
        Sh = self.Sh[self.f]
        Hs = self.hits.astype(np.float64) if self.hits is not None else np.zeros_like(Sh)
        gcode = self.g.astype(np.int64) * self.K + code
        K0, self.K = self.K, self.G * self.K
        A = np.empty((self.G * K0, len(METHODS), self.T, 2))
        den = self._tab(gcode, np.ones_like(Sh), None)
        A[:, 0, :, 0] = self._tab(gcode, Sh, None)
        A[:, 1, :, 0] = self._tab(gcode, Hs, None)
        A[:, 2, :, 0] = self._tab(gcode, np.maximum(Sh, Hs), None)
        A[:, 3, :, 0] = self._tab(gcode, Hs * (1 - Sh), None)
        A[:, 0, :, 1] = A[:, 1, :, 1] = A[:, 2, :, 1] = den
        A[:, 3, :, 1] = self._tab(gcode, 1 - Sh, None)
        self.K = K0
        return A.reshape(self.G, K0, len(METHODS), self.T, 2)

    def run(self, draws, seed):
        A = self.protein_tensor()
        ones = np.ones(self.G)
        point = {"family": self.family_counts(ones), "protein": np.tensordot(ones, A, 1)}
        rng = np.random.default_rng(seed)
        k = self.nblocks
        boot = {lv: np.empty((draws,) + point[lv].shape, dtype=np.float32) for lv in point}
        for i in range(draws):
            wb = np.bincount(rng.integers(0, k, k), minlength=k)
            wg = wb[self.gb].astype(float)
            boot["family"][i] = self.family_counts(wg)
            boot["protein"][i] = np.tensordot(wg, A, 1)
        return point, boot


def filters(nl, lenbins):
    """Named boolean masks over the (lenbin, genome bin, block bin) strata."""
    L, Gb, Bb = np.meshgrid(np.arange(nl), np.arange(NG), np.arange(NB), indexing="ij")
    f = {"all": L >= 0,
         "focal_genomes_eq1": Gb == 0, "focal_genomes_ge2": Gb >= 1, "focal_genomes_ge5": Gb >= 2,
         "focal_genomes_ge10": Gb >= 3,
         "focal_blocks_eq1": Bb == 0, "focal_blocks_ge2": Bb >= 1, "focal_blocks_ge3": Bb >= 2,
         "focal_blocks_ge5": Bb >= 3, "focal_blocks_ge10": Bb >= 4,
         "widespread_ge10genomes_ge3blocks": (Gb >= 3) & (Bb >= 2)}
    for i in range(nl):
        hi = lenbins[i + 1]
        lab = f"len_{lenbins[i]}_{hi if hi < 10**8 else 'inf'}"
        f[lab] = L == i
        f[lab + "_focal_genomes_ge2"] = (L == i) & (Gb >= 1)
    return {k: v.ravel() for k, v in f.items()}


def summarise(point, boot, flt, targets, key, methods, nl):
    """Tidy fraction and contrast rows from stratum counts."""
    rows, crows, srows = [], [], []
    for lv in point:
        P = point[lv].reshape(nl * NG * NB, 4, len(METHODS), len(targets), 2)
        B = boot[lv].reshape((boot[lv].shape[0],) + P.shape)
        for fname, mask in flt.items():
            p = P[mask].sum(0)                       # class, M, T, 2
            b = B[:, mask].sum(1).astype(np.float64)  # draws, class, M, T, 2
            fr, bfr = ratio(p[..., 0], p[..., 1]), ratio(b[..., 0], b[..., 1])
            ci = q95(bfr)
            for mi, m in enumerate(METHODS):
                if m not in methods:
                    continue
                for ti, t in enumerate(targets):
                    for ci_, c in enumerate(CLASSES):
                        rows.append({**key, "level": lv, "filter": fname, "method": m, "target": t,
                                     "class": c, "n": p[ci_, mi, ti, 1], "n_shared": p[ci_, mi, ti, 0],
                                     "fraction_shared": fr[ci_, mi, ti],
                                     "ci_low": ci[0, ci_, mi, ti], "ci_high": ci[1, ci_, mi, ti]})
                    for cj in (1, 2, 3):
                        d = fr[0, mi, ti] - fr[cj, mi, ti]
                        bd = bfr[:, 0, mi, ti] - bfr[:, cj, mi, ti]
                        r = ratio(fr[0, mi, ti], fr[cj, mi, ti])
                        br = ratio(bfr[:, 0, mi, ti], bfr[:, cj, mi, ti])
                        ru = ratio(1 - fr[0, mi, ti], 1 - fr[cj, mi, ti])
                        bru = ratio(1 - bfr[:, 0, mi, ti], 1 - bfr[:, cj, mi, ti])
                        qd, qr, qu = q95(bd), q95(br), q95(bru)
                        crows.append({**key, "level": lv, "filter": fname, "method": m, "target": t,
                                      "comparison_class": CLASSES[cj],
                                      "fraction_dark_exported": fr[0, mi, ti],
                                      "fraction_comparison": fr[cj, mi, ti],
                                      "difference": d, "difference_ci_low": qd[0], "difference_ci_high": qd[1],
                                      "ratio_shared": float(r), "ratio_shared_ci_low": qr[0],
                                      "ratio_shared_ci_high": qr[1],
                                      "ratio_unshared": float(ru), "ratio_unshared_ci_low": qu[0],
                                      "ratio_unshared_ci_high": qu[1],
                                      "bootstrap_p_two_sided": float(min(1.0, 2 * min((bd <= 0).mean(), (bd >= 0).mean())))})
        # Direct standardisation of each comparison class to the dark-exported distribution over
        # length bin x focal-genome bin (block bins pooled).
        S = P.reshape(nl, NG, NB, 4, len(METHODS), len(targets), 2).sum(2).reshape(nl * NG, 4, len(METHODS), len(targets), 2)
        SB = B.reshape((B.shape[0], nl, NG, NB, 4, len(METHODS), len(targets), 2)).sum(3).reshape(
            (B.shape[0], nl * NG, 4, len(METHODS), len(targets), 2)).astype(np.float64)
        gsel = {"all": np.ones(nl * NG, bool),
                "focal_genomes_ge2": (np.arange(nl * NG) % NG) >= 1}
        for fname, gm_ in gsel.items():
            for mi, m in enumerate(METHODS[:3]):
                if m not in methods:
                    continue
                for ti, t in enumerate(targets):
                    for cj in (1, 2, 3):
                        def std(X):
                            de_n, de_d = X[..., 0, mi, ti, 0], X[..., 0, mi, ti, 1]
                            c_n, c_d = X[..., cj, mi, ti, 0], X[..., cj, mi, ti, 1]
                            use = (de_d > 0) & (c_d > 0) & gm_
                            w = np.where(use, de_d, 0.0)
                            adj_c = np.nansum(w * np.nan_to_num(ratio(c_n, c_d)), -1) / w.sum(-1)
                            de = np.where(use, de_n, 0.0).sum(-1) / w.sum(-1)
                            cov = w.sum(-1) / np.where(gm_, de_d, 0.0).sum(-1)
                            return de, adj_c, cov
                        de, ac, cov = std(S)
                        bde, bac, _ = std(SB)
                        qd, qr = q95(bde - bac), q95(ratio(bde, bac))
                        srows.append({**key, "level": lv, "filter": fname, "method": m, "target": t,
                                      "comparison_class": CLASSES[cj],
                                      "fraction_dark_exported_in_common_strata": float(de),
                                      "fraction_comparison_standardised": float(ac),
                                      "dark_exported_coverage_of_common_strata": float(cov),
                                      "difference": float(de - ac), "difference_ci_low": qd[0],
                                      "difference_ci_high": qd[1],
                                      "ratio_shared": float(ratio(de, ac)), "ratio_shared_ci_low": qr[0],
                                      "ratio_shared_ci_high": qr[1]})
    return rows, crows, srows


def split_half(fo, real_cols, real_names, mode, size, reps, seed, key):
    """EQ-internal baseline: hold out random EQ blocks as a pseudo-reference and compute the sharing
    of the remaining (focal) EQ families with it, next to their sharing with the real references."""
    rng = np.random.default_rng(seed)
    k = fo.nblocks
    gsize = np.bincount(fo.gb, minlength=k)
    names = ["eq_pseudo_reference"] + real_names
    gmin = [1, 2, 5, 10]
    out = np.empty((reps, len(gmin), 4, 1 + len(names)))
    meta = np.empty((reps, 2))
    for r in range(reps):
        perm = rng.permutation(k)
        if mode == "blocks":
            ref = perm[:size]
        else:
            ref = perm[:int(np.searchsorted(np.cumsum(gsize[perm]), size) + 1)]
        isref = np.zeros(k, bool)
        isref[ref] = True
        wr = isref[fo.gb].astype(float)
        wf = 1.0 - wr
        meta[r] = [isref.sum(), wr.sum()]
        n, d, p, cls, ng, nb, wgene = fo.fam_state(wf)
        inref = np.bincount(fo.f, wr[fo.g], fo.F) > 0
        for gi, gm_ in enumerate(gmin):
            sel = (n > 0) & (ng >= gm_)
            out[r, gi, :, 0] = np.bincount(cls[sel], minlength=4)
            out[r, gi, :, 1] = np.bincount(cls[sel], inref[sel].astype(float), 4)
            for j, col in enumerate(real_cols):
                out[r, gi, :, 2 + j] = np.bincount(cls[sel], fo.Sh[sel, col], 4)
    rows = []
    fr = ratio(out[..., 1:], out[..., :1])            # reps, gmin, class, names
    for gi, gm_ in enumerate(gmin):
        for j, nm in enumerate(names):
            for c in range(4):
                x = fr[:, gi, c, j]
                rows.append({**key, "split": f"{mode}_{size}", "reps": reps,
                             "pseudo_ref_blocks_mean": meta[:, 0].mean(),
                             "pseudo_ref_genomes_mean": meta[:, 1].mean(),
                             "filter": f"focal_genomes_ge{gm_}", "target": nm, "class": CLASSES[c],
                             "stat": "fraction_shared", "n_focal_families_mean": out[:, gi, c, 0].mean(),
                             "mean": np.nanmean(x), "p2.5": q95(x)[0], "p97.5": q95(x)[1]})
            for cj in (1, 2, 3):
                dd = fr[:, gi, 0, j] - fr[:, gi, cj, j]
                rr = ratio(fr[:, gi, 0, j], fr[:, gi, cj, j])
                for st, x in (("difference_dark_exported_minus_class", dd),
                              ("ratio_dark_exported_over_class", rr)):
                    rows.append({**key, "split": f"{mode}_{size}", "reps": reps,
                                 "pseudo_ref_blocks_mean": meta[:, 0].mean(),
                                 "pseudo_ref_genomes_mean": meta[:, 1].mean(),
                                 "filter": f"focal_genomes_ge{gm_}", "target": nm, "class": CLASSES[cj],
                                 "stat": st, "n_focal_families_mean": out[:, gi, cj, 0].mean(),
                                 "mean": np.nanmean(x), "p2.5": q95(x)[0], "p97.5": q95(x)[1]})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--families", nargs="+", required=True, help="label=path of MMseqs createtsv output")
    ap.add_argument("--search", type=Path, default=None, help="eq_vs_ext_by_group.tsv")
    ap.add_argument("--max-seqs", type=int, default=4000)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--split-reps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--expect-genes", type=int, default=3617011)
    ap.add_argument("--receipt-extra", type=str, default="{}")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(a.root / "analysis1_export"))
    import export_comparison as ref          # reference loader, bins, cohort names

    BINS = ref.BINS
    nl = len(BINS) - 1
    gm = pd.read_csv(a.root / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
    gm = gm[~gm.cohort.isin(ref.EXCLUDED)].copy()
    ani = pd.read_csv(a.root / "source-metadata/ani_groups/ani_95_groups_1957.tsv", sep="\t", dtype=str)
    gm = gm.merge(ani[["genome_id", "ani95_group_id"]], on="genome_id", how="left", validate="1:1")
    if gm.ani95_group_id.isna().any():
        raise SystemExit("genomes without ANI group")
    gm["cohort_short"] = gm.cohort.map(ref.SHORT)
    gm["block"] = np.where(gm.tax_genus.fillna("g__").isin(["g__", "", "missing", "none"]),
                           "ani::" + gm.ani95_group_id, gm.tax_genus)
    gm = gm.reset_index(drop=True)
    eq_orders = set(gm.loc[gm.cohort_short == "eq", "tax_order"]) - {"o__", "missing", "none"}

    df = ref.load(a.root)
    df = df[df.genome_id.isin(gm.genome_id)].reset_index(drop=True)
    checks = {"n_genes": int(len(df)), "n_genomes": int(len(gm)),
              "genes_match_expected": bool(len(df) == a.expect_genes)}
    if set(gm.genome_id) - set(df.genome_id):
        raise SystemExit("genomes without genes")
    extra = pd.read_csv(a.root / "mmseqs/harmonized_uniref50_per_gene.tsv.gz", sep="\t",
                        usecols=["gene_id", "direct_uniref50_max_pident"])
    df = df.merge(extra, on="gene_id", validate="1:1")
    gidx = pd.Index(gm.genome_id)
    df["gi"] = gidx.get_indexer(df.genome_id)
    assert (df.gi >= 0).all()
    gi = df.gi.to_numpy()

    c = gm.cohort_short
    so = gm.tax_order.isin(eq_orders)
    main3 = gm.tax_order.isin(eq_orders - {"o__RBG-16-68-12"})
    ext = c != "eq"
    unb = c.isin(["maeke35", "smag354", "spire1084"])
    tsets = {
        "sameorder_external247": (c == "external247") & so,
        "sameorder_smag": (c == "smag354") & so,
        "sameorder_spire": (c == "spire1084") & so,
        "sameorder_smag_spire": unb & so,
        "sameorder_smag_spire_excl_RBG-16-68-12": unb & main3,
        "sameorder_all_external": ext & so,
        "any_external": ext,
        "external247": c == "external247",
        "maeke35": c == "maeke35",
        "smag354": c == "smag354",
        "spire1084": c == "spire1084",
        "smag_spire_maeke": unb,
        "other_order_external": ext & ~so,
    }
    tnames = list(tsets)
    tmat = np.column_stack([tsets[t].to_numpy() for t in tnames])     # genomes x T
    ng_genes = np.bincount(gi, minlength=len(gm))
    pd.DataFrame([{"target": t, "n_genomes": int(tsets[t].sum()),
                   "n_blocks": int(gm.block[tsets[t]].nunique()),
                   "n_proteins": int(ng_genes[tsets[t].to_numpy()].sum()),
                   "cohorts": ",".join(sorted(set(c[tsets[t]])))} for t in tnames]).to_csv(
        a.out / "target_sets.tsv", sep="\t", index=False)

    # ---- circularity: how many proteins have a (near-)self hit in UniRef50 ----
    mp = df.direct_uniref50_max_pident.to_numpy()
    circ = []
    for scope, gmask in (("all_genomes", np.ones(len(gm), bool)), ("eq_named_orders", so.to_numpy())):
        for coh in ["eq", "external247", "maeke35", "smag354", "spire1084"]:
            sel = ((c == coh).to_numpy() & gmask)[gi]
            for gs, s2 in (("all_genes", sel), ("complete_interior", sel & df.interior.to_numpy())):
                if s2.sum() == 0:
                    continue
                x = mp[s2]
                circ.append({"scope": scope, "cohort": coh, "gene_set": gs, "n_proteins": int(s2.sum()),
                             "frac_dark": float(df.dark.to_numpy()[s2].mean()),
                             "frac_max_pident_ge90": float((x >= 90).mean()),
                             "frac_max_pident_ge95": float((x >= 95).mean()),
                             "frac_max_pident_ge99": float((x >= 99).mean()),
                             "frac_max_pident_eq100": float((x >= 100).mean())})
    circ = pd.DataFrame(circ)
    # Coverage-aware version from the raw hit table (identity >= 95 % over >= 80 % of the query).
    cov_status = "not_run"
    try:
        man = pd.read_csv(a.root / "genecall/harmonized_gene_manifest.tsv.gz", sep="\t",
                          usecols=["gene_id", "canonical_seq_sha256", "aa_length"])
        man = man[man.gene_id.isin(df.gene_id)]
        hlen = man.drop_duplicates("canonical_seq_sha256").set_index("canonical_seq_sha256").aa_length
        glen = man.set_index("gene_id").aa_length
        best = {}
        nrow = 0
        for ch in pd.read_csv(a.root / "mmseqs/harmonized_vs_uniref50_s75_e1e3.m8.gz", sep="\t", header=None,
                              usecols=[0, 2, 3], names=["q", "pid", "aln"], chunksize=5_000_000):
            nrow += len(ch)
            ch = ch[ch.pid >= 95]
            ln = ch.q.map(hlen)
            if ln.isna().any():
                ln = ln.fillna(ch.q.map(glen))
            ch = ch[(ch.aln / ln) >= 0.8]
            best.update(dict.fromkeys(ch.q.unique(), True))
        keys = pd.Index(list(best))
        selfhit = man.canonical_seq_sha256.isin(keys) | man.gene_id.isin(keys)
        sh = pd.Series(selfhit.to_numpy(), index=man.gene_id).reindex(df.gene_id).to_numpy()
        df["selfhit"] = sh
        add = []
        for _, r in circ.iterrows():
            gmask = np.ones(len(gm), bool) if r.scope == "all_genomes" else so.to_numpy()
            s2 = ((c == r.cohort).to_numpy() & gmask)[gi]
            if r.gene_set == "complete_interior":
                s2 = s2 & df.interior.to_numpy()
            add.append(float(sh[s2].mean()))
        circ["frac_pident_ge95_alncov_ge80"] = add
        cov_status = f"ok ({nrow} hit rows read)"
    except Exception as e:  # the raw table is optional
        cov_status = f"failed: {e!r}"
        df["selfhit"] = mp >= 95
    circ.to_csv(a.out / "circularity_uniref50_selfhits.tsv", sep="\t", index=False)
    print(circ.to_string(), flush=True)

    # ---- search hits per EQ protein and target set ----
    eq_g = (c == "eq").to_numpy()
    is_eq = eq_g[gi]
    eq_idx = np.flatnonzero(is_eq)
    hits_eq, search_info = None, {"status": "not_provided"}
    if a.search is not None:
        s = pd.read_csv(a.search, sep="\t", header=None,
                        names=["query", "group", "n", "bits", "evalue", "pident", "qcov", "tcov", "best_target", "tot"])
        grp = (gm.cohort_short + "|" + gm.tax_order)
        pos = pd.Index(df.gene_id.to_numpy()[eq_idx]).get_indexer(s["query"])
        if (pos < 0).any():
            raise SystemExit("search queries that are not EQ genes")
        if not ((s.qcov >= 0.5 - 1e-6).all() and (s.evalue <= 1e-5).all()):
            raise SystemExit("search rows outside thresholds")
        hits_eq = np.zeros((len(eq_idx), len(tnames)), bool)
        for ti, t in enumerate(tnames):
            gset = set(grp[tsets[t]])
            other = set(grp[~tsets[t]])
            if gset & other:
                raise SystemExit(f"target {t} is not a union of cohort|order groups")
            hits_eq[pos[s.group.isin(gset).to_numpy()], ti] = True
        tot = s.drop_duplicates("query").tot
        search_info = {"status": "loaded", "rows": int(len(s)), "queries_with_hit": int(s["query"].nunique()),
                       "eq_queries": int(len(eq_idx)),
                       "queries_at_max_seqs_cap": int((tot >= a.max_seqs).sum()),
                       "file_sha256": sha256(a.search)}
        del s
    print("search", search_info, flush=True)

    block_codes, block_names = pd.factorize(gm.block)
    eq_orders_l = sorted(eq_orders)
    all_rows, all_crows, all_srows, split_rows, count_rows, order_rows, ext_rows = [], [], [], [], [], [], []
    fam_info = {}
    for spec in a.families:
        label, path = spec.split("=", 1)
        ft = pd.read_csv(path, sep="\t", header=None, names=["rep", "gene_id"], dtype=str)
        ok_unique = bool(ft.gene_id.is_unique)
        ok_set = bool(len(ft) == len(df) and set(ft.gene_id) == set(df.gene_id))
        fam_info[label] = {"path": path, "sha256": sha256(path), "n_rows": int(len(ft)),
                           "n_families": int(ft.rep.nunique()),
                           "every_gene_in_exactly_one_family": ok_unique and ok_set}
        print(label, fam_info[label], flush=True)
        if not (ok_unique and ok_set):
            raise SystemExit(f"identifier closure failed for {label}")
        rep_codes, rep_names = pd.factorize(ft.rep)
        fam = pd.Series(rep_codes, index=ft.gene_id).reindex(df.gene_id).to_numpy()
        nfam = len(rep_names)
        fam_info[label]["n_singleton_families"] = int((np.bincount(fam, minlength=nfam) == 1).sum())
        del ft
        # family x target membership
        gene_t = tmat[gi]                                    # genes x T
        Sh = np.zeros((nfam, len(tnames)), bool)
        for ti in range(len(tnames)):
            Sh[np.unique(fam[gene_t[:, ti]]), ti] = True
        eq_fam = np.zeros(nfam, bool)
        eq_fam[np.unique(fam[is_eq])] = True
        ShEQ = eq_fam[:, None]

        for gene_set in ("all_genes", "complete_interior"):
            key = {"identity": label, "gene_set": gene_set}
            gsel = is_eq if gene_set == "all_genes" else (is_eq & df.interior.to_numpy())
            sub = np.flatnonzero(gsel)
            # focal genome index among EQ genomes
            eq_gen = np.flatnonzero(eq_g)
            remap = -np.ones(len(gm), int)
            remap[eq_gen] = np.arange(len(eq_gen))
            bcodes, _ = pd.factorize(block_codes[eq_gen])
            h = hits_eq[np.searchsorted(eq_idx, sub)] if hits_eq is not None else None
            fo = Focal(fam[sub], remap[gi[sub]], df.dark.to_numpy()[sub], df.sp.to_numpy()[sub],
                       df.aa_length.to_numpy()[sub], bcodes, Sh, h, BINS)
            point, boot = fo.run(a.draws, a.seed)
            flt = filters(nl, BINS)
            methods = METHODS if h is not None else ["cluster"]
            r, cr, sr = summarise(point, boot, flt, tnames, {**key, "focal": "eq"}, methods, nl)
            all_rows += r; all_crows += cr; all_srows += sr
            del boot
            print("done forward", key, flush=True)

            full = fo.full
            cls = full["cls"]
            # counts per class
            for ci, cn in enumerate(CLASSES):
                m = cls == ci
                count_rows.append({**key, "class": cn, "n_families": int(m.sum()),
                                   "n_eq_proteins_in_families_of_class": int(full["n"][m].sum()),
                                   "n_eq_proteins_own_class": int((fo.own_cls == ci).sum()),
                                   "n_families_single_eq_protein": int((m & (full["n"] == 1)).sum()),
                                   "n_families_ge2_eq_genomes": int((m & (full["ng"] >= 2)).sum()),
                                   "n_families_ge10genomes_ge3blocks": int((m & (full["ng"] >= 10) & (full["nb"] >= 3)).sum())})
            # split-half baseline (EQ-internal)
            real = ["sameorder_external247", "sameorder_smag_spire", "sameorder_all_external", "any_external"]
            rc = [tnames.index(x) for x in real]
            nb67 = int(gm.block[tsets["sameorder_external247"]].nunique())
            n67 = int(tsets["sameorder_external247"].sum())
            split_rows += split_half(fo, rc, real, "blocks", nb67, a.split_reps, a.seed + 1, key)
            split_rows += split_half(fo, rc, real, "genomes", n67, a.split_reps, a.seed + 2, key)
            print("done split", key, flush=True)

            # order distribution and family table
            ords = gm.tax_order.to_numpy()[gi[sub]]
            ocode = np.array([eq_orders_l.index(o) if o in eq_orders else len(eq_orders_l) for o in ords])
            O = np.bincount(fo.f * (len(eq_orders_l) + 1) + ocode,
                            minlength=fo.F * (len(eq_orders_l) + 1)).reshape(fo.F, -1)
            n_ord = (O[:, :len(eq_orders_l)] > 0).sum(1)
            Hs = fo.search_shared(np.ones(len(sub)))
            ti_any = tnames.index("any_external")
            restricted_cluster = fo.Sh[:, ti_any] == 0
            restricted_either = restricted_cluster & (Hs[:, ti_any] == 0)
            wide = (full["ng"] >= 10) & (full["nb"] >= 3)
            for ci, cn in enumerate(CLASSES):
                m = cls == ci
                for oi, on in enumerate(eq_orders_l + ["unassigned_order"]):
                    po = m & (O[:, oi] > 0)
                    order_rows.append({**key, "class": cn, "table": "families_with_members_in_order", "level": on,
                                       "n_families": int(po.sum()), "n_eq_proteins_in_order": int(O[m, oi].sum()),
                                       "n_eq_restricted_cluster": int((po & restricted_cluster).sum()),
                                       "n_eq_restricted_cluster_and_search": int((po & restricted_either).sum()),
                                       "n_widespread": int((po & wide).sum()),
                                       "n_widespread_eq_restricted_cluster": int((po & wide & restricted_cluster).sum()),
                                       "n_widespread_eq_restricted_cluster_and_search": int((po & wide & restricted_either).sum())})
                for k_ in range(0, len(eq_orders_l) + 1):
                    po = m & (n_ord == k_)
                    order_rows.append({**key, "class": cn, "table": "n_named_orders_spanned", "level": str(k_),
                                       "n_families": int(po.sum()), "n_eq_proteins_in_order": int(full["n"][po].sum()),
                                       "n_eq_restricted_cluster": int((po & restricted_cluster).sum()),
                                       "n_eq_restricted_cluster_and_search": int((po & restricted_either).sum()),
                                       "n_widespread": int((po & wide).sum()),
                                       "n_widespread_eq_restricted_cluster": int((po & wide & restricted_cluster).sum()),
                                       "n_widespread_eq_restricted_cluster_and_search": int((po & wide & restricted_either).sum())})
                order_rows.append({**key, "class": cn, "table": "total", "level": "all",
                                   "n_families": int(m.sum()), "n_eq_proteins_in_order": int(full["n"][m].sum()),
                                   "n_eq_restricted_cluster": int((m & restricted_cluster).sum()),
                                   "n_eq_restricted_cluster_and_search": int((m & restricted_either).sum()),
                                   "n_widespread": int((m & wide).sum()),
                                   "n_widespread_eq_restricted_cluster": int((m & wide & restricted_cluster).sum()),
                                   "n_widespread_eq_restricted_cluster_and_search": int((m & wide & restricted_either).sum())})
            ftab = pd.DataFrame({"family_rep": rep_names[fo.fam_global], "class": np.array(CLASSES)[cls],
                                 "n_eq_proteins": full["n"].astype(int), "n_eq_dark": full["d"].astype(int),
                                 "n_eq_signal_peptide": full["p"].astype(int),
                                 "n_eq_genomes": full["ng"].astype(int), "n_eq_blocks": full["nb"].astype(int),
                                 "median_eq_aa_length": fo.med_len})
            for oi, on in enumerate(eq_orders_l + ["unassigned_order"]):
                ftab["n_eq_" + on] = O[:, oi]
            for ti, t in enumerate(tnames):
                ftab["cluster_member_in_" + t] = fo.Sh[:, ti].astype(int)
                if h is not None:
                    ftab["search_hit_in_" + t] = Hs[:, ti].astype(int)
            ftab.to_csv(a.out / f"eq_family_table_{label}_{gene_set}.tsv.gz", sep="\t", index=False)
            if gene_set == "all_genes":
                ftab[(ftab["class"] == "dark_exported") & (ftab.n_eq_genomes >= 10) & (ftab.n_eq_blocks >= 3)
                     ].sort_values("n_eq_genomes", ascending=False).to_csv(
                    a.out / f"dark_exported_widespread_families_{label}.tsv", sep="\t", index=False)

            # properties of the external members of shared families (all genes only)
            if gene_set == "all_genes":
                cls_g = -np.ones(nfam, int)
                cls_g[fo.fam_global] = cls
                dk, spv, shv = df.dark.to_numpy(), df.sp.to_numpy(), df.selfhit.to_numpy().astype(bool)
                for ti, t in enumerate(tnames):
                    em = gene_t[:, ti] & (cls_g[fam] >= 0)
                    ef = fam[em]
                    n_e = np.bincount(ef, minlength=nfam)
                    d_e = np.bincount(ef, dk[em].astype(float), nfam)
                    p_e = np.bincount(ef, spv[em].astype(float), nfam)
                    for ci, cn in enumerate(CLASSES):
                        fm = (cls_g == ci) & (n_e > 0)
                        pm = em & (cls_g[fam] == ci)
                        ext_rows.append({**key, "target": t, "eq_family_class": cn,
                                         "n_shared_families": int(fm.sum()),
                                         "n_external_members": int(pm.sum()),
                                         "frac_external_members_dark": float(dk[pm].mean()) if pm.any() else np.nan,
                                         "frac_external_members_signal_peptide": float(spv[pm].mean()) if pm.any() else np.nan,
                                         "frac_external_members_uniref50_selfhit": float(shv[pm].mean()) if pm.any() else np.nan,
                                         "frac_families_external_majority_dark": float((2 * d_e[fm] > n_e[fm]).mean()) if fm.any() else np.nan,
                                         "frac_families_external_majority_signal_peptide": float((2 * p_e[fm] > n_e[fm]).mean()) if fm.any() else np.nan,
                                         "frac_families_external_majority_dark_and_exported": float(((2 * d_e[fm] > n_e[fm]) & (2 * p_e[fm] > n_e[fm])).mean()) if fm.any() else np.nan})

            # reverse direction: external same-order families found in EQ (cluster membership only)
            for rname in ("sameorder_external247", "sameorder_smag_spire",
                          "sameorder_smag_spire_excl_RBG-16-68-12", "sameorder_all_external"):
                rg = tsets[rname].to_numpy()
                rsel = rg[gi] if gene_set == "all_genes" else (rg[gi] & df.interior.to_numpy())
                rsub = np.flatnonzero(rsel)
                rgen = np.flatnonzero(rg)
                rmap = -np.ones(len(gm), int)
                rmap[rgen] = np.arange(len(rgen))
                rb, _ = pd.factorize(block_codes[rgen])
                # the EQ side keeps all EQ genes, as the external side does in the forward direction
                ro = Focal(fam[rsub], rmap[gi[rsub]], df.dark.to_numpy()[rsub], df.sp.to_numpy()[rsub],
                           df.aa_length.to_numpy()[rsub], rb, ShEQ, None, BINS)
                pt, bt = ro.run(a.draws, a.seed)
                r, cr, sr = summarise(pt, bt, flt, ["eq_all"], {**key, "focal": rname}, ["cluster"], nl)
                all_rows += r; all_crows += cr; all_srows += sr
                del bt
            print("done reverse", key, flush=True)

    def w(rows, name):
        pd.DataFrame(rows).to_csv(a.out / name, sep="\t", index=False)
    w(all_rows, "sharing_fractions.tsv")
    w(all_crows, "sharing_contrasts.tsv")
    w(all_srows, "sharing_contrasts_standardised_length_prevalence.tsv")
    w(split_rows, "eq_internal_split_baseline.tsv")
    w(count_rows, "class_counts.tsv")
    w(order_rows, "order_distribution.tsv")
    w(ext_rows, "shared_family_external_member_properties.tsv")

    receipt = {
        "status": "COMPLETE", "seed": a.seed, "draws": a.draws, "split_reps": a.split_reps,
        "block": "named genus, else ANI95 group", "eq_named_orders": eq_orders_l,
        "length_bins": BINS[:-1], "focal_genome_bins": [1] + GB_EDGES, "focal_block_bins": [1] + BB_EDGES,
        "checks": checks, "families": fam_info, "search": search_info,
        "uniref50_selfhit_coverage_table": cov_status,
        "inputs": {p: sha256(a.root / p) for p in [
            "metadata/harmonized_genome_manifest.tsv",
            "source-metadata/ani_groups/ani_95_groups_1957.tsv",
            "genecall/harmonized_gene_manifest.tsv.gz",
            "mmseqs/harmonized_uniref50_per_gene.tsv.gz",
            "signalp/harmonized_signalp6_predictions.tsv.gz",
            "analysis1_export/export_comparison.py"]},
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "slurm_job_id_analysis": os.environ.get("SLURM_JOB_ID"),
        **json.loads(a.receipt_extra),
    }
    (a.out / "receipt.json").write_text(json.dumps(receipt, indent=1))
    print("ANALYSIS_DONE", flush=True)


if __name__ == "__main__":
    main()
