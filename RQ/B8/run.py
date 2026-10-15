#!/usr/bin/env python3
"""
RQ-B8: Is the dark proteome temporally stable across the 5 expeditions, or does it
turn over between trips?

Run from repo root:  python3 RQ/B8/run.py

Inputs  (results/darkmatter_ecology/):
  env_matrix.tsv         sample -> site, compartment, trip, trip_num, lat, lon, ...
  darkfam_abundance.tsv  family_id + 551 per-sample CoverM abundance columns (~92 MB, streamed)
  darkfam_meta.tsv       family_id, n_members, n_samples_catalog (not required here)

Outputs (RQ/B8/):
  fig_temporal.pdf / .png
  turnover_by_class.tsv
  persistence_histogram.tsv
  permanova.tsv
  summary.json

Dependencies: stdlib + numpy + scipy + matplotlib
"""
import os, sys, csv, json
from collections import defaultdict, Counter
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ECO  = os.path.join(ROOT, "results", "darkmatter_ecology")
ENVF = os.path.join(ECO, "env_matrix.tsv")
ABUF = os.path.join(ECO, "darkfam_abundance.tsv")

RNG = np.random.default_rng(42)
N_PERM = 999

# ----------------------------------------------------------------------------- env
def load_env():
    env = {}
    with open(ENVF) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            env[row["sample"]] = {
                "site": row["site"],
                "compartment": row["compartment"],
                "trip": int(row["trip_num"]),
            }
    return env

# ------------------------------------------------------------------- stream matrix
def load_abundance(samples_keep):
    """Stream the big table; return (family_ids, matrix[nfam, nsamp], col_samples)."""
    with open(ABUF) as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = header[1:]
        keep_idx = [i for i, c in enumerate(cols) if c in samples_keep]
        col_samples = [cols[i] for i in keep_idx]
        fam_ids = []
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            fam_ids.append(parts[0])
            vals = parts[1:]
            rows.append([float(vals[i]) for i in keep_idx])
    mat = np.asarray(rows, dtype=np.float64)          # [nfam, nsamp]
    return fam_ids, mat, col_samples

# --------------------------------------------------------------------- PERMANOVA
def permanova(D, labels, n_perm=N_PERM, rng=RNG):
    """One-way PERMANOVA (Anderson 2001) on a squareform distance matrix D."""
    labels = np.asarray(labels)
    n = D.shape[0]
    groups = np.unique(labels)
    a = len(groups)
    D2 = D ** 2
    ss_total = D2[np.triu_indices(n, 1)].sum() / n

    def ss_within(lab):
        ssw = 0.0
        for g in groups:
            idx = np.where(lab == g)[0]
            ng = len(idx)
            if ng > 1:
                sub = D2[np.ix_(idx, idx)]
                ssw += sub[np.triu_indices(ng, 1)].sum() / ng
        return ssw

    ssw_obs = ss_within(labels)
    ssa_obs = ss_total - ssw_obs
    F_obs = (ssa_obs / (a - 1)) / (ssw_obs / (n - a))

    ge = 1  # count observed
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        ssw = ss_within(perm)
        ssa = ss_total - ssw
        F = (ssa / (a - 1)) / (ssw / (n - a))
        if F >= F_obs:
            ge += 1
    p = ge / (n_perm + 1)
    R2 = ssa_obs / ss_total
    return F_obs, R2, p

# ------------------------------------------------------------------------- main
def main():
    env = load_env()
    print(f"[env] {len(env)} samples")
    by_trip = defaultdict(lambda: [0, set()])
    for s, m in env.items():
        by_trip[m["trip"]][0] += 1
        by_trip[m["trip"]][1].add(m["site"])
    print("[env] samples / sites per trip:")
    for t in sorted(by_trip):
        print(f"   trip{t}: {by_trip[t][0]:4d} samples  {len(by_trip[t][1]):3d} sites")

    fam_ids, mat, cols = load_abundance(set(env))
    nfam, nsamp = mat.shape
    print(f"[abund] {nfam} dark families x {nsamp} samples")
    meta = [env[c] for c in cols]
    trip = np.array([m["trip"] for m in meta])
    site = np.array([m["site"] for m in meta])

    # ---- temporal persistence: trips in which each family is detected (>0) ----
    pres = mat > 0
    trips_present = np.zeros(nfam, dtype=int)
    detected_in = defaultdict(set)
    for j in range(nsamp):
        t = trip[j]
        col = pres[:, j]
        for i in np.where(col)[0]:
            detected_in[i].add(t)
    for i in range(nfam):
        trips_present[i] = len(detected_in.get(i, ()))
    # families never detected in any kept sample -> 0 trips
    hist = Counter(trips_present.tolist())
    n_detected = int((trips_present > 0).sum())
    n_all5 = int((trips_present == 5).sum())
    n_single = int((trips_present == 1).sum())
    print(f"[persist] detected families: {n_detected}")
    print(f"[persist] in all 5 trips (core): {n_all5} "
          f"({100*n_all5/n_detected:.1f}% of detected)")
    print(f"[persist] trip-restricted (1 trip): {n_single} "
          f"({100*n_single/n_detected:.1f}% of detected)")

    with open(os.path.join(HERE, "persistence_histogram.tsv"), "w") as fh:
        fh.write("n_trips\tn_families\tpct_of_detected\n")
        for k in range(0, 6):
            v = hist.get(k, 0)
            pct = 100 * v / n_detected if k > 0 else float("nan")
            fh.write(f"{k}\t{v}\t{pct:.3f}\n")

    # ---- community turnover: Bray-Curtis on depth-normalised sample vectors ----
    col_sum = mat.sum(axis=0)
    col_sum[col_sum == 0] = 1.0
    rel = (mat / col_sum).T                 # [nsamp, nfam], rows sum to 1
    print("[turnover] computing Bray-Curtis pairwise ...")
    bc = pdist(rel, metric="braycurtis")
    D = squareform(bc)

    # classify each unordered sample pair
    iu = np.triu_indices(nsamp, 1)
    same_site = site[iu[0]] == site[iu[1]]
    same_trip = trip[iu[0]] == trip[iu[1]]
    dvals = bc
    classes = {
        "within-site_within-trip (replicates)": dvals[same_site & same_trip],
        "within-site_across-trip (temporal)":   dvals[same_site & ~same_trip],
        "across-site_within-trip (spatial)":    dvals[~same_site & same_trip],
        "across-site_across-trip":              dvals[~same_site & ~same_trip],
    }
    with open(os.path.join(HERE, "turnover_by_class.tsv"), "w") as fh:
        fh.write("class\tn_pairs\tmean_BC\tmedian_BC\tsd_BC\n")
        for k, v in classes.items():
            fh.write(f"{k}\t{len(v)}\t{np.mean(v):.4f}\t{np.median(v):.4f}\t{np.std(v):.4f}\n")
    print("[turnover] mean Bray-Curtis by class:")
    for k, v in classes.items():
        print(f"   {k:42s} n={len(v):7d}  mean={np.mean(v):.3f}  median={np.median(v):.3f}")

    # Key test: are repeat-site pairs over time MORE similar than different-site pairs?
    temporal = classes["within-site_across-trip (temporal)"]
    spatial_t = classes["across-site_across-trip"]
    u, p_ws = mannwhitneyu(temporal, spatial_t, alternative="less")
    print(f"[test] within-site-across-trip vs across-site-across-trip: "
          f"MWU p(less)={p_ws:.2e}  (mean {np.mean(temporal):.3f} vs {np.mean(spatial_t):.3f})")

    # ---- PERMANOVA: trip vs site as grouping factors ----
    print(f"[permanova] {N_PERM} permutations ...")
    F_trip, R2_trip, p_trip = permanova(D, trip)
    F_site, R2_site, p_site = permanova(D, site)
    print(f"[permanova] trip:  pseudo-F={F_trip:.2f}  R2={R2_trip:.4f}  p={p_trip:.4f}")
    print(f"[permanova] site:  pseudo-F={F_site:.2f}  R2={R2_site:.4f}  p={p_site:.4f}")

    with open(os.path.join(HERE, "permanova.tsv"), "w") as fh:
        fh.write("factor\tn_groups\tpseudo_F\tR2\tp_value\tn_perm\n")
        fh.write(f"trip\t{len(np.unique(trip))}\t{F_trip:.4f}\t{R2_trip:.5f}\t{p_trip:.4f}\t{N_PERM}\n")
        fh.write(f"site\t{len(np.unique(site))}\t{F_site:.4f}\t{R2_site:.5f}\t{p_site:.4f}\t{N_PERM}\n")

    # ----------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    order = [
        ("within-site_within-trip (replicates)", "repeat\n(same site,\nsame trip)"),
        ("within-site_across-trip (temporal)",   "temporal\n(same site,\ndiff trip)"),
        ("across-site_within-trip (spatial)",    "spatial\n(diff site,\nsame trip)"),
        ("across-site_across-trip",              "diff site,\ndiff trip"),
    ]
    data = [classes[k] for k, _ in order]
    labels = [lab for _, lab in order]
    cols_box = ["#4c72b0", "#55a868", "#c44e52", "#8172b3"]
    bp = ax1.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], cols_box):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax1.set_ylabel("Bray-Curtis dissimilarity (dark families)")
    ax1.set_title("Community turnover: time vs space")
    ax1.tick_params(axis="x", labelsize=8)
    ax1.text(0.02, 0.02,
             f"PERMANOVA trip R2={R2_trip:.3f} (p={p_trip:.3f})\n"
             f"PERMANOVA site R2={R2_site:.3f} (p={p_site:.3f})",
             transform=ax1.transAxes, fontsize=8, va="bottom",
             bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.8))

    ks = list(range(1, 6))
    vs = [hist.get(k, 0) for k in ks]
    bars = ax2.bar(ks, vs, color="#4c72b0", alpha=0.85)
    ax2.set_xlabel("Number of expeditions in which family is detected")
    ax2.set_ylabel("Number of dark families")
    ax2.set_title("Temporal persistence of dark-protein families")
    ax2.set_xticks(ks)
    for k, v in zip(ks, vs):
        ax2.text(k, v, f"{100*v/n_detected:.0f}%", ha="center", va="bottom", fontsize=8)
    ax2.text(0.97, 0.97,
             f"core (all 5 trips): {n_all5:,}\n({100*n_all5/n_detected:.1f}% of {n_detected:,})",
             transform=ax2.transAxes, fontsize=9, ha="right", va="top",
             bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.8))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_temporal.pdf"))
    fig.savefig(os.path.join(HERE, "fig_temporal.png"), dpi=150)
    print("[fig] wrote fig_temporal.pdf / .png")

    # ----------------------------------------------------------------- summary
    summary = {
        "n_samples": int(nsamp),
        "n_dark_families_total": int(nfam),
        "n_dark_families_detected": n_detected,
        "samples_per_trip": {int(t): by_trip[t][0] for t in sorted(by_trip)},
        "sites_per_trip": {int(t): len(by_trip[t][1]) for t in sorted(by_trip)},
        "persistence": {
            "core_all5_trips": n_all5,
            "core_all5_pct_of_detected": round(100 * n_all5 / n_detected, 2),
            "trip_restricted_1trip": n_single,
            "trip_restricted_pct": round(100 * n_single / n_detected, 2),
            "histogram_n_families_by_n_trips": {int(k): int(hist.get(k, 0)) for k in range(0, 6)},
        },
        "turnover_mean_BC": {k: round(float(np.mean(v)), 4) for k, v in classes.items()},
        "test_temporal_vs_spatial_MWU_p_less": float(p_ws),
        "permanova": {
            "trip": {"pseudo_F": round(F_trip, 3), "R2": round(R2_trip, 4), "p": p_trip},
            "site": {"pseudo_F": round(F_site, 3), "R2": round(R2_site, 4), "p": p_site},
        },
    }
    with open(os.path.join(HERE, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("[done] summary.json written")

if __name__ == "__main__":
    main()
