#!/usr/bin/env python3
"""RQ-B7. Dispersal limitation vs environmental filtering for dark-gene biogeography.

How much of dark-gene community structure is pure-space (dispersal limitation) vs
pure-environment (environmental filtering) vs shared?

This script:
  1. Re-extracts the per-block PERMANOVA (adonis) R2 / p-values and the distance-decay
     Mantel r from results/darkmatter_ecology/phase_b_summary.txt (computed by
     code/phase_b_ecology.py) into a clean table.
  2. Recomputes a proper Borcard distance-based variance partition (pure space | pure
     environment | shared | residual) directly on the real dark-family abundance matrix,
     using the self-contained PERMANOVA in code/eco_stats.py (no scikit-bio / vegan).
     Two environment blocks are partitioned against space (lat/lon):
       - CLIMATE  (temp, precip, aridity)  on the climate-complete subset
       - GEOCHEM  (XRF elements)           on the XRF subset

Run from repo root:   python3 RQ/B7/run.py
Outputs (this dir):   b7_blocks.tsv, b7_varpart.tsv, fig_varpart.pdf/.png
"""
import os, sys, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
import phase_b_ecology as pb          # noqa: E402
import eco_stats as es                # noqa: E402

ECO = os.path.join(ROOT, "results", "darkmatter_ecology")
ENVF = os.path.join(ECO, "env_matrix.tsv")
ABUND = os.path.join(ECO, "darkfam_abundance.tsv")
SUMMARY = os.path.join(ECO, "phase_b_summary.txt")
PERMS = 499


# ---------------------------------------------------------------- 1. parse summary
def parse_summary(path):
    """Pull the per-block adonis rows and the Mantel distance-decay r from phase_b_summary.txt."""
    blocks, mantel_r, mantel_p = [], None, None
    in_blocks = False
    for ln in open(path):
        s = ln.rstrip("\n")
        if s.startswith("=== PERMANOVA"):
            in_blocks = True; continue
        if s.startswith("=== variance") or s.startswith("=== distance"):
            in_blocks = False
        if in_blocks:
            m = re.match(r"^(\S.*?)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", s)
            if m and m.group(1) != "block":
                blocks.append((m.group(1), int(m.group(2)),
                               float(m.group(3)), float(m.group(4)), float(m.group(5))))
        m2 = re.search(r"r=([+\-][\d.]+)\s+p=([\d.]+)", s)
        if m2:
            mantel_r, mantel_p = float(m2.group(1)), float(m2.group(2))
    return blocks, mantel_r, mantel_p


# ---------------------------------------------------------------- 2. Borcard partition
def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def zscore(rows):
    X = np.array(rows, float)
    mu, sd = np.nanmean(X, 0), np.nanstd(X, 0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def borcard(D, Xspace, Xenv, perms=PERMS):
    """2-block distance-based variance partition (Borcard et al. 1992) via PERMANOVA R2.
    ab = R2(space alone); bc = R2(env alone); abc = R2(space+env combined).
    pure_space = abc-bc ; pure_env = abc-ab ; shared = ab+bc-abc ; residual = 1-abc."""
    ab = es.adonis(D, Xspace, perms=perms)["R2_marginal"]
    bc = es.adonis(D, Xenv, perms=perms)["R2_marginal"]
    comb = es.adonis(D, np.hstack([Xspace, Xenv]), perms=perms)
    abc = comb["R2_marginal"]
    # significance of each unique fraction = partial adonis (term | other block)
    p_space = es.adonis(D, Xspace, covariate=Xenv, perms=perms)["p"]
    p_env = es.adonis(D, Xenv, covariate=Xspace, perms=perms)["p"]
    return {
        "space_alone": ab, "env_alone": bc, "combined": abc,
        "pure_space": abc - bc, "pure_env": abc - ab,
        "shared": ab + bc - abc, "residual": 1 - abc,
        "p_pure_space": p_space, "p_pure_env": p_env,
    }


def main():
    print("=== RQ-B7: dispersal limitation vs environmental filtering ===\n")

    # ---- 1. tabulate phase_b numbers ----
    blocks, mantel_r, mantel_p = parse_summary(SUMMARY)
    print("--- PERMANOVA (adonis) per environmental block  [from phase_b_summary.txt] ---")
    print(f"{'block':<20}{'n':>5}{'R2':>8}{'pseudoF':>9}{'p':>8}")
    bt = open(os.path.join(HERE, "b7_blocks.tsv"), "w")
    bt.write("block\tn\tR2\tpseudoF\tp\n")
    for name, n, r2, f, p in blocks:
        print(f"{name:<20}{n:>5}{r2:>8.3f}{f:>9.2f}{p:>8.3f}")
        bt.write(f"{name}\t{n}\t{r2}\t{f}\t{p}\n")
    bt.close()
    print(f"\nDistance-decay Mantel (community vs geography): r={mantel_r:+.3f} p={mantel_p:.3f}")
    print("  -> significant positive distance-decay = dispersal limitation present.\n")

    # ---- 2. recompute clean Borcard space/env/shared partition on the real matrix ----
    env = pb.load_env(ENVF)
    samples, fams, M = pb.load_matrix(ABUND)
    keep = [i for i, s in enumerate(samples) if s in env]
    samples = [samples[i] for i in keep]; M = M[keep]
    prev = (M > 0).sum(0); M = M[:, prev >= pb.MIN_PREV]
    print(f"recomputing on real matrix: {M.shape[0]} samples x {M.shape[1]} families "
          f"(prev>={pb.MIN_PREV})")

    lat = np.array([fnum(env[s].get("lat")) for s in samples], float)
    lon = np.array([fnum(env[s].get("lon")) for s in samples], float)

    vt = open(os.path.join(HERE, "b7_varpart.tsv"), "w")
    vt.write("env_block\tn_samples\tpure_space\tpure_env\tshared\tresidual\t"
             "p_pure_space\tp_pure_env\tspace_alone\tenv_alone\tcombined\n")

    partitions = {}
    specs = {
        "climate": (["temp", "precip", "aridity_index"],
                    lambda s: fnum(env[s].get("temp")) is not None
                    and fnum(env[s].get("precip")) is not None),
        "geochem": (["V", "Fe", "Ca", "K", "Si", "Ti", "Mn", "salinity_proxy"],
                    lambda s: env[s].get("has_xrf") == "1"),
    }
    for label, (cols, sel) in specs.items():
        idx = [i for i, s in enumerate(samples) if sel(s)]
        if len(idx) < 30:
            continue
        sub = M[idx]
        D = es.bray_curtis(es.hellinger(sub))
        Xspace = zscore([[lat[i], lon[i]] for i in idx])
        Xenv = zscore([[fnum(env[samples[i]].get(c) or 0) for c in cols] for i in idx])
        part = borcard(D, Xspace, Xenv)
        partitions[label] = (len(idx), part)
        print(f"\n--- Borcard variance partition: SPACE vs {label.upper()} "
              f"(n={len(idx)}) ---")
        print(f"  pure space  (dispersal)    R2={part['pure_space']:.3f}  "
              f"(p={part['p_pure_space']:.3f})")
        print(f"  pure {label:<8} (filtering) R2={part['pure_env']:.3f}  "
              f"(p={part['p_pure_env']:.3f})")
        print(f"  shared space+{label:<8}      R2={part['shared']:.3f}")
        print(f"  residual (unexplained)     R2={part['residual']:.3f}")
        vt.write(f"{label}\t{len(idx)}\t{part['pure_space']:.4f}\t{part['pure_env']:.4f}\t"
                 f"{part['shared']:.4f}\t{part['residual']:.4f}\t{part['p_pure_space']:.4f}\t"
                 f"{part['p_pure_env']:.4f}\t{part['space_alone']:.4f}\t"
                 f"{part['env_alone']:.4f}\t{part['combined']:.4f}\n")
    vt.close()

    make_figure(partitions, mantel_r)
    print("\nwritten: b7_blocks.tsv, b7_varpart.tsv, fig_varpart.pdf/.png")
    return partitions, blocks, mantel_r


def make_figure(partitions, mantel_r):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[fig] matplotlib unavailable: {e}")
        return
    if not partitions:
        return
    labels = list(partitions)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    comps = ["pure_space", "shared", "pure_env"]
    colors = {"pure_space": "#3b6ea5", "shared": "#bdbdbd", "pure_env": "#c0392b"}
    names = {"pure_space": "pure space (dispersal)", "shared": "shared",
             "pure_env": "pure environment (filtering)"}
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for c in comps:
        vals = np.array([max(partitions[l][1][c], 0) for l in labels])
        ax.bar(x, vals, bottom=bottom, color=colors[c], label=names[c], width=0.55)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([f"space vs {l}\n(n={partitions[l][0]})" for l in labels])
    ax.set_ylabel("explained R$^2$ (PERMANOVA)")
    ax.set_title("Dark-gene biogeography: dispersal vs environmental filtering")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.text(0.01, 0.97, f"distance-decay Mantel r = {mantel_r:+.3f}",
            transform=ax.transAxes, fontsize=8, va="top", style="italic")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_varpart.pdf"))
    fig.savefig(os.path.join(HERE, "fig_varpart.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
