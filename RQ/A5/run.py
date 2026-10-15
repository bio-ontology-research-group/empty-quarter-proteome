#!/usr/bin/env python3
"""RQ-A5: Do dark families encode genuinely novel enzymatic chemistry?

Analyses the EC/pathway landscape of confidently decoded dark protein families
(C.V.T decode engine; FESNov-selected novel families, rodriguezdelrio2024fesnov).

Run from repo root:  python3 RQ/A5/run.py
Stdlib + matplotlib only.
"""
import csv, os, statistics
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.getcwd()
INP = os.path.join(REPO, "results/darkmatter_ecology/family_ec.tsv")
OUT = os.path.join(REPO, "RQ/A5")
os.makedirs(OUT, exist_ok=True)

EC_CLASS = {
    "1": "Oxidoreductases (1.x)", "2": "Transferases (2.x)",
    "3": "Hydrolases (3.x)", "4": "Lyases (4.x)",
    "5": "Isomerases (5.x)", "6": "Ligases (6.x)",
    "7": "Translocases (7.x)",
}


def parse_ecs(s):
    """best_ec is comma-separated; each token may carry '/'-separated alternatives."""
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        for alt in tok.split("/"):
            alt = alt.strip()
            if alt:
                out.append(alt)
    return out


def is_specific(ec):
    """Fully resolved 4-level EC (no missing '-' levels)."""
    parts = ec.split(".")
    return len(parts) == 4 and all(p and p != "-" for p in parts)


def parse_pathways(s):
    out = []
    for tok in s.split(","):
        tok = tok.strip().strip("|").strip()
        if tok:
            out.append(tok)
    return out


def main():
    rows = list(csv.DictReader(open(INP), delimiter="\t"))
    n = len(rows)
    conf = sorted(float(r["C_confidence"]) for r in rows)

    # --- 1. confidence distribution + cutoff ------------------------------
    q = lambda p: conf[int(len(conf) * p)]
    Q3 = q(0.75)
    CUTOFF = round(Q3, 3)  # top-quartile cutoff (defensible elbow)
    print(f"[conf] n={n} min={conf[0]:.4f} median={statistics.median(conf):.4f} "
          f"max={conf[-1]:.4f}")
    print(f"[conf] Q1={q(0.25):.4f} Q2={q(0.5):.4f} Q3={Q3:.4f} "
          f"P90={q(0.9):.4f} P95={q(0.95):.4f}")
    print(f"[conf] CUTOFF (top quartile) = {CUTOFF}")

    conf_rows = [r for r in rows if float(r["C_confidence"]) >= CUTOFF]
    # confident AND decoded to >=1 fully specified EC
    conf_specific = [r for r in conf_rows
                     if any(is_specific(e) for e in parse_ecs(r["best_ec"]))]
    print(f"[conf] families >= cutoff: {len(conf_rows)} "
          f"({100*len(conf_rows)/n:.1f}%)")
    print(f"[conf] of those, with a fully-specified EC: {len(conf_specific)} "
          f"({100*len(conf_specific)/n:.1f}% of all)")

    # --- 2. EC-class / pathway / core-secondary landscape -----------------
    ec_class = Counter()
    ec_full = Counter()
    pw = Counter()
    cs = Counter()
    for r in conf_specific:
        cs[r["core_secondary"]] += 1
        seen_class = set()
        for e in parse_ecs(r["best_ec"]):
            if is_specific(e):
                cls = e.split(".")[0]
                if cls not in seen_class:
                    ec_class[EC_CLASS.get(cls, f"EC {cls}.x")] += 1
                    seen_class.add(cls)
                ec_full[e] += 1
        for p in parse_pathways(r["pathways"]):
            pw[p] += 1

    print("\n[EC class] (families counted once per class)")
    for k, v in ec_class.most_common():
        print(f"  {k:28s} {v:6d}  ({100*v/len(conf_specific):.1f}%)")
    print(f"\n[core vs secondary] {dict(cs)}")
    print("\n[top pathways]")
    for k, v in pw.most_common(15):
        print(f"  {k:18s} {v}")

    # write EC-class table
    with open(os.path.join(OUT, "ec_class_landscape.tsv"), "w") as f:
        f.write("ec_class\tn_families\tpct_of_confident\n")
        for k, v in ec_class.most_common():
            f.write(f"{k}\t{v}\t{100*v/len(conf_specific):.2f}\n")

    # write pathway table
    with open(os.path.join(OUT, "pathway_counts.tsv"), "w") as f:
        f.write("pathway\tn_confident_families\n")
        for k, v in pw.most_common():
            f.write(f"{k}\t{v}\n")

    # write full EC table
    with open(os.path.join(OUT, "ec_full_counts.tsv"), "w") as f:
        f.write("ec\tn_confident_families\n")
        for k, v in ec_full.most_common():
            f.write(f"{k}\t{v}\n")

    # --- 3. candidate novel-chemistry shortlist ---------------------------
    # Compelling = confident + fully-specified EC + has a named reaction +
    # mapped to a pathway (pathway-completion candidate). Rank by confidence.
    # All these families are 'dark' (no canonical UniRef50 ortholog) by construction.
    cands = []
    for r in conf_specific:
        ecs = [e for e in parse_ecs(r["best_ec"]) if is_specific(e)]
        pws = parse_pathways(r["pathways"])
        if r["rxn_name"].strip() and pws:
            cands.append((float(r["C_confidence"]), r["family_id"],
                          ",".join(ecs), r["rxn_name"].strip(),
                          ";".join(pws), r["core_secondary"]))
    cands.sort(reverse=True)
    with open(os.path.join(OUT, "candidate_novel_chemistry.tsv"), "w") as f:
        f.write("family_id\tbest_ec\trxn_name\tpathways\tC_confidence\tcore_secondary\n")
        for c, fid, ec, rxn, p, scs in cands:
            f.write(f"{fid}\t{ec}\t{rxn}\t{p}\t{c:.4f}\t{scs}\n")
    print(f"\n[shortlist] {len(cands)} confident + named-rxn + pathway-mapped "
          f"candidate families -> candidate_novel_chemistry.tsv")
    print("[shortlist] top 15:")
    for c, fid, ec, rxn, p, scs in cands[:15]:
        print(f"  {c:.3f}  {ec:20s} {rxn[:45]:45s} {p[:30]}")

    # --- 4. figure --------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    items = ec_class.most_common()
    labels = [k.split(" (")[0] for k, _ in items]
    vals = [v for _, v in items]
    ax1.barh(range(len(labels)), vals, color="#3b6ea5")
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.invert_yaxis()
    ax1.set_xlabel("Confident dark families")
    ax1.set_title(f"EC-class landscape\n(confident decodes, C>={CUTOFF}, n={len(conf_specific)})",
                  fontsize=10)
    for i, v in enumerate(vals):
        ax1.text(v, i, f" {v}", va="center", fontsize=8)

    top_pw = pw.most_common(15)
    plabels = [k for k, _ in top_pw]
    pvals = [v for _, v in top_pw]
    ax2.barh(range(len(plabels)), pvals, color="#a5673b")
    ax2.set_yticks(range(len(plabels)))
    ax2.set_yticklabels(plabels, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("Confident dark families")
    ax2.set_title("Top decoded pathways", fontsize=10)
    fig.suptitle("RQ-A5  Enzymatic-chemistry landscape of decoded dark families",
                 fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_decode_landscape.pdf"))
    fig.savefig(os.path.join(OUT, "fig_decode_landscape.png"), dpi=150)
    print("\n[fig] fig_decode_landscape.pdf/png written")

    # stash key numbers for result.md
    with open(os.path.join(OUT, "_summary.txt"), "w") as f:
        f.write(f"n_total={n}\ncutoff={CUTOFF}\nn_confident={len(conf_rows)}\n"
                f"n_confident_specific={len(conf_specific)}\n"
                f"core={cs.get('core',0)}\nsecondary={cs.get('secondary',0)}\n"
                f"n_candidates={len(cands)}\n")
        for k, v in ec_class.most_common():
            f.write(f"ecclass\t{k}\t{v}\n")


if __name__ == "__main__":
    main()
