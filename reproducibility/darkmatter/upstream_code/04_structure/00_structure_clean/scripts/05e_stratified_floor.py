#!/usr/bin/env python3
"""05e — confidence-standardised known-fold floor for §1.3.

Why. The headline 10.8x compares two pooled rates on cohorts that are not comparable on the one
covariate that drives the measurement. Gate-passing controls are far more confident than the
candidates (mean pLDDT 0.867 vs 0.775) and the control false-dual-novel rate falls steeply with
confidence, so the pooled floor is measured mostly where errors are rarest and the crude ratio is
inflated. The PI's own claim_boundary.txt asks for this matching before interpretation.

What it does. Direct standardisation: the control rate is computed within pLDDT strata (and within
pLDDT x length strata) and reweighted to the candidate covariate distribution, giving the floor the
same instrument would produce on cohorts shaped like ours. Reports crude and standardised floors,
the resulting multipliers, and bootstrap intervals. Nothing is reweighted on the candidate side;
the candidate rate is its own observed value throughout.

Both TM normalisations are reported, because the paper's claim rests on their contrast:
  qTM   dual no-hit = conf_qualified_dual_novel (as shipped)
  alnTM dual no-hit = best_pdb_alntm < TM and best_afdb_alntm < TM, on gate-passing models

Inputs (read-only): 05_validate/per_rep_structures.tsv, 05_validate/per_control_structures.tsv
Output: 05_validate/stratified_floor.tsv, and a summary to stdout + 05_validate/stratified_floor.txt

Usage: python3 05e_stratified_floor.py
"""
import os, csv, math, random

HERE = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(HERE, "..", "05_validate")
TM = 0.5
PLDDT_BINS = [(0.70, 0.75), (0.75, 0.80), (0.80, 0.85), (0.85, 0.90), (0.90, 1.01)]
LEN_BINS = [(100, 150), (150, 250), (250, 400), (400, 10 ** 9)]
NBOOT = 2000
SEED = 1


def load(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("confidence_qualified") != "1":
                continue
            rows.append({
                "plddt": float(r["mean_plddt"]),
                "length": int(r["length"]),
                "qtm_novel": int(r["conf_qualified_dual_novel"]),
                "alntm_novel": int(float(r["best_pdb_alntm"]) < TM and float(r["best_afdb_alntm"]) < TM),
            })
    return rows


def pbin(x):
    for i, (lo, hi) in enumerate(PLDDT_BINS):
        if lo <= x < hi:
            return i
    return len(PLDDT_BINS) - 1 if x >= PLDDT_BINS[-1][0] else 0


def lbin(x):
    for i, (lo, hi) in enumerate(LEN_BINS):
        if lo <= x < hi:
            return i
    return len(LEN_BINS) - 1


def rate(rows, key):
    return (sum(r[key] for r in rows) / len(rows)) if rows else float("nan")


def standardise(cand, ctrl, key, keyfn):
    """Control rate reweighted to the candidate distribution over keyfn strata.

    Strata with no control model contribute nothing and their candidate weight is removed from the
    denominator, so the result is a floor over the covered part of the candidate distribution; the
    uncovered weight is reported so it cannot hide.
    """
    cw, cr, covered = {}, {}, 0.0
    for r in cand:
        cw[keyfn(r)] = cw.get(keyfn(r), 0) + 1
    by = {}
    for r in ctrl:
        by.setdefault(keyfn(r), []).append(r)
    num = 0.0
    for k, w in cw.items():
        if k in by:
            num += (w / len(cand)) * rate(by[k], key)
            covered += w / len(cand)
    for k, v in by.items():
        cr[k] = (len(v), rate(v, key))
    return (num / covered if covered else float("nan")), covered, cr, cw


def boot(cand, ctrl, key, keyfn, n=NBOOT, seed=SEED):
    rnd = random.Random(seed)
    out = []
    for _ in range(n):
        c = [cand[rnd.randrange(len(cand))] for _ in range(len(cand))]
        k = [ctrl[rnd.randrange(len(ctrl))] for _ in range(len(ctrl))]
        f, cov, _, _ = standardise(c, k, key, keyfn)
        if f and f == f and f > 0:
            out.append(rate(c, key) / f)
    out.sort()
    if not out:
        return (float("nan"), float("nan"))
    return (out[int(0.025 * len(out))], out[int(0.975 * len(out))])


def main():
    cand = load(os.path.join(V, "per_rep_structures.tsv"))
    ctrl = load(os.path.join(V, "per_control_structures.tsv"))
    lines = []

    def say(s):
        print(s, flush=True)
        lines.append(s)

    say(f"gate-passing candidates {len(cand)}   controls {len(ctrl)}")
    say(f"mean pLDDT  candidates {sum(r['plddt'] for r in cand)/len(cand):.4f}   "
        f"controls {sum(r['plddt'] for r in ctrl)/len(ctrl):.4f}")
    say(f"mean length candidates {sum(r['length'] for r in cand)/len(cand):.1f}   "
        f"controls {sum(r['length'] for r in ctrl)/len(ctrl):.1f}")

    rows = []
    for key, label in (("qtm_novel", "query-normalised"), ("alntm_novel", "alignment-normalised")):
        c_rate, f_crude = rate(cand, key), rate(ctrl, key)
        say(f"\n--- {label} (TM {TM}) ---")
        say(f"candidate rate {c_rate:.4f} ({sum(r[key] for r in cand)}/{len(cand)})   "
            f"crude floor {f_crude:.4f} ({sum(r[key] for r in ctrl)}/{len(ctrl)})   "
            f"crude ratio {c_rate/f_crude:.2f}x")
        rows.append(dict(normalisation=label, adjustment="crude", floor=round(f_crude, 6),
                         multiplier=round(c_rate / f_crude, 3), candidate_rate=round(c_rate, 6),
                         covered_weight=1.0, lo="", hi=""))
        for adj, fn in (("pLDDT", lambda r: pbin(r["plddt"])),
                        ("pLDDT x length", lambda r: (pbin(r["plddt"]), lbin(r["length"])))):
            f, cov, cr, cw = standardise(cand, ctrl, key, fn)
            lo, hi = boot(cand, ctrl, key, fn)
            say(f"standardised on {adj:15s} floor {f:.4f}  ratio {c_rate/f:.2f}x "
                f"(95% CI {lo:.2f}-{hi:.2f})  candidate weight covered {cov:.3f}")
            rows.append(dict(normalisation=label, adjustment=adj, floor=round(f, 6),
                             multiplier=round(c_rate / f, 3), candidate_rate=round(c_rate, 6),
                             covered_weight=round(cov, 4), lo=round(lo, 3), hi=round(hi, 3)))
        if key == "qtm_novel":
            # recompute on pLDDT alone for the display, so the strata keys are the bin indices
            _, _, cr_p, cw_p = standardise(cand, ctrl, key, lambda r: pbin(r["plddt"]))
            say("control false-call rate by pLDDT stratum, against the candidate share:")
            for i, (lo_, hi_) in enumerate(PLDDT_BINS):
                n, rr = cr_p.get(i, (0, float("nan")))
                say(f"  {lo_:.2f}-{hi_:.2f}  control n={n:5d} rate={rr:.4f}   "
                    f"candidate n={cw_p.get(i, 0):5d} ({100*cw_p.get(i, 0)/len(cand):.1f} %)")

    with open(os.path.join(V, "stratified_floor.tsv"), "w") as o:
        head = ["normalisation", "adjustment", "candidate_rate", "floor", "multiplier",
                "covered_weight", "lo", "hi"]
        o.write("\t".join(head) + "\n")
        for r in rows:
            o.write("\t".join(str(r[k]) for k in head) + "\n")
    with open(os.path.join(V, "stratified_floor.txt"), "w") as o:
        o.write("\n".join(lines) + "\n")
    say("\nwrote 05_validate/stratified_floor.{tsv,txt}")


if __name__ == "__main__":
    main()
