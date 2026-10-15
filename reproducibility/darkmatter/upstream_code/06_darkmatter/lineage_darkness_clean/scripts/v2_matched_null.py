#!/usr/bin/env python3
"""v2_matched_null — is Thermoplasmatota darkness a rarity or quality artefact?

Maeke 2025 (Microbiome, 10.1186/s40168-025-02140-8) reports the same pattern in a sister
Thermoplasmatota class and attributes it to MAG rarity: seldom-recovered lineages are recovered
badly, and badly recovered genomes look unknown. The contrast computed here previously (0.310 vs
0.125) does not answer that, because its comparators are all bacterial and better assembled, so it
confounds phylum with domain and with genome quality.

Design: 1:1 nearest-neighbour matching without replacement on standardised covariates, against two
pools reported separately. The archaeal pool is the control that matters, since it removes the
domain term. Covariates: completeness, contamination, log10 corrected size, log10 N50, log10 mean
abundance, prevalence. Rarity enters through the last two, quality through the rest.

Inputs (read-only): results/v2_reps.tsv and the workbench abundance matrix.
Outputs: results/v2_matched_null.tsv, results/v2_matched_null_summary.txt
"""
import os, csv, math, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
ABUND = ("/ibex/scratch/projects/c2014/rob/metagenomics/10_thermal/abundance/"
         "rep_relative_abundance.tsv")
TARGET = "Thermoplasmatota"
COVS = ["completeness", "contamination", "log_size", "log_n50", "log_abund", "prevalence"]
CALIPERS = [0.25, 0.5, 1.0, 2.0]


def load_reps():
    out = {}
    with open(os.path.join(RES, "v2_reps.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            comp, ln = float(r["completeness"]), float(r["length"])
            if comp <= 0:
                continue
            out[r["genome"]] = {
                "genome": r["genome"],
                "phylum": (r["phylum"] or "").replace("p__", ""),
                "domain": (r["domain"] or "").replace("d__", ""),
                "completeness": comp,
                "contamination": float(r["contamination"]),
                "log_size": math.log10(ln / (comp / 100.0)),
                "log_n50": math.log10(max(float(r["N50"]), 1)),
                "dark": float(r["dark_fraction"]),
            }
    return out


def add_abundance(reps):
    seen = 0
    with open(ABUND) as fh:
        ncol = len(fh.readline().rstrip("\n").split("\t")) - 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = f[0][:-3] if f[0].endswith(".fa") else f[0]
            if g not in reps:
                continue
            seen += 1
            # the matrix carries literal "nan" cells (failed sample columns); treat as absent
            v = []
            for x in f[1:]:
                try:
                    y = float(x)
                except ValueError:
                    y = 0.0
                v.append(y if y == y else 0.0)
            reps[g]["prevalence"] = sum(1 for x in v if x > 0)
            reps[g]["log_abund"] = math.log10(sum(v) / ncol + 1e-6)
    for g in [g for g, r in reps.items() if "prevalence" not in r]:
        del reps[g]
    return seen


def scales(rows):
    return {c: (st.pstdev([r[c] for r in rows]) or 1.0) for c in COVS}


def dist(a, b, s):
    return math.sqrt(sum(((a[c] - b[c]) / s[c]) ** 2 for c in COVS) / len(COVS))


def match(tgt, pool, s, caliper):
    cand = []
    for t in tgt:
        for p in pool:
            d = dist(t, p, s)
            if d <= caliper:
                cand.append((d, t["genome"], p["genome"]))
    cand.sort()
    byg = {r["genome"]: r for r in tgt + pool}
    used, done, pairs = set(), set(), []
    for d, tg, pg in cand:
        if tg in done or pg in used:
            continue
        done.add(tg); used.add(pg)
        pairs.append((byg[tg], byg[pg], d))
    return pairs


def smd(a, b, c, s):
    return (st.mean(x[c] for x in a) - st.mean(x[c] for x in b)) / s[c]


def wilcoxon(d):
    d = [x for x in d if x != 0]
    n = len(d)
    if n < 10:
        return float("nan")
    order = sorted(range(n), key=lambda i: abs(d[i]))
    ranks = [0.0] * n
    for rank, i in enumerate(order, 1):
        ranks[i] = rank
    wp = sum(ranks[i] for i in range(n) if d[i] > 0)
    mu, sd = n * (n + 1) / 4, math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (wp - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def main():
    reps = load_reps()
    add_abundance(reps)
    rows = list(reps.values())
    tgt = [r for r in rows if r["phylum"] == TARGET]
    pools = {"archaea": [r for r in rows if r["domain"] == "Archaea" and r["phylum"] != TARGET],
             "bacteria": [r for r in rows if r["domain"] == "Bacteria"]}
    out, pairs_all = [], []

    def say(s):
        print(s, flush=True)
        out.append(s)

    say(f"{TARGET} {len(tgt)} reps; pools " + ", ".join(f"{k} {len(v)}" for k, v in pools.items()))
    say(f"unmatched median dark: {TARGET} {st.median(r['dark'] for r in tgt):.4f}  "
        + "  ".join(f"{k} {st.median(r['dark'] for r in v):.4f}" for k, v in pools.items()))

    for name, pool in pools.items():
        s = scales(tgt + pool)
        say(f"\n=== matched against {name} ===")
        say("  unmatched imbalance (SMD): " +
            "  ".join(f"{c} {smd(tgt, pool, c, s):+.2f}" for c in COVS))
        for cal in CALIPERS:
            pairs = match(tgt, pool, s, cal)
            if len(pairs) < 30:
                say(f"  caliper {cal}: only {len(pairs)} pairs, skipped")
                continue
            a, b = [p[0] for p in pairs], [p[1] for p in pairs]
            diff = [x["dark"] - y["dark"] for x, y in zip(a, b)]
            worst = max(abs(smd(a, b, c, s)) for c in COVS)
            say(f"  caliper {cal}: {len(pairs)} pairs of {len(tgt)}   "
                f"median {st.median(x['dark'] for x in a):.4f} vs "
                f"{st.median(x['dark'] for x in b):.4f}   "
                f"ratio {st.median(x['dark'] for x in a)/st.median(x['dark'] for x in b):.2f}x   "
                f"paired median {st.median(diff):+.4f}  darker in "
                f"{sum(1 for x in diff if x > 0)}/{len(diff)}  p={wilcoxon(diff):.1e}  "
                f"worst residual SMD {worst:.2f}")
            if cal == 1.0:
                say("    balance after matching: " +
                    "  ".join(f"{c} {smd(a, b, c, s):+.2f}" for c in COVS))
                for t, c_, d in pairs:
                    pairs_all.append(dict(pool=name, target=t["genome"], control=c_["genome"],
                                          control_phylum=c_["phylum"], distance=round(d, 4),
                                          target_dark=t["dark"], control_dark=c_["dark"],
                                          target_completeness=t["completeness"],
                                          control_completeness=c_["completeness"],
                                          target_prevalence=t["prevalence"],
                                          control_prevalence=c_["prevalence"]))

    if pairs_all:
        with open(os.path.join(RES, "v2_matched_null.tsv"), "w") as o:
            head = list(pairs_all[0].keys())
            o.write("\t".join(head) + "\n")
            for r in pairs_all:
                o.write("\t".join(str(r[k]) for k in head) + "\n")
    with open(os.path.join(RES, "v2_matched_null_summary.txt"), "w") as o:
        o.write("\n".join(out) + "\n")
    say("\nwrote results/v2_matched_null.tsv and _summary.txt")


if __name__ == "__main__":
    main()
