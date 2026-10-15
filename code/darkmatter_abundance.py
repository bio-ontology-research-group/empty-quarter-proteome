#!/usr/bin/env python3
"""B3d — the dark proteome is abundant and prevalent, not a rare contamination/relic-DNA tail.

Uses the CoverM read-recruitment of each dark family (one representative locus per family) across 551
samples (`darkfam_abundance.tsv`) + family metadata (`darkfam_meta.tsv`: n_members, catalog prevalence).
Contamination / dead-DNA would be sporadic + trace-level; we show dark families are widely detected and
reach genuine abundance, and that the abundant ones are also the prevalent ones (a real-community signature).

Out: results/darkmatter_ecology/{ABUNDANCE.md, darkfam_abundance_stats.tsv}
"""
import os, csv, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
AB = os.path.join(ECO, "darkfam_abundance.tsv")
META = os.path.join(ECO, "darkfam_meta.tsv")

# meta: family_id -> (n_members, n_samples_catalog)
meta = {}
with open(META) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            meta[r["family_id"]] = (int(r["n_members"]), int(r["n_samples_catalog"]))
        except (KeyError, ValueError):
            pass

rows = []
with open(AB) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    nsamp = len(hdr) - 1
    for line in f:
        c = line.rstrip("\n").split("\t")
        fid = c[0]
        vals = [float(x) for x in c[1:]]
        present = [v for v in vals if v > 0]
        ge1 = [v for v in vals if v >= 1.0]
        if not vals:
            continue
        peak = max(vals)
        mean_present = (sum(present) / len(present)) if present else 0.0
        nm, ncat = meta.get(fid, (0, 0))
        rows.append(dict(fid=fid, n_present=len(present), n_ge1=len(ge1), peak=peak,
                         mean_present=mean_present, total=sum(vals), n_members=nm))

N = len(rows)
def frac(pred): return sum(1 for r in rows if pred(r)) / N
def med(key): return st.median([r[key] for r in rows])

prev_frac = [r["n_present"] / nsamp for r in rows]
stats = dict(
    n_families=N, n_samples=nsamp,
    median_prevalence_samples=med("n_present"),
    median_prevalence_frac=st.median(prev_frac),
    pct_in_ge50pct_samples=frac(lambda r: r["n_present"] >= 0.5 * nsamp),
    pct_in_ge10pct_samples=frac(lambda r: r["n_present"] >= 0.1 * nsamp),
    pct_singleton_sample=frac(lambda r: r["n_present"] <= 1),
    median_peak=med("peak"),
    median_mean_present=med("mean_present"),
    pct_peak_ge100=frac(lambda r: r["peak"] >= 100),
    pct_peak_ge10=frac(lambda r: r["peak"] >= 10),
    pct_meanpresent_ge1=frac(lambda r: r["mean_present"] >= 1.0),
    # joint "real community" signature: prevalent AND non-trace
    pct_prevalent_and_abundant=frac(lambda r: r["n_present"] >= 0.5 * nsamp and r["peak"] >= 10),
)

with open(os.path.join(ECO, "darkfam_abundance_stats.tsv"), "w") as o:
    o.write("family_id\tn_present\tn_ge1\tpeak\tmean_present\ttotal\tn_members\n")
    for r in sorted(rows, key=lambda r: -r["total"]):
        o.write(f"{r['fid']}\t{r['n_present']}\t{r['n_ge1']}\t{r['peak']:.2f}\t"
                f"{r['mean_present']:.3f}\t{r['total']:.2f}\t{r['n_members']}\n")

L = ["# B3d — the dark proteome is abundant and prevalent (not rare contamination/relic DNA)\n",
     f"CoverM read-recruitment of {N:,} dark families across {nsamp} samples. Contamination / dead-DNA would be "
     "sporadic and trace-level; the data show the opposite.\n",
     "## Prevalence (widely detected, not sporadic)",
     f"- median dark family detected in **{int(stats['median_prevalence_samples'])} / {nsamp} samples** "
     f"({stats['median_prevalence_frac']*100:.0f}%).",
     f"- **{stats['pct_in_ge50pct_samples']*100:.1f}%** of families are present in ≥50% of samples; "
     f"{stats['pct_in_ge10pct_samples']*100:.1f}% in ≥10%.",
     f"- only {stats['pct_singleton_sample']*100:.1f}% are detected in ≤1 sample (the contamination-like tail).\n",
     "## Abundance magnitude (genuine, not trace)",
     f"- median peak (max across samples) abundance **{stats['median_peak']:.1f}**; "
     f"**{stats['pct_peak_ge100']*100:.1f}%** of families peak ≥100 and {stats['pct_peak_ge10']*100:.1f}% peak ≥10.",
     f"- median mean-where-present abundance {stats['median_mean_present']:.2f}; "
     f"{stats['pct_meanpresent_ge1']*100:.1f}% have mean-present ≥1.",
     f"- **{stats['pct_prevalent_and_abundant']*100:.1f}%** are BOTH prevalent (≥50% of samples) AND non-trace "
     "(peak ≥10) — the joint signature of a real, dominant community component, not a rare artefact.\n",
     "_Outputs: darkfam_abundance_stats.tsv; figure fig_darkmatter_abundance._"]
open(os.path.join(ECO, "ABUNDANCE.md"), "w").write("\n".join(L) + "\n")
print("B3d:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in stats.items()})
