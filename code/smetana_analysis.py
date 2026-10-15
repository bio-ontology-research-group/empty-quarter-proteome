#!/usr/bin/env python3
"""Cross-feeding analysis of the SMETANA pilot (16 per-sample MAG communities).

Three questions:
  (1) Is the community metabolically interdependent? -> community MIP (cooperation potential) and
      MRO (resource overlap/competition) from SMETANA global mode.
  (2) Are the atmospheric "life-on-air" guild MAGs the metabolic DONORS (feeders) of the community?
      -> per-MAG donor vs receiver strength from SMETANA detailed mode; guild vs non-guild test.
  (3) Does interdependence track the aridity/compartment gradient? -> MIP/MRO vs aridity, compartment.

CONTROL (the key reviewer concern): MAG incompleteness -> gap-filled false auxotrophies -> inflated
apparent cross-feeding. We regress per-MAG receiver strength on CheckM2 completeness; a strong
positive slope would mean the cross-feeding signal is a completeness artifact, not biology.

In:  results/smetana/smetana_out/global/*_global.tsv      (community  medium  size  mip  mro)
     results/smetana/smetana_out/detailed/*_detailed.tsv  (community medium receiver donor compound scs mus mps smetana)
     results/smetana/pilot_mags.tsv  pilot_samples.tsv
     workbench/results/soil_gtdbtk/checkm2_all_samples.tsv  (completeness)
Out: results/smetana/{community_summary.tsv, mag_donor_receiver.tsv, crossfeed_top.tsv, SMETANA.md}
"""
import os, csv, glob
from collections import defaultdict
from statistics import mean, median

HERE = os.path.dirname(os.path.abspath(__file__))
SM = os.path.join(HERE, "..", "results", "smetana")
GLOB = os.path.join(SM, "smetana_out", "global")
DET = os.path.join(SM, "smetana_out", "detailed")
CK = os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", "checkm2_all_samples.tsv")

def sp(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

def spearman(xs, ys):
    pts = [(x, y) for x, y in zip(xs, ys) if x == x and y == y]
    n = len(pts)
    if n < 4: return float("nan"), n
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1): r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank([p[0] for p in pts]), rank([p[1] for p in pts])
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return (num / den if den else float("nan")), n

def mwu(a, b):
    """Mann-Whitney U, return (U, p_two_sided approx, p_one_sided a>b)."""
    a = [x for x in a if x == x]; b = [x for x in b if x == x]
    if len(a) < 3 or len(b) < 3: return float("nan"), float("nan"), float("nan")
    allv = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = {}; i = 0; vals = [v for v, _ in allv]
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]: j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1): ranks[k] = avg
        i = j + 1
    Ra = sum(ranks[k] for k, (_, g) in enumerate(allv) if g == 0)
    na, nb = len(a), len(b)
    Ua = Ra - na * (na + 1) / 2.0
    mu = na * nb / 2.0; sd = (na * nb * (na + nb + 1) / 12.0) ** 0.5
    if sd == 0: return Ua, float("nan"), float("nan")
    z = (Ua - mu) / sd
    # normal approx
    import math
    p_two = math.erfc(abs(z) / 2 ** 0.5)
    p_one = 0.5 * math.erfc(z / 2 ** 0.5)   # a > b
    return Ua, p_two, p_one

# ---- metadata ----
mag_meta = {}   # mag_id -> dict(sample, is_guild, compartment)
for r in csv.DictReader(open(os.path.join(SM, "pilot_mags.tsv")), delimiter="\t"):
    mag_meta[r["mag_id"]] = dict(sample=r["sample"], is_guild=int(r["is_guild"]), comp=r["compartment"])
samp_env = {r["sample"]: r for r in csv.DictReader(open(os.path.join(SM, "pilot_samples.tsv")), delimiter="\t")}
compl = {}      # mag_id -> completeness
for r in csv.DictReader(open(CK), delimiter="\t"):
    compl[f"{r['sample']}__{r['Name']}"] = sp(r["Completeness"])

# ---- global MIP/MRO ----
comm = {}       # sample -> dict(mip, mro, size)
for f in glob.glob(os.path.join(GLOB, "*_global.tsv")) + glob.glob(os.path.join(GLOB, "*.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        s = r.get("community", "").strip()
        if not s: continue
        comm[s] = dict(mip=sp(r.get("mip")), mro=sp(r.get("mro")), size=sp(r.get("size")))

# ---- detailed cross-feeding edges ----
# per MAG: donor strength (sum smetana as donor), receiver strength (sum as receiver)
donor = defaultdict(float); receiver = defaultdict(float); n_as_donor = defaultdict(int)
edges = []      # (community, donor, receiver, compound, smetana)
comm_xf = defaultdict(float)   # sample -> total smetana
for f in glob.glob(os.path.join(DET, "*_detailed.tsv")) + glob.glob(os.path.join(DET, "*.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        sc = sp(r.get("smetana"))
        if not (sc == sc): continue
        d = r.get("donor", "").strip(); rc = r.get("receiver", "").strip()
        cpd = r.get("compound", "").strip(); cm = r.get("community", "").strip()
        if not d or not rc: continue
        donor[d] += sc; receiver[rc] += sc; n_as_donor[d] += (sc > 0.1)
        comm_xf[cm] += sc
        if sc >= 0.1: edges.append((cm, d, rc, cpd, sc))

# ---- per-MAG donor/receiver table ----
rows = []
for mag, m in mag_meta.items():
    ds, rs = donor.get(mag, 0.0), receiver.get(mag, 0.0)
    tot = ds + rs
    rows.append(dict(mag=mag, sample=m["sample"], comp=m["comp"], is_guild=m["is_guild"],
                     completeness=compl.get(mag, float("nan")), donor=ds, receiver=rs,
                     donor_frac=(ds / tot if tot > 0 else float("nan"))))
rows.sort(key=lambda r: -r["donor"])
with open(os.path.join(SM, "mag_donor_receiver.tsv"), "w") as o:
    o.write("mag\tsample\tcompartment\tis_guild\tcompleteness\tdonor_strength\treceiver_strength\tdonor_frac\n")
    for r in rows:
        o.write(f"{r['mag']}\t{r['sample']}\t{r['comp']}\t{r['is_guild']}\t"
                f"{r['completeness']:.1f}\t{r['donor']:.3f}\t{r['receiver']:.3f}\t"
                f"{r['donor_frac'] if r['donor_frac']!=r['donor_frac'] else round(r['donor_frac'],3)}\n")

# ---- community summary ----
csum = []
for s, e in samp_env.items():
    c = comm.get(s, {});
    mags_here = [r for r in rows if r["sample"] == s]
    csum.append(dict(sample=s, comp=e["compartment"], arid=sp(e["aridity_index"]),
                     guild_share=sp(e["guild_share"]), n_mags=len(mags_here),
                     n_guild=sum(r["is_guild"] for r in mags_here),
                     mip=c.get("mip", float("nan")), mro=c.get("mro", float("nan")),
                     total_xf=comm_xf.get(s, 0.0)))
csum.sort(key=lambda r: r["sample"])
with open(os.path.join(SM, "community_summary.tsv"), "w") as o:
    o.write("sample\tcompartment\taridity\tguild_share\tn_mags\tn_guild\tMIP\tMRO\ttotal_crossfeed\n")
    for r in csum:
        o.write(f"{r['sample']}\t{r['comp']}\t{r['arid']:.2f}\t{round(r['guild_share'],3)}\t{r['n_mags']}\t"
                f"{r['n_guild']}\t{r['mip']}\t{r['mro']}\t{round(r['total_xf'],2)}\n")

# top cross-feeding edges
edges.sort(key=lambda e: -e[4])
with open(os.path.join(SM, "crossfeed_top.tsv"), "w") as o:
    o.write("community\tdonor\treceiver\tcompound\tsmetana\tdonor_is_guild\n")
    for cm, d, rc, cpd, sc in edges[:200]:
        o.write(f"{cm}\t{d}\t{rc}\t{cpd}\t{round(sc,3)}\t{mag_meta.get(d,{}).get('is_guild','')}\n")

# ---- tests ----
g_donor = [r["donor"] for r in rows if r["is_guild"]]
ng_donor = [r["donor"] for r in rows if not r["is_guild"]]
g_df = [r["donor_frac"] for r in rows if r["is_guild"]]
ng_df = [r["donor_frac"] for r in rows if not r["is_guild"]]
U1, p2_1, p1_1 = mwu(g_donor, ng_donor)         # guild donor strength > non-guild?
U2, p2_2, p1_2 = mwu(g_df, ng_df)               # guild donor fraction > non-guild?
# completeness control
rho_rc, n_rc = spearman([r["completeness"] for r in rows], [r["receiver"] for r in rows])
rho_dc, n_dc = spearman([r["completeness"] for r in rows], [r["donor"] for r in rows])
# gradient
rho_mip_ar, _ = spearman([r["arid"] for r in csum], [r["mip"] for r in csum])
rho_mro_ar, _ = spearman([r["arid"] for r in csum], [r["mro"] for r in csum])
rho_mip_gs, _ = spearman([r["guild_share"] for r in csum], [r["mip"] for r in csum])

valid_mip = [r["mip"] for r in csum if r["mip"] == r["mip"]]
valid_mro = [r["mro"] for r in csum if r["mro"] == r["mro"]]

ml = ["# SMETANA cross-feeding pilot — 16 per-sample Empty Quarter MAG communities\n",
      f"GEMs: gapseq (auto-predicted media) on HQ single-sample MAGs. Communities with SMETANA "
      f"global output: {len(comm)}/16; with detailed cross-feeding: "
      f"{len(set(e[0] for e in edges))}/16. Cross-feeding edges (smetana>=0.1): {len(edges)}.\n",
      "## (1) Are communities interdependent? (global MIP/MRO)\n",
      f"MIP (cooperation potential, # metabolites producible only by the community): "
      f"median {median(valid_mip):.1f} (range {min(valid_mip):.0f}-{max(valid_mip):.0f}) "
      f"across {len(valid_mip)} communities." if valid_mip else "No global MIP output yet.",
      f"MRO (resource overlap / competition): median {median(valid_mro):.2f} "
      f"(range {min(valid_mro):.2f}-{max(valid_mro):.2f})." if valid_mro else "",
      "\n## (2) Are the atmospheric guild MAGs the metabolic donors?\n",
      f"Guild MAGs (n={len(g_donor)}) vs non-guild (n={len(ng_donor)}):",
      f"- donor strength: guild median {median(g_donor):.3f} vs non-guild {median(ng_donor):.3f}; "
      f"MWU one-sided (guild>non-guild) p={p1_1:.3f}." if g_donor and ng_donor else "- insufficient data.",
      f"- donor FRACTION (donor/(donor+receiver)): guild {median([x for x in g_df if x==x] or [float('nan')]):.3f} "
      f"vs non-guild {median([x for x in ng_df if x==x] or [float('nan')]):.3f}; MWU one-sided p={p1_2:.3f}.",
      "\n## (3) Gradient — does interdependence track aridity / guild content?\n",
      f"- MIP vs aridity: Spearman rho={rho_mip_ar:.2f}",
      f"- MRO vs aridity: Spearman rho={rho_mro_ar:.2f}",
      f"- MIP vs guild_share: Spearman rho={rho_mip_gs:.2f}",
      "\n## CONTROL — is cross-feeding just a MAG-incompleteness artifact?\n",
      f"Per-MAG receiver strength vs CheckM2 completeness: Spearman rho={rho_rc:.2f} (n={n_rc}). "
      f"Donor strength vs completeness: rho={rho_dc:.2f}. A strongly NEGATIVE receiver~completeness "
      f"rho would mean incomplete MAGs look artificially auxotrophic (uptake-biased) — i.e. the "
      f"cross-feeding signal is partly a gap-fill artifact. |rho|<~0.3 supports a biological signal. "
      f"All MAGs here are HQ (>=90% complete), which already minimises this.",
      "\n## Caveats",
      "- GEMs see only the DECODABLE metabolism; the dark proteome (the paper's core) is invisible "
      "to SMETANA, so this describes the known metabolic fraction only.",
      "- Medium = SMETANA 'complete' environment (all exchanges available) => an UPPER BOUND on "
      "cross-feeding capacity; a defined desert-minimal medium would tighten (reduce) it.",
      "- gapseq auto-predicts each MAG's gap-fill medium independently; community context is imposed "
      "only at the SMETANA step, not during reconstruction.",
      "- Pilot = 16 of 547 samples; significance here is a feasibility/effect-size read, not the final test."]
open(os.path.join(SM, "SMETANA.md"), "w").write("\n".join(x for x in ml if x) + "\n")

print(f"communities global={len(comm)} detailed={len(set(e[0] for e in edges))}; edges={len(edges)}")
print(f"guild donor median={median(g_donor):.3f} vs non-guild={median(ng_donor):.3f} (MWU p1={p1_1:.3f})" if g_donor and ng_donor else "no donor data")
print(f"completeness control: receiver~compl rho={rho_rc:.2f}, donor~compl rho={rho_dc:.2f}")
print(f"MIP~aridity rho={rho_mip_ar:.2f}; MIP~guild_share rho={rho_mip_gs:.2f}")
if valid_mip: print(f"MIP median={median(valid_mip):.1f}; MRO median={median(valid_mro):.2f}" if valid_mro else "")
