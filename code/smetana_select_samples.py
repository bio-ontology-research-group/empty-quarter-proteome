#!/usr/bin/env python3
"""Pick representative samples for the SMETANA cross-feeding pilot.

We want a small, balanced panel that lets us ask: (a) is the community metabolically
interdependent (cross-feeding) vs. independent generalists, and (b) are the atmospheric
"life-on-air" MAGs the metabolic donors? To answer both we need samples that:
  - span the three compartments (Surface UV-exposed / Deep shielded / PlantRoot rhizosphere),
  - span the aridity gradient (so we can see whether interdependence tracks aridity),
  - contrast guild-RICH vs guild-POOR communities (the guild-as-feeder test),
  - have enough HQ MAGs to form a real community but few enough to model fast (cap on build cost).

Community = HQ MAGs only (CheckM2 completeness>=90, contamination<5) — gating on completeness
is the single most important control against gap-fill-inflated false cross-feeding.

In:  workbench/results/soil_gtdbtk/checkm2_all_samples.tsv  (per-MAG quality)
     results/darkmatter_ecology/env_matrix.tsv               (compartment, aridity, geochem)
     results/atmospheric/guild_core_mags.txt                 (core guild MAG ids)
     results/atmospheric/abundance/guild_abundance_by_sample.tsv (guild read-fraction)
Out: results/smetana/pilot_samples.tsv         (the chosen panel + why)
     results/smetana/pilot_mags.tsv            (MAG -> sample -> is_guild, the GEM build list)
"""
import os, csv
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
WB = os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", "checkm2_all_samples.tsv")
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology", "env_matrix.tsv")
ATMO = os.path.join(HERE, "..", "results", "atmospheric")
OUT = os.path.join(HERE, "..", "results", "smetana")
os.makedirs(OUT, exist_ok=True)

COMP_MIN, CONTAM_MAX = 90.0, 5.0
SIZE_LO, SIZE_HI = 12, 45          # community size window: big enough for signal, small enough to build
N_PER_CELL = 2                      # samples per (compartment x guild-tertile) cell where possible

# ---- HQ MAGs per sample ----
hq = defaultdict(list)             # sample -> [SemiBin_N,...]
for r in csv.DictReader(open(WB), delimiter="\t"):
    try:
        if float(r["Completeness"]) >= COMP_MIN and float(r["Contamination"]) < CONTAM_MAX:
            hq[r["sample"]].append(r["Name"])
    except (ValueError, KeyError):
        continue

# ---- env (compartment, aridity, geochem availability) ----
env = {}
for r in csv.DictReader(open(ECO), delimiter="\t"):
    env[r["sample"]] = r

# ---- guild MAGs per sample + guild read-fraction ----
guild_mags = defaultdict(set)      # sample -> {SemiBin_N guild}
for ln in open(os.path.join(ATMO, "guild_core_mags.txt")):
    mid = ln.strip()
    if "__" in mid:
        s, b = mid.split("__", 1)
        guild_mags[s].add(b)

guild_frac = {}                    # sample -> guild read fraction (full/genus/any)
for r in csv.DictReader(open(os.path.join(ATMO, "abundance", "guild_abundance_by_sample.tsv")), delimiter="\t"):
    if r["table"] == "full" and r["rank"] == "genus" and r["metric"] == "any":
        s = r["sample"]
        s = s.split("_", 1)[1] if s[:1] == "e" and "_" in s else s
        try: guild_frac[s] = float(r["guild_fraction"])
        except ValueError: pass

# ---- candidate samples: in size window, known compartment, known aridity ----
cand = []
for s, mags in hq.items():
    e = env.get(s)
    if not e or e.get("compartment") not in ("Surface", "Deep", "PlantRoot"):
        continue
    if not (SIZE_LO <= len(mags) <= SIZE_HI):
        continue
    try: arid = float(e["aridity_index"])
    except (ValueError, KeyError): continue
    gm = guild_mags.get(s, set()) & set(mags)        # guild MAGs that are also HQ
    cand.append(dict(sample=s, comp=e["compartment"], n_hq=len(mags), arid=arid,
                     n_guild_hq=len(gm), guild_share=len(gm) / len(mags),
                     guild_frac=guild_frac.get(s, float("nan")), has_xrf=e.get("has_xrf", "0"),
                     V=e.get("V", ""), Na=e.get("Na", ""), mags=mags, guild=gm))

# guild tertiles by guild_share (within candidates)
shares = sorted(c["guild_share"] for c in cand)
def tertile(v):
    n = len(shares)
    if v <= shares[n // 3]: return "low"
    if v <= shares[2 * n // 3]: return "high" if v >= shares[2 * n // 3] else "mid"
    return "high"
for c in cand:
    c["guild_tier"] = tertile(c["guild_share"])

# ---- balanced pick: for each (compartment x guild_tier), take up to N_PER_CELL samples,
#      preferring XRF-available + aridity spread (alternate extremes), then more guild MAGs ----
cells = defaultdict(list)
for c in cand:
    cells[(c["comp"], c["guild_tier"])].append(c)

chosen = []
for (comp, tier), lst in sorted(cells.items()):
    # rank: XRF present first, then by |aridity - median| to spread, then community size sweet-spot
    amed = sorted(x["arid"] for x in lst)[len(lst) // 2]
    lst.sort(key=lambda x: (x["has_xrf"] != "1", -abs(x["arid"] - amed),
                            -x["n_guild_hq"] if tier != "low" else x["n_guild_hq"]))
    chosen.extend(lst[:N_PER_CELL])

# de-dup + cap total ~12, keep compartment balance
seen, final = set(), []
for c in sorted(chosen, key=lambda x: (x["comp"], x["guild_tier"], -x["n_hq"])):
    if c["sample"] in seen: continue
    seen.add(c["sample"]); final.append(c)

# ---- write panel ----
cols = ["sample", "comp", "n_hq", "n_guild_hq", "guild_share", "guild_tier", "guild_frac",
        "arid", "has_xrf", "V", "Na"]
with open(os.path.join(OUT, "pilot_samples.tsv"), "w") as o:
    o.write("\t".join(cols) + "\n")
    for c in final:
        o.write("\t".join(str(round(c[k], 4) if isinstance(c[k], float) else c[k]) for k in cols) + "\n")

with open(os.path.join(OUT, "pilot_mags.tsv"), "w") as o:
    o.write("sample\tmag\tmag_id\tis_guild\tcompartment\n")
    for c in final:
        for b in sorted(c["mags"]):
            o.write(f"{c['sample']}\t{b}\t{c['sample']}__{b}\t{int(b in c['guild'])}\t{c['comp']}\n")

n_mags = sum(len(c["mags"]) for c in final)
print(f"candidates in window: {len(cand)}; chosen: {len(final)} samples, {n_mags} HQ MAGs to model")
print(f"{'sample':10s} {'comp':10s} {'nHQ':>4s} {'guildHQ':>7s} {'share':>6s} {'tier':>4s} {'arid':>5s} xrf")
for c in final:
    print(f"{c['sample']:10s} {c['comp']:10s} {c['n_hq']:4d} {c['n_guild_hq']:7d} "
          f"{c['guild_share']:.2f}   {c['guild_tier']:>4s} {c['arid']:.2f}  {c['has_xrf']}")
by_comp = defaultdict(int)
for c in final: by_comp[c["comp"]] += 1
print("by compartment:", dict(by_comp))
