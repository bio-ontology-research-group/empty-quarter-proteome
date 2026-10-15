#!/usr/bin/env python3
"""Bin the RAK-PR cultured genomes against the 5,438 soil MAG species representatives.

Bins (as in W/07_aquatic_paper/07_endemicity):
  shared  : ANI >= 95 and aligned fraction >= 30
  related : 80 <= ANI < 95 (or ANI >= 95 with AF < 30)
  none    : ANI < 80, or no hit at all at the -s 75 screen

Forward direction (per cultured genome) uses the query aligned fraction; reverse direction
(per MAG representative) uses the reference aligned fraction, i.e. each genome is judged on how
much of itself aligns.

Outputs: cultured_best_hit.tsv, mag_best_cultured_hit.tsv, shared_species.tsv, summary.tsv
"""
import csv, os, statistics

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
ABUND = os.path.join(W, "10_thermal/abundance/rep_relative_abundance.tsv")
SKANI = os.path.join(P, "results/skani_search.tsv")

ANI_SHARED, ANI_RELATED, AF_MIN = 95.0, 80.0, 30.0


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def binof(ani, af):
    if ani >= ANI_SHARED and af >= AF_MIN:
        return "shared"
    if ani >= ANI_RELATED:
        return "related"
    return "none"


cult = read_tsv(os.path.join(P, "results/cultured_meta.tsv"))
mag = read_tsv(os.path.join(P, "results/mag_meta.tsv"))
by_path = {r["path"]: r for r in cult}
mag_tax = {r["mag_rep"]: r["mag_gtdb"] for r in mag}
print(f"cultured genomes: {len(cult)}   MAG reps: {len(mag)}")

# ---- parse the pair table ---------------------------------------------------------------------
pairs = []
for r in read_tsv(SKANI):
    q = by_path[r["Query_file"]]["genome_id"]
    ref = os.path.basename(r["Ref_file"])[:-3]      # strip .fa
    assert ref in mag_tax, ref
    pairs.append({"cultured": q, "mag_rep": ref, "ani": float(r["ANI"]),
                  "af_ref": float(r["Align_fraction_ref"]),
                  "af_query": float(r["Align_fraction_query"])})
print(f"pairs at >= 75 % ANI: {len(pairs)}")

# ---- forward: best hit per cultured genome -----------------------------------------------------
best_c = {}
for p in pairs:
    cur = best_c.get(p["cultured"])
    if cur is None or p["ani"] > cur["ani"] or (p["ani"] == cur["ani"] and p["af_query"] > cur["af_query"]):
        best_c[p["cultured"]] = p

frows = []
for c in cult:
    p = best_c.get(c["genome_id"])
    b = binof(p["ani"], p["af_query"]) if p else "none"
    frows.append({
        "cultured_genome": c["genome_id"], "batch": c["batch"],
        "source_sample": c["source_sample"], "site": c["site"],
        "quality_tier": c["quality_tier"], "cultured_gtdb": c["cultured_gtdb"],
        "best_mag_rep": p["mag_rep"] if p else "",
        "ani": f"{p['ani']:.2f}" if p else "",
        "af_query": f"{p['af_query']:.2f}" if p else "",
        "af_ref": f"{p['af_ref']:.2f}" if p else "",
        "mag_gtdb": mag_tax[p["mag_rep"]] if p else "",
        "bin": b,
    })
frows.sort(key=lambda r: (-float(r["ani"] or 0), r["cultured_genome"]))
with open(os.path.join(P, "results/cultured_best_hit.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(frows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(frows)

# ---- reverse: best cultured hit per MAG representative -----------------------------------------
best_m = {}
for p in pairs:
    cur = best_m.get(p["mag_rep"])
    if cur is None or p["ani"] > cur["ani"] or (p["ani"] == cur["ani"] and p["af_ref"] > cur["af_ref"]):
        best_m[p["mag_rep"]] = p

cult_by_id = {r["genome_id"]: r for r in cult}
rrows = []
for m in mag:
    p = best_m.get(m["mag_rep"])
    if p is None:
        continue                       # only the MAG reps that have any cultured relative
    c = cult_by_id[p["cultured"]]
    rrows.append({
        "mag_rep": m["mag_rep"], "mag_gtdb": m["mag_gtdb"], "domain": m["domain"],
        "best_cultured_genome": p["cultured"], "batch": c["batch"], "site": c["site"],
        "cultured_gtdb": c["cultured_gtdb"],
        "ani": f"{p['ani']:.2f}", "af_ref": f"{p['af_ref']:.2f}",
        "af_query": f"{p['af_query']:.2f}",
        "bin": binof(p["ani"], p["af_ref"]),
    })
rrows.sort(key=lambda r: -float(r["ani"]))
with open(os.path.join(P, "results/mag_best_cultured_hit.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rrows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rrows)
n_mag_shared = sum(1 for r in rrows if r["bin"] == "shared")
print(f"MAG reps with any cultured relative >= 75 %: {len(rrows)}   of them shared species: {n_mag_shared}")

# ---- shared species pairs, with soil abundance of the MAG --------------------------------------
shared_pairs = [p for p in pairs if p["ani"] >= ANI_SHARED and p["af_query"] >= AF_MIN]
needed = {p["mag_rep"] for p in shared_pairs}
ab = {}
with open(ABUND) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    n_samples = len(hdr) - 1
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[0] not in needed:
            continue
        v = [float(x) for x in parts[1:]]
        nz = [x for x in v if x > 0]
        ab[parts[0]] = {
            "median_all": statistics.median(v),
            "median_nonzero": statistics.median(nz) if nz else 0.0,
            "max": max(v), "n_gt0": len(nz),
        }
print(f"abundance file: {n_samples} samples; MAG reps in shared pairs with an abundance row: "
      f"{len(ab)} / {len(needed)}")

srows = []
for p in sorted(shared_pairs, key=lambda x: (-x["ani"], x["cultured"])):
    c = cult_by_id[p["cultured"]]
    a = ab.get(p["mag_rep"])
    srows.append({
        "cultured_genome": p["cultured"], "batch": c["batch"], "site": c["site"],
        "quality_tier": c["quality_tier"], "cultured_gtdb": c["cultured_gtdb"],
        "mag_rep": p["mag_rep"], "mag_gtdb": mag_tax[p["mag_rep"]],
        "ani": f"{p['ani']:.2f}", "af_query": f"{p['af_query']:.2f}", "af_ref": f"{p['af_ref']:.2f}",
        "mag_median_rel_abund_all_samples": f"{a['median_all']:.6f}" if a else "NA",
        "mag_median_rel_abund_nonzero": f"{a['median_nonzero']:.6f}" if a else "NA",
        "mag_max_rel_abund": f"{a['max']:.6f}" if a else "NA",
        "mag_n_samples_gt0": a["n_gt0"] if a else "NA",
        "n_samples_total": n_samples,
    })
with open(os.path.join(P, "results/shared_species.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(srows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(srows)
print(f"shared-species pairs (ANI >= 95, AF >= 30): {len(srows)}   "
      f"distinct cultured genomes: {len({r['cultured_genome'] for r in srows})}   "
      f"distinct MAG species: {len({r['mag_rep'] for r in srows})}")

# ---- summary -------------------------------------------------------------------------------------
def block(group, label, subset):
    n = len(subset)
    d = {b: sum(1 for r in subset if r["bin"] == b) for b in ("shared", "related", "none")}
    anis = [float(r["ani"]) for r in subset if r["ani"]]
    return {"group": group, "label": label, "n_genomes": n,
            "shared": d["shared"], "related": d["related"], "none": d["none"],
            "pct_shared": f"{100*d['shared']/n:.1f}" if n else "0.0",
            "pct_related": f"{100*d['related']/n:.1f}" if n else "0.0",
            "pct_none": f"{100*d['none']/n:.1f}" if n else "0.0",
            "n_with_any_hit": len(anis),
            "median_best_ani_of_hits": f"{statistics.median(anis):.2f}" if anis else "NA"}

srow = [block("overall", "all cultured genomes", frows)]
for b in sorted({r["batch"] for r in frows}):
    srow.append(block("batch", b, [r for r in frows if r["batch"] == b]))
for s in sorted({r["site"] for r in frows}, key=lambda x: (x == "unassigned", x)):
    srow.append(block("site", s, [r for r in frows if r["site"] == s]))
srow.append({"group": "reverse", "label": "MAG species reps (n=5438)", "n_genomes": len(mag),
             "shared": n_mag_shared,
             "related": sum(1 for r in rrows if r["bin"] == "related"),
             "none": len(mag) - sum(1 for r in rrows if r["bin"] in ("shared", "related")),
             "pct_shared": f"{100*n_mag_shared/len(mag):.2f}",
             "pct_related": f"{100*sum(1 for r in rrows if r['bin']=='related')/len(mag):.2f}",
             "pct_none": f"{100*(len(mag)-len([r for r in rrows if r['bin'] in ('shared','related')]))/len(mag):.2f}",
             "n_with_any_hit": len(rrows), "median_best_ani_of_hits":
                 f"{statistics.median([float(r['ani']) for r in rrows]):.2f}" if rrows else "NA"})

with open(os.path.join(P, "results/summary.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(srow[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(srow)
for r in srow:
    print(f"{r['group']:8s} {r['label']:26s} n={r['n_genomes']:5d}  shared={r['shared']:4d} "
          f"({r['pct_shared']:>5s} %)  related={r['related']:4d}  none={r['none']:5d}")
print("done")
