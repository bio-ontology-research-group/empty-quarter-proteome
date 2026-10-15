#!/usr/bin/env python3
"""RQ-F2 -- BiG-SLICE GCF novelty analysis (runs on IBEX after clustering).

Reads a BiG-SLICE v2 result database (data.db) produced by co-clustering all
Empty Quarter antiSMASH-8 BGC regions (dataset 'eq') with MIBiG reference BGCs
(dataset 'mibig') into Gene Cluster Families (GCFs; BIRCH centroids).

Novelty rule
------------
A GCF is *novel* (no characterized analog) if it contains >=1 EQ BGC and 0
MIBiG BGCs. The headline is the fraction of EQ-containing GCFs that are novel.

Outputs (into <outdir>):
  result_summary.json / result_summary.tsv  -- headline numbers
  gcf_composition.tsv                        -- per-GCF eq/mibig counts + novel flag
  novelty_by_chemclass.tsv                   -- by broad biosynthetic class
  novelty_by_product.tsv                     -- by antiSMASH product (chem_subclass)

Usage: python 87_bgc_novelty_analyze.py <data.db> <outdir>
"""
import sqlite3
import sys
import os
import json
import csv
from collections import defaultdict

db_path, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
con = sqlite3.connect(db_path)
con.row_factory = sqlite3.Row
cur = con.cursor()

# --- latest clustering run ---
cl = cur.execute(
    "SELECT id, threshold, num_centroids, clustering_method "
    "FROM clustering ORDER BY id DESC LIMIT 1").fetchone()
clid = cl["id"]
threshold = cl["threshold"]
num_centroids = cl["num_centroids"]
method = cl["clustering_method"]

# --- dataset ids: reference = any dataset whose name contains 'mibig'; rest = EQ ---
ds = {r["name"]: r["id"] for r in cur.execute("SELECT id, name FROM dataset")}
mibig_ids = {i for n, i in ds.items() if "mibig" in n.lower()}
eq_ids = {i for n, i in ds.items() if i not in mibig_ids}
eq_names = [n for n, i in ds.items() if i in eq_ids]
mibig_names = [n for n, i in ds.items() if i in mibig_ids]

# --- per-dataset totals + contig-edge ---
totals = {}
for r in cur.execute(
        "SELECT d.name name, COUNT(*) c, "
        "SUM(CASE WHEN b.on_contig_edge=1 THEN 1 ELSE 0 END) edge "
        "FROM bgc b JOIN dataset d ON d.id=b.dataset_id GROUP BY d.name"):
    totals[r["name"]] = {"n_bgc": r["c"], "n_contig_edge": r["edge"] or 0}

# --- assignment: bgc -> best-rank GCF for this clustering ---
assign = {}  # bgc_id -> (gcf_id, dataset_id, rank)
for r in cur.execute(
        "SELECT m.bgc_id bid, m.gcf_id gid, m.rank rk, b.dataset_id ds "
        "FROM gcf_membership m "
        "JOIN gcf g ON g.id=m.gcf_id AND g.clustering_id=? "
        "JOIN bgc b ON b.id=m.bgc_id", (clid,)):
    bid = r["bid"]
    if bid not in assign or r["rk"] < assign[bid][2]:
        assign[bid] = (r["gid"], r["ds"], r["rk"])

# --- per-GCF composition ---
gcf_eq = defaultdict(int)
gcf_mibig = defaultdict(int)
for bid, (gid, dsid, rk) in assign.items():
    if dsid in eq_ids:
        gcf_eq[gid] += 1
    elif dsid in mibig_ids:
        gcf_mibig[gid] += 1

all_gcfs = set(gcf_eq) | set(gcf_mibig)
eq_gcfs = [g for g in all_gcfs if gcf_eq[g] > 0]
novel_gcfs = set(g for g in eq_gcfs if gcf_mibig[g] == 0)
shared_gcfs = set(g for g in eq_gcfs if gcf_mibig[g] > 0)

n_eq_gcf = len(eq_gcfs)
n_novel_gcf = len(novel_gcfs)
n_shared_gcf = len(shared_gcfs)
pct_novel_gcf = 100.0 * n_novel_gcf / n_eq_gcf if n_eq_gcf else 0.0

# --- EQ BGCs in novel vs shared GCFs ---
eq_assign = {b: g for b, (g, d, rk) in assign.items() if d in eq_ids}
eq_bgc_total = len(eq_assign)
eq_bgc_in_novel = sum(1 for g in eq_assign.values() if g in novel_gcfs)
pct_eq_bgc_novel = 100.0 * eq_bgc_in_novel / eq_bgc_total if eq_bgc_total else 0.0

# --- products / superclasses per EQ bgc ---
prod = defaultdict(set)   # bgc_id -> {(subclass, class)}
eq_id_list = ",".join(str(i) for i in sorted(eq_ids))
for r in cur.execute(
        "SELECT bc.bgc_id bid, cs.name sub, cc.name cls "
        "FROM bgc_class bc "
        "JOIN chem_subclass cs ON cs.id=bc.chem_subclass_id "
        "JOIN chem_class cc ON cc.id=cs.class_id "
        "JOIN bgc b ON b.id=bc.bgc_id "
        "WHERE b.dataset_id IN (" + eq_id_list + ")"):
    prod[r["bid"]].add((r["sub"], r["cls"]))

# by broad class and by product (subclass)
class_bgc = defaultdict(set)       # class -> set(bgc_id)
class_bgc_novel = defaultdict(set)
class_gcf = defaultdict(set)       # class -> set(gcf_id)
class_gcf_novel = defaultdict(set)
sub_bgc = defaultdict(set)
sub_bgc_novel = defaultdict(set)
sub_gcf = defaultdict(set)
sub_gcf_novel = defaultdict(set)

for bid, gid in eq_assign.items():
    is_novel = gid in novel_gcfs
    for sub, cls in prod.get(bid, {("Unknown", "Unknown")}):
        class_bgc[cls].add(bid)
        class_gcf[cls].add(gid)
        sub_bgc[sub].add(bid)
        sub_gcf[sub].add(gid)
        if is_novel:
            class_bgc_novel[cls].add(bid)
            class_gcf_novel[cls].add(gid)
            sub_bgc_novel[sub].add(bid)
            sub_gcf_novel[sub].add(gid)


def write_class_table(path, keys, bgc, bgc_novel, gcf, gcf_novel):
    rows = []
    for k in keys:
        nb = len(bgc[k])
        nbn = len(bgc_novel[k])
        ng = len(gcf[k])
        ngn = len(gcf_novel[k])
        rows.append({
            "class": k,
            "n_eq_bgc": nb,
            "n_eq_bgc_in_novel_gcf": nbn,
            "pct_eq_bgc_novel": round(100.0 * nbn / nb, 2) if nb else 0.0,
            "n_gcf": ng,
            "n_novel_gcf": ngn,
            "pct_gcf_novel": round(100.0 * ngn / ng, 2) if ng else 0.0,
        })
    rows.sort(key=lambda x: -x["n_eq_bgc"])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["class"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


cls_rows = write_class_table(
    os.path.join(outdir, "novelty_by_chemclass.tsv"),
    class_bgc.keys(), class_bgc, class_bgc_novel, class_gcf, class_gcf_novel)
write_class_table(
    os.path.join(outdir, "novelty_by_product.tsv"),
    sub_bgc.keys(), sub_bgc, sub_bgc_novel, sub_gcf, sub_gcf_novel)

# --- per-GCF composition table ---
with open(os.path.join(outdir, "gcf_composition.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gcf_id", "n_eq_bgc", "n_mibig_bgc", "novel"])
    for g in sorted(all_gcfs):
        w.writerow([g, gcf_eq[g], gcf_mibig[g],
                    1 if (gcf_eq[g] > 0 and gcf_mibig[g] == 0) else 0])

# --- headline summary ---
n_eq_bgc_tot = sum(totals.get(n, {}).get("n_bgc", 0) for n in eq_names)
n_eq_edge_tot = sum(totals.get(n, {}).get("n_contig_edge", 0) for n in eq_names)
n_mibig_tot = sum(totals.get(n, {}).get("n_bgc", 0) for n in mibig_names)
summary = {
    "clustering_id": clid,
    "clustering_method": method,
    "threshold": threshold,
    "num_centroids_total": num_centroids,
    "eq_datasets": eq_names,
    "mibig_datasets": mibig_names,
    "n_eq_bgc": n_eq_bgc_tot,
    "n_eq_bgc_contig_edge": n_eq_edge_tot,
    "n_mibig_bgc": n_mibig_tot,
    "n_eq_bgc_clustered": eq_bgc_total,
    "n_gcf_eq": n_eq_gcf,
    "n_gcf_novel": n_novel_gcf,
    "n_gcf_shared_with_mibig": n_shared_gcf,
    "pct_gcf_novel": round(pct_novel_gcf, 2),
    "n_eq_bgc_in_novel_gcf": eq_bgc_in_novel,
    "pct_eq_bgc_in_novel_gcf": round(pct_eq_bgc_novel, 2),
}
with open(os.path.join(outdir, "result_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
with open(os.path.join(outdir, "result_summary.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    for k, v in summary.items():
        w.writerow([k, v])

print(json.dumps(summary, indent=2))
print("\nTop classes by EQ BGC count:")
for r in cls_rows[:12]:
    print("  {class:24s} n_bgc={n_eq_bgc:7d}  novel_gcf={pct_gcf_novel:6.2f}%  "
          "bgc_in_novel={pct_eq_bgc_novel:6.2f}%".format(**r))
con.close()
