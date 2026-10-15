#!/usr/bin/env python3
"""Build the query/reference lists and the metadata joins for the RAK-PR vs soil-MAG skani test.

Writes into P/results and P/scripts only. E and W are read only.

Outputs
  scripts/cultured_paths.txt   one absolute FASTA path per catalogue genome
  scripts/mag_paths.txt        one absolute FASTA path per soil MAG species representative
  results/cultured_meta.tsv    genome_id, batch, source_sample, site, quality_tier, gtdb_taxonomy, path
  results/mag_meta.tsv         mag_rep, gtdb_taxonomy (r232), domain
"""
import csv, os, re, sys

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
E = "/ibex/scratch/projects/c2014/alelopezv/eqmgc"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"

CATALOGUE = os.path.join(E, "catalogue.tsv")
INPUTS = os.path.join(E, "inputs.tsv")
META = [os.path.join(E, "data/metadata", f) for f in
        ("metadata_enrichments.tsv", "metadata_isolates.tsv", "metadata_thermophilic.tsv")]
MAGDIR = os.path.join(W, "05_refine/_drep/out/dereplicated_genomes")
GTDB_BAC = os.path.join(W, "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv")
GTDB_AR = os.path.join(W, "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv")


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---- 1. sample -> site, from the three RAK-PR metadata tables -------------------------------
site_of = {}
for m in META:
    for row in read_tsv(m):
        hit = re.search(r"Site\s+(\d+)", row.get("geo_loc_name", "") or "")
        if hit:
            site_of[row["sample_name"]] = hit.group(1)
print(f"metadata samples with a site: {len(site_of)}")

# ---- 2. genome_id -> FASTA path, from inputs.tsv ---------------------------------------------
path_of = {r["genome_id"]: r["genome_path"] for r in read_tsv(INPUTS)}
print(f"inputs.tsv rows: {len(path_of)}")

# ---- 3. the 722 catalogue genomes -------------------------------------------------------------
cat = read_tsv(CATALOGUE)
print(f"catalogue genomes: {len(cat)}")

rows, missing = [], []
for r in cat:
    gid = r["genome_id"]
    # single-bin enrichments are registered in inputs.tsv as <sample>.1 but named <sample> in the
    # catalogue (verified: identical CheckM2/BUSCO values for CO-7/CO-7.1, Ht-24, Ht-35, NaB-5, NaB-P-2)
    p = path_of.get(gid) or path_of.get(gid + ".1")
    if not p or not os.path.isfile(p):
        missing.append(gid)
        continue
    rows.append({
        "genome_id": gid,
        "batch": r["batch"],
        "source_sample": r["source_sample"],
        "site": site_of.get(r["source_sample"], "unassigned"),
        "quality_tier": r["quality_tier"],
        "cultured_gtdb": r["gtdb_taxonomy"],
        "path": p,
    })

print(f"FASTAs located: {len(rows)} / {len(cat)}   missing: {missing}")
by_batch = {}
for r in rows:
    by_batch[r["batch"]] = by_batch.get(r["batch"], 0) + 1
print("per batch:", by_batch)
by_site = {}
for r in rows:
    by_site[r["site"]] = by_site.get(r["site"], 0) + 1
print("per site:", by_site)

os.makedirs(os.path.join(P, "results"), exist_ok=True)
with open(os.path.join(P, "results/cultured_meta.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
with open(os.path.join(P, "scripts/cultured_paths.txt"), "w") as f:
    for r in rows:
        f.write(r["path"] + "\n")

# ---- 4. the soil MAG species representatives ---------------------------------------------------
mags = sorted(x for x in os.listdir(MAGDIR) if x.endswith(".fa"))
print(f"MAG species representatives: {len(mags)}")
with open(os.path.join(P, "scripts/mag_paths.txt"), "w") as f:
    for m in mags:
        f.write(os.path.join(MAGDIR, m) + "\n")

tax_of, dom_of = {}, {}
for src, dom in ((GTDB_BAC, "Bacteria"), (GTDB_AR, "Archaea")):
    for r in read_tsv(src):
        tax_of[r["user_genome"]] = r["classification"]
        dom_of[r["user_genome"]] = dom
print(f"gtdbtk r232 rows: {len(tax_of)}")

mrows, no_tax = [], 0
for m in mags:
    stem = m[:-3]                       # strip .fa
    t = tax_of.get(stem) or tax_of.get(m)
    if t is None:
        no_tax += 1
        t = "unclassified"
    mrows.append({"mag_rep": stem, "mag_gtdb": t,
                  "domain": dom_of.get(stem, dom_of.get(m, "NA"))})
print(f"MAG reps with no gtdbtk row: {no_tax}")
with open(os.path.join(P, "results/mag_meta.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["mag_rep", "mag_gtdb", "domain"], delimiter="\t",
                       lineterminator="\n")
    w.writeheader(); w.writerows(mrows)

print("done")
