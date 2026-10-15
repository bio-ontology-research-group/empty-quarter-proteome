#!/bin/bash
# fam_9 — phylum attribution of the admitted dark families (E4 successor; §1.4 display support).
# Materializes the ad-hoc join run 2026-08-14 (results/fam_phylum.tsv): member protein → contig →
# refined bin (contig_to_bin_refined) → dRep species cluster (Cdb) → winner rep (Wdb) → GTDB R232
# phylum (05_refine/_gtdbtk_r232). No name translation needed — all current-style names.
# Rerunnable; overwrites results/fam_phylum.tsv atomically.
#   sbatch scripts/fam_9_phylum.sh   (light: <10 min, <8G)
#
#SBATCH --job-name=ale_fam_phylum
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam9_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
set -euo pipefail
python3 - <<'PY'
import gzip, csv, os
M = "/ibex/scratch/projects/c2014/rob/metagenomics"
B = f"{M}/06_darkmatter/darkfam_basis_clean"
G = f"{M}/05_refine/_gtdbtk_r232"

cdb = {}
for row in csv.DictReader(open(f"{M}/05_refine/_drep/out/data_tables/Cdb.csv")):
    cdb[row["genome"][:-3]] = row["secondary_cluster"]
win = {}
for row in csv.DictReader(open(f"{M}/05_refine/_drep/out/data_tables/Wdb.csv")):
    win[row["cluster"]] = row["genome"][:-3]
tax = {}
for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    for row in csv.DictReader(open(f"{G}/{f}"), delimiter="\t"):
        phy = [x for x in row["classification"].split(";") if x.startswith("p__")]
        tax[row["user_genome"]] = phy[0][3:] if phy else "unclassified"
assert len(tax) == 5438, f"taxonomy rows {len(tax)} != 5438"
print(f"[fam9] bins {len(cdb)}, clusters {len(win)}, taxonomy {len(tax)}", flush=True)

c2p = {}
with gzip.open(f"{M}/04_protein_catalog/contig_to_bin_refined.tsv.gz", "rt") as fh:
    for line in fh:
        c, _, b = line.rstrip("\n").split("\t")
        cl = cdb.get(b)
        if cl:
            p = tax.get(win.get(cl, ""), None)
            if p: c2p.setdefault(c, set()).add(p)

fam_phy, n_nomap = {}, 0
for line in open(f"{B}/results/member_fam_final.tsv"):
    m, f = line.rstrip("\n").split("\t")
    ps = c2p.get(m[:m.rfind("_")])
    if ps: fam_phy.setdefault(f, set()).update(ps)
    else: n_nomap += 1

with open(f"{B}/results/fam_phylum.tsv.part", "w") as out:
    out.write("family_rep\tn_phyla\tphyla\n")
    for f, v in sorted(fam_phy.items()):
        out.write(f"{f}\t{len(v)}\t{','.join(sorted(v))}\n")
os.replace(f"{B}/results/fam_phylum.tsv.part", f"{B}/results/fam_phylum.tsv")
single = sum(1 for v in fam_phy.values() if len(v) == 1)
print(f"[fam9] attributed {len(fam_phy)} families; members w/o refined map {n_nomap}")
print(f"[fam9] single-phylum {single} ({100*single/len(fam_phy):.1f}%)")
PY
echo "[$(date -Iseconds)] fam_9 done"
