#!/bin/bash
# REBUILD step 2 — recurrence / persistence tables (§1.4 legs 1+4) + final family admission.
# One streaming pass over all_cluster90.tsv (753M rows): for every dark cluster90 rep, collect the
# sample set of ALL its protein instances; aggregate to family via fam_cluster.tsv. Sample→site/trip
# from results/sample_site_trip.tsv (built from the sequencing-results directory tree; 551 covered).
# Admission: >=3 diverged members AND >=3 samples AND >=2 sites. Paralog-only clusters are dropped;
# recurrent single-variant dark (1-2 members, >=3 samples/>=2 sites) reported as a named class.
#
# Submit: sbatch scripts/fam_2_recurrence.sh
#
#SBATCH --job-name=ale_fam_recur
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam2_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
CLU=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/all_linclust90/all_cluster90.tsv
python3 - "$B" "$CLU" <<'PY'
import sys, gzip
from collections import defaultdict
B, CLU = sys.argv[1], sys.argv[2]
R = B + "/results"

# sample -> (site_idx, trip_idx)
samp_meta, sites, trips = {}, {}, {}
for line in open(R + "/sample_site_trip.tsv"):
    s, site, trip = line.rstrip("\n").split("\t")
    samp_meta[s] = (sites.setdefault(site, len(sites)), trips.setdefault(trip, len(trips)))
samp_idx = {}

# FIX 2026-08-13 (review finding): restrict instance counting to the 551 canonical samples.
# The catalogue holds 3 non-canonical assemblies (48PRR1, T8DDr2, T8DDr3) whose instances
# inflated n_samples (no site/trip bits) and could wrongly satisfy the >=3-samples gate.
canon = set(l.strip() for l in open(
    "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology/canon_samples.txt"))
print(f"[step2] canonical filter: {len(canon)} samples", flush=True)

# dark rep -> family index; family member counts
fam_idx, fam_members = {}, []
rep_fam = {}
for line in open(R + "/fam_cluster.tsv"):
    f, m = line.rstrip("\n").split("\t")
    fi = fam_idx.get(f)
    if fi is None:
        fi = fam_idx[f] = len(fam_members); fam_members.append(0)
    fam_members[fi] += 1
    rep_fam[m] = fi
nfam = len(fam_members)
print(f"[step2] families(all clusters): {nfam}  dark reps: {len(rep_fam)}", flush=True)

fam_s = [0]*nfam; fam_site = [0]*nfam; fam_trip = [0]*nfam; fam_inst = [0]*nfam
rep_ids = list(rep_fam.keys())
rep_pos = {r: i for i, r in enumerate(rep_ids)}
rep_s = [0]*len(rep_ids); rep_inst = [0]*len(rep_ids)

n = 0
with open(CLU) as fh:
    for line in fh:
        rep, _, member = line.partition("\t")
        n += 1
        if n % 100000000 == 0:
            print(f"[step2] {n//1000000}M rows", flush=True)
        fi = rep_fam.get(rep)
        if fi is None:
            continue
        member = member.rstrip("\n")
        s = member[:member.find("::")]
        if s not in canon:
            continue
        si = samp_idx.setdefault(s, len(samp_idx))
        bit = 1 << si
        fam_s[fi] |= bit; fam_inst[fi] += 1
        p = rep_pos[rep]
        rep_s[p] |= bit; rep_inst[p] += 1
        meta = samp_meta.get(s)
        if meta:
            fam_site[fi] |= 1 << meta[0]
            fam_trip[fi] |= 1 << meta[1]

print(f"[step2] pass done: {n} rows; distinct samples seen: {len(samp_idx)}", flush=True)

inv_fam = {v: k for k, v in fam_idx.items()}
adm, par_only, young = 0, 0, 0
with open(R + "/fam_recurrence.tsv", "w") as out, open(R + "/FAMILIES_FINAL.txt", "w") as fin:
    out.write("family_rep\tn_members\tn_instances\tn_samples\tn_sites\tn_trips\n")
    for fi in range(nfam):
        ns, nsite, ntrip = bin(fam_s[fi]).count("1"), bin(fam_site[fi]).count("1"), bin(fam_trip[fi]).count("1")
        out.write(f"{inv_fam[fi]}\t{fam_members[fi]}\t{fam_inst[fi]}\t{ns}\t{nsite}\t{ntrip}\n")
        if fam_members[fi] >= 3:
            if ns >= 3 and nsite >= 2:
                adm += 1; fin.write(inv_fam[fi] + "\n")
            else:
                par_only += 1
        elif ns >= 3 and nsite >= 2:
            young += 1
with gzip.open(R + "/darkrep_recurrence.tsv.gz", "wt") as out:
    out.write("dark_rep\tn_instances\tn_samples\n")
    for i, r in enumerate(rep_ids):
        out.write(f"{r}\t{rep_inst[i]}\t{bin(rep_s[i]).count('1')}\n")
print(f"[step2] ADMITTED families (>=3 members & >=3 samples & >=2 sites): {adm}")
print(f"[step2] dropped >=3-member clusters failing instance gate (paralog-only): {par_only}")
print(f"[step2] recurrent single/double-variant dark (young class): {young}")
print("[step2] done")
PY
