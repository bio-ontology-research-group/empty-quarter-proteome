#!/bin/bash
# darkclean_5 (rewrite) — restrict the catalogue to proteins that map to an HQ+MQ MAG (genome-resolved,
# BOTH binnings): the clean subset that replaces the OLD single-sample ladder 51.9M proteins -> 23.4M reps
# -> 2.84M UniRef50-dark. Rationale: the 753M assembly proteins include unbinned/fragment/noise contigs;
# only proteins on a binned (HQ+MQ) contig are genome-resolved.
#
# METHOD (the two properties differ):
#   * MAG-membership = PER-PROTEIN -> checked on EVERY member (a cluster's members sit on different contigs,
#     some binned, some not; the rep alone is not representative).
#   * dark = CLUSTER property -> taken from the rep (90/80 clusters = same family; only reps were UniRef50-searched).
# A cluster is "MAG-mapped" if ANY member is on an HQ+MQ-binned contig. One member-level pass over
# all_cluster90.tsv (grouped by rep) yields BOTH:
#   (1) mag_dark_summary.tsv   — clean global totals (<-> 51.9M / 23.4M / 2.84M)
#   (2) per_bin_dark.tsv.gz     — per-MAG (per-species = is_rep subset) total/dark/fraction  [free in the same pass]
#
# Submit (after darkclean_4 + Phase-1 provenance table): sbatch .../darkmatter_clean/scripts/darkclean_5_mag.sh
#
#SBATCH --job-name=velazqal_darkclean5
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean5_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean5_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate instrain

python3 - <<'PY'
import gzip, csv, os
from collections import defaultdict
M="/ibex/scratch/projects/c2014/rob/metagenomics"
CAT=f"{M}/04_protein_catalog"
OUT=f"{CAT}/darkmatter_clean/05_mag_dark"; os.makedirs(OUT, exist_ok=True)
PROV=f"{M}/results/catalogue_provenance.tsv"                         # Phase-1: 37,949 HQ+MQ bins + rep_genome
C2B=f"{CAT}/contig_to_bin.tsv.gz"                                    # unrefined: BOTH binnings
DARK=f"{CAT}/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz"   # 69.9M dark protein reps
CLUSTER=f"{CAT}/all_linclust90/all_cluster90.tsv"                    # 753M rep<TAB>member, grouped by rep

def strip(g): return g[:-6] if g.endswith(".fa.gz") else g

# HQ+MQ bins + winner (is_rep) set, from the Phase-1 provenance table
hqmq=set(); winners=set()
with open(PROV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        g=strip(row["genome"]); hqmq.add(g)
        if g==strip(row["rep_genome"]): winners.add(g)
print(f"[info] HQ+MQ bins {len(hqmq):,}; winners {len(winners):,}", flush=True)

# contig -> HQ+MQ bin(s) (a contig can be in one single-sample AND one self-supervised bin)
pref={"single_sample":"ss","self_supervised":"selfsup"}
c2b=defaultdict(list)
with gzip.open(C2B,"rt") as fh:
    for line in fh:
        c,strat,sb=line.rstrip("\n").split("\t")
        b=f"{pref[strat]}__{c.split('::')[0]}__{sb}"
        if b in hqmq: c2b[c].append(b)
print(f"[info] HQ+MQ-binned contigs {len(c2b):,}", flush=True)

# dark protein-rep set
dark=set()
with gzip.open(DARK,"rt") as fh:
    for line in fh:
        v=line.rstrip("\n")
        if v: dark.add(v)
print(f"[info] dark reps {len(dark):,}", flush=True)

# one member-level pass, grouped by rep
tot=defaultdict(int); drk=defaultdict(int)          # per-bin (view 2)
mag_prot=mag_dark_prot=mag_reps=mag_dark_reps=0      # global (view 1)
cur=None; cur_dark=False; cur_in_mag=False; n=0
dark_ids=gzip.open(f"{OUT}/mag_dark_rep_ids.txt.gz","wt")   # genome-resolved dark rep IDs (for 00_structure_clean)
with open(CLUSTER) as fh:
    for line in fh:
        n+=1
        if n % 100_000_000 == 0: print(f"[info] {n:,} members scanned", flush=True)
        t=line.find("\t")
        if t<0: continue
        rep=line[:t]; member=line[t+1:].rstrip("\n")
        if rep!=cur:
            if cur is not None and cur_in_mag:
                mag_reps+=1
                if cur_dark:
                    mag_dark_reps+=1
                    dark_ids.write(cur+"\n")
            cur=rep; cur_dark=(rep in dark); cur_in_mag=False
        bins=c2b.get(member.rsplit("_",1)[0])        # protein's contig -> HQ+MQ bin(s)
        if bins:
            mag_prot+=1                               # unique protein (counted once)
            if cur_dark: mag_dark_prot+=1
            cur_in_mag=True
            for b in bins:                            # per-bin (a shared protein counts in each genome)
                tot[b]+=1
                if cur_dark: drk[b]+=1
    if cur is not None and cur_in_mag:               # close last cluster
        mag_reps+=1
        if cur_dark:
            mag_dark_reps+=1
            dark_ids.write(cur+"\n")
dark_ids.close()

# (1) clean global totals
with open(f"{OUT}/mag_dark_summary.tsv","w") as s:
    s.write("metric\tvalue\told_ladder_analog\n")
    s.write(f"mag_mapped_proteins\t{mag_prot}\t51.9M\n")
    s.write(f"mag_mapped_reps\t{mag_reps}\t23.4M\n")
    s.write(f"mag_mapped_dark_reps\t{mag_dark_reps}\t2.84M\n")
    s.write(f"dark_fraction_reps\t{mag_dark_reps/mag_reps if mag_reps else 0:.6f}\t-\n")
    s.write(f"mag_mapped_dark_proteins\t{mag_dark_prot}\t-\n")
    s.write(f"dark_fraction_proteins\t{mag_dark_prot/mag_prot if mag_prot else 0:.6f}\t-\n")

# (2) per-bin / per-species census (free in the same pass)
with gzip.open(f"{OUT}/per_bin_dark.tsv.gz","wt") as o:
    o.write("genome\tis_rep\ttotal_proteins\tdark_proteins\tdark_fraction\n")
    for b in sorted(tot):
        tt=tot[b]; dd=drk[b]
        o.write(f"{b}\t{int(b in winners)}\t{tt}\t{dd}\t{dd/tt if tt else 0:.6f}\n")

print(f"[info] MAG-mapped: {mag_prot:,} proteins / {mag_reps:,} reps; "
      f"dark reps {mag_dark_reps:,} ({100*mag_dark_reps/max(mag_reps,1):.2f}%); "
      f"per-bin rows {len(tot):,}", flush=True)
print("DONE_DARKCLEAN5")
PY
