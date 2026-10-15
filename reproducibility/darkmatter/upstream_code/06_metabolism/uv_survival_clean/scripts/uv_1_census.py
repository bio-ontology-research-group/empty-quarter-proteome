#!/usr/bin/env python3
"""Step 1 -- UV/DNA-repair marker census on the CURRENT 5,438-representative catalogue.

Reuses, verbatim, the 28-KO damage-response panel defined by
  scripts/ibex/124_uv_survival_census.sh  ->  06_metabolism/uv_survival/uv_markers.tsv
and looks it up in the per-representative KO table built by
  scripts/eco_probes/gather_rep_ko.py    ->  04_functional/refs/rep_ko.tsv   (5,438 rows)

Nothing is re-annotated. This is the like-for-like replacement of the old census, which
ran over the retired 15,469 single-sample MAG set (06_metabolism/uv_survival/).

Out (results/):
  group_prevalence_5438.tsv   per functional group: n species, % of catalogue (+ old 15,469 value)
  marker_prevalence_5438.tsv  per KO: gene, group, n species, % of catalogue
  rep_traits.tsv              per representative: phylum, per-group flags, the three focal
                              traits (photolyase / NHEJ / carotenoid), n distinct panel KOs
                              -- this is the join key for steps 2-3.
"""
import csv, os, collections

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
U    = f"{M}/06_metabolism/uv_survival_clean"
RES  = f"{U}/results"
PANEL = f"{M}/06_metabolism/uv_survival/uv_markers.tsv"       # reused verbatim
# repko_v2 (04_functional/repko_clean, 2026-08-12): 0 unannotated, 690 reps had spurious KOs shed
# (subset-representative annotation collision) and the 168 empty rows recovered / re-annotated.
REPKO = f"{M}/04_functional/repko_clean/results/rep_ko_v2.tsv"
GTDB  = [f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
         f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"]
os.makedirs(RES, exist_ok=True)

# focal traits: the variable ones that carry the compartment test (old prevalence in brackets)
# K01859 is EXCLUDED from the inherited panel: labelled "phrA photolyase/cryptochrome" in
# uv_markers.tsv but it is chalcone isomerase (EC 5.5.1.6) per kofam_ko_list. Verified 2026-08-13:
# 1 carrier in 5,438, which also carries K01669 -> zero effect on any count.
EXCLUDE    = {"K01859"}
PHOTOLYASE = {"K01669"}               # photoreactivation      (old 17.5%)
NHEJ       = {"K10979", "K01971"}     # ku + ligD              (old 78.4%)
CAROT      = {"K02291", "K10027"}     # crtB + crtI            (old 58.5%)

# old 15,469-MAG values, for the correction table (06_metabolism/uv_survival/group_prevalence.tsv)
OLD = {"base_excision_repair": 99.5, "nucleotide_excision_repair": 99.4, "SOS_translesion": 94.2,
       "dna_protection_oxidative": 91.9, "dsb_nhej": 78.4, "photoprotective_pigment": 58.5,
       "mismatch_repair": 48.8, "photoreactivation": 17.5}

# ---- panel -------------------------------------------------------------------
ko_gene, ko_group = {}, {}
for r in csv.reader(open(PANEL), delimiter="\t"):
    if len(r) >= 3 and r[0] not in EXCLUDE:
        ko_gene[r[0]], ko_group[r[0]] = r[1], r[2]
print(f"[census] panel: {len(ko_gene)} KOs across {len(set(ko_group.values()))} functional groups")

# ---- taxonomy ----------------------------------------------------------------
rep_phylum = {}
for p in GTDB:
    for r in csv.DictReader(open(p), delimiter="\t"):
        lin = r["classification"]
        rep_phylum[r["user_genome"]] = next((t[3:] for t in lin.split(";") if t.startswith("p__")), "")
print(f"[census] GTDB classifications: {len(rep_phylum)}")

# ---- per-representative KO content ------------------------------------------
rep_kos, rep_nko = {}, {}
for line in open(REPKO):
    c = line.rstrip("\n").split("\t")
    if not c or not c[0]:
        continue
    allk = set(c[1].split(",")) if len(c) > 1 and c[1] else set()
    rep_nko[c[0]] = len(allk)
    rep_kos[c[0]] = allk & set(ko_gene)                     # keep only panel KOs

# repko_v2 already carries every fix; the step-0 patch is redundant and skipped.

# 168 representatives carry NO eggNOG annotation at all (empty KO field) -- they are
# unannotated, not gene-poor. Including them in the denominator deflates every prevalence
# by 5,270/5,438 = 3.1%. Prevalence is therefore reported over the ANNOTATED catalogue.
unann = {r for r, n in rep_nko.items() if n == 0}
rep_kos = {r: k for r, k in rep_kos.items() if r not in unann}
N_ALL, N = len(rep_nko), len(rep_kos)
print(f"[census] representatives: {N_ALL} total, {len(unann)} unannotated (0 KOs) -> denominator N={N}")
missing_tax = [r for r in rep_kos if r not in rep_phylum]
if missing_tax:
    print(f"[census] WARNING {len(missing_tax)} reps without GTDB classification, e.g. {missing_tax[:3]}")

# ---- prevalence --------------------------------------------------------------
group_reps  = collections.defaultdict(set)
marker_reps = collections.defaultdict(set)
for rep, kos in rep_kos.items():
    for ko in kos:
        marker_reps[ko].add(rep)
        group_reps[ko_group[ko]].add(rep)

with open(f"{RES}/group_prevalence_5438.tsv", "w") as o:
    o.write("functional_group\tn_species\tpct_of_catalogue\told_pct_15469\tdelta_pp\n")
    print(f"\n{'group':28s} {'n':>6s} {'% of 5,438':>11s}   {'old % (15,469)':>14s}")
    for g, s in sorted(group_reps.items(), key=lambda kv: -len(kv[1])):
        pct = len(s) / N * 100
        old = OLD.get(g, float("nan"))
        o.write(f"{g}\t{len(s)}\t{pct:.1f}\t{old:.1f}\t{pct-old:+.1f}\n")
        print(f"{g:28s} {len(s):6d} {pct:10.1f}%   {old:13.1f}%")

with open(f"{RES}/marker_prevalence_5438.tsv", "w") as o:
    o.write("KO\tgene\tfunctional_group\tn_species\tpct_of_catalogue\n")
    for ko, s in sorted(marker_reps.items(), key=lambda kv: -len(kv[1])):
        o.write(f"{ko}\t{ko_gene[ko]}\t{ko_group[ko]}\t{len(s)}\t{len(s)/N*100:.1f}\n")

# any marker at all
any_marker = sum(1 for kos in rep_kos.values() if kos)
print(f"\n[census] species carrying >=1 panel marker: {any_marker}/{N} ({any_marker/N*100:.1f}%)")
print(f"[census] (on the raw 5,438 denominator this would read {any_marker/N_ALL*100:.1f}% -- the "
      f"{len(unann)} unannotated reps are excluded, not counted as negative)")

# ---- per-representative trait table (join key for steps 2-3) -----------------
# all 5,438 rows with an explicit `annotated` flag: downstream steps must treat the 168
# unannotated representatives as UNKNOWN (drop from numerator AND denominator), never as negative.
groups = sorted(group_reps)
with open(f"{RES}/rep_traits.tsv", "w") as o:
    o.write("rep\tannotated\tphylum\tn_panel_kos\tphotolyase\tnhej\tcarotenoid\t" + "\t".join(groups) + "\n")
    for rep in sorted(rep_nko):
        ann = rep not in unann
        kos = rep_kos.get(rep, set())
        if ann:
            flags = [str(int(bool(kos & PHOTOLYASE))), str(int(bool(kos & NHEJ))), str(int(bool(kos & CAROT)))]
            flags += [str(int(any(ko_group[k] == g for k in kos))) for g in groups]
        else:
            flags = ["NA"] * (3 + len(groups))
        o.write(f"{rep}\t{int(ann)}\t{rep_phylum.get(rep,'')}\t{len(kos) if ann else 'NA'}\t"
                + "\t".join(flags) + "\n")

for name, s in (("photolyase", PHOTOLYASE), ("NHEJ", NHEJ), ("carotenoid", CAROT)):
    n = sum(1 for kos in rep_kos.values() if kos & s)
    print(f"[census] focal trait {name:11s}: {n:5d} species ({n/N*100:.1f}%)")
print(f"\n[census] wrote {RES}/{{group_prevalence_5438,marker_prevalence_5438,rep_traits}}.tsv")
