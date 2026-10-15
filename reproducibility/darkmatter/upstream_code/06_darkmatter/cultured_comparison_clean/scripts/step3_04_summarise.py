#!/usr/bin/env python3
"""step3_04 - build the two step-three tables from the finished searches.

  results/step3_dark_with_cultured_homolog.tsv
      how many of the 4,292,446 MAG-dark reps have a homolog in a cultured genome,
      overall, by best-hit identity bin, by shared-species status of the cultured
      genome, and restricted to the 1,466,943 reps >= 100 aa.

  results/step3_cultured_dark_rate.tsv
      per cultured genome, the dark rate under the MAG rule (no UniRef50 hit at
      -s 7.5 -e 1e-3 --max-seqs 10), with the matching MAG rep's dark fraction.
"""
import glob
import gzip
import os
from collections import defaultdict

from scipy.stats import binomtest, wilcoxon

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
E = "/ibex/scratch/projects/c2014/alelopezv/eqmgc"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
LIN = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/results"
# UniRef50 hit tables, one per query chunk, written by step3_03 in this folder.
UNIREF_HITS = f"{P}/tmp/uniref_hits"

N_DARK_REPS = 4292446
GE100_IDS = f"{W}/04_structure/00_structure_clean/03_search/target_ids.txt"

# ---------------------------------------------------------------- shared species
shared_mag = {}          # cultured genome -> (mag_rep, ani) best ANI
with open(f"{P}/results/shared_species.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ic, im, ia = h.index("cultured_genome"), h.index("mag_rep"), h.index("ani")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        g, m, a = f[ic], f[im], float(f[ia])
        if g not in shared_mag or a > shared_mag[g][1]:
            shared_mag[g] = (m, a)
print(f"[info] cultured genomes sharing a species with a MAG: {len(shared_mag)}")

with open(f"{P}/results/step3_cultured_counts.tsv") as fh:
    n_cultured_genomes = sum(1 for _ in fh) - 1
print(f"[info] cultured genomes with a proteome: {n_cultured_genomes}")

# ------------------------------------------------------------ part 1: DIAMOND
# DIAMOND ran without a coverage filter; the reciprocal >= 50 % rule is applied here so the
# relaxed count (no coverage requirement) stays available as a sensitivity.
COV = 50.0
best_pid = {}            # dark rep -> best pident among reciprocal-50 % cultured hits
hit_shared = set()       # dark reps with >=1 such hit in a shared-species cultured genome
hit_any = set()          # dark reps with >=1 reciprocal-50 % cultured hit
hit_relaxed = set()      # dark reps with >=1 cultured hit at e-5, any coverage
hit_shared_ge90 = set()  # ... at >= 90 % identity: same species, same protein
cult_subject_hit = set() # cultured proteins recovered as a MAG-dark rep (reciprocal 50 %)
cult_genomes_hit = set()
n_rows = 0
for m8 in sorted(glob.glob(f"{P}/tmp/dia_hits/hits_[0-9]*.m8")):
    with open(m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s, pid = f[0], f[1], float(f[2])
            qcov, scov = float(f[8]), float(f[9])
            n_rows += 1
            hit_relaxed.add(q)
            if qcov < COV or scov < COV:
                continue
            hit_any.add(q)
            if q not in best_pid or pid > best_pid[q]:
                best_pid[q] = pid
            cg = s.split("@")[0]
            cult_genomes_hit.add(cg)
            cult_subject_hit.add(s)
            if cg in shared_mag:
                hit_shared.add(q)
                if pid >= 90:
                    hit_shared_ge90.add(q)
print(f"[info] diamond rows {n_rows:,}; dark reps with a cultured homolog {len(hit_any):,}")

# matched control: the 106,042 dark family reps and the 106,042 length- and size-matched
# ANNOTATED family reps, searched against the same cultured DB in the same run.
ctrl_hit = {"dark": set(), "ctrl": set()}
fam = {}
for tag, path in (("dark", f"{W}/06_darkmatter/darkfam_basis_clean/results/FAMILIES_FINAL.txt"),
                  ("ctrl", f"{W}/06_darkmatter/darkfam_basis_clean/results/CONTROL_FAMILIES.txt")):
    with open(path) as fh:
        for line in fh:
            fam[line.strip()] = tag
ctrl_m8 = f"{P}/tmp/dia_hits/hits_ctrl.m8"
if os.path.exists(ctrl_m8):
    with open(ctrl_m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if float(f[8]) < COV or float(f[9]) < COV:
                continue
            t = fam.get(f[0])
            if t:
                ctrl_hit[t].add(f[0])
n_fam = {"dark": sum(1 for v in fam.values() if v == "dark"),
         "ctrl": sum(1 for v in fam.values() if v == "ctrl")}
print(f"[info] control: dark families {len(ctrl_hit['dark'])}/{n_fam['dark']}, "
      f"annotated families {len(ctrl_hit['ctrl'])}/{n_fam['ctrl']}")

ge100 = set()
with open(GE100_IDS) as fh:
    for line in fh:
        ge100.add(line.strip())
n_ge100 = len(ge100)
hit_ge100 = sum(1 for q in hit_any if q in ge100)


def frac(n, d):
    return f"{n / d:.6f}" if d else "NA"


bins = {"lt30": 0, "30_50": 0, "50_70": 0, "70_90": 0, "ge90": 0}
for q, pid in best_pid.items():
    if pid < 30:
        bins["lt30"] += 1
    elif pid < 50:
        bins["30_50"] += 1
    elif pid < 70:
        bins["50_70"] += 1
    elif pid < 90:
        bins["70_90"] += 1
    else:
        bins["ge90"] += 1

rows = [
    ("mag_dark_reps_total", N_DARK_REPS, N_DARK_REPS),
    ("with_cultured_homolog", len(hit_any), N_DARK_REPS),
    ("with_cultured_hit_any_coverage", len(hit_relaxed), N_DARK_REPS),
    ("best_identity_lt30", bins["lt30"], N_DARK_REPS),
    ("best_identity_30_50", bins["30_50"], N_DARK_REPS),
    ("best_identity_50_70", bins["50_70"], N_DARK_REPS),
    ("best_identity_70_90", bins["70_90"], N_DARK_REPS),
    ("best_identity_ge90", bins["ge90"], N_DARK_REPS),
    ("homolog_in_shared_species_genome", len(hit_shared), N_DARK_REPS),
    ("homolog_in_shared_species_genome_ge90pct", len(hit_shared_ge90), N_DARK_REPS),
    ("homolog_only_in_unshared_genome", len(hit_any) - len(hit_shared), N_DARK_REPS),
    ("dark_reps_ge100aa", n_ge100, n_ge100),
    ("ge100aa_with_cultured_homolog", hit_ge100, n_ge100),
    ("cultured_genomes_contributing", len(cult_genomes_hit), n_cultured_genomes),
    ("control_dark_family_reps_with_homolog", len(ctrl_hit["dark"]), n_fam["dark"]),
    ("control_annotated_family_reps_with_homolog", len(ctrl_hit["ctrl"]), n_fam["ctrl"]),
]
with open(f"{P}/results/step3_dark_with_cultured_homolog.tsv", "w") as out:
    out.write("metric\tn\tdenominator\tfraction\n")
    for k, n, d in rows:
        out.write(f"{k}\t{n}\t{d}\t{frac(n, d)}\n")

# ---------------------------------------------------- part 2: cultured dark rate
hypo = set()
with open(f"{P}/tmp/cultured_hypo_ids.txt") as fh:
    for line in fh:
        hypo.add(line.strip())

# A missing or truncated UniRef50 chunk would silently turn all of its proteins dark, so there
# must be one hit table per query chunk and every chunk must report a plausible hit rate.
n_query_chunks = len(glob.glob(f"{P}/tmp/cultured_chunks/chunk_*.faa"))
chunks = sorted(glob.glob(f"{UNIREF_HITS}/hits_chunk_*.tsv"))
assert len(chunks) == n_query_chunks, (
    f"expected {n_query_chunks} UniRef50 hit tables, one per query chunk, "
    f"found {len(chunks)} in {UNIREF_HITS}")
uniref_hit = set()
for tsv in chunks:
    before = len(uniref_hit)
    with open(tsv) as fh:
        for line in fh:
            uniref_hit.add(line.split("\t", 1)[0])
    idx = os.path.basename(tsv)[len("hits_chunk_"):-len(".tsv")]
    n_q = sum(1 for _ in open(f"{P}/tmp/cultured_chunks/chunk_{idx}.faa") if _[0] == ">")
    rate = (len(uniref_hit) - before) / n_q
    print(f"[info] {os.path.basename(tsv)}: {rate:.4f} of {n_q} queries hit")
    assert rate > 0.5, f"{tsv} hit rate {rate:.4f} is implausible - rerun that chunk"
print(f"[info] cultured proteins with a UniRef50 hit: {len(uniref_hit):,}")

n_prot = defaultdict(int)
n_dark = defaultdict(int)
n_hyp = defaultdict(int)
n_hyp_dark = defaultdict(int)
dark_ids_path = f"{P}/tmp/cultured_dark_ids.txt"
dark_ids_set = set()
with open(f"{P}/tmp/cultured_all.faa") as fh, open(dark_ids_path, "w") as dout:
    for line in fh:
        if line[0] != ">":
            continue
        pid = line[1:].split()[0]
        g = pid.split("@")[0]
        n_prot[g] += 1
        ish = pid in hypo
        if ish:
            n_hyp[g] += 1
        if pid not in uniref_hit:
            n_dark[g] += 1
            dark_ids_set.add(pid)
            dout.write(pid + "\n")
            if ish:
                n_hyp_dark[g] += 1

# cultured taxonomy + batch
meta = {}
with open(f"{E}/catalogue.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ig, ib, it = h.index("genome_id"), h.index("batch"), h.index("gtdb_taxonomy")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        tax = {p[:3]: p[3:] for p in f[it].split(";")}
        meta[f[ig]] = (f[ib], tax.get("p__", ""), tax.get("g__", ""))

# MAG rep dark fraction. v2_rep_dark.tsv is the native-key rebuild covering all 5,438 dRep
# winners (reps_dark.tsv is the earlier renamed join and lost 79 of them).
mag_dark, mag_nprot, mag_ndark = {}, {}, {}
mag_dark_src = f"{LIN}/v2_rep_dark.tsv"
if not os.path.exists(mag_dark_src):
    mag_dark_src = f"{LIN}/reps_all.tsv"
with open(mag_dark_src) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ig = h.index("rep") if "rep" in h else h.index("genome")
    idf = h.index("dark_fraction")
    inp = h.index("n_proteins") if "n_proteins" in h else None
    ind = h.index("n_dark") if "n_dark" in h else None
    for line in fh:
        f = line.rstrip("\n").split("\t")
        mag_dark[f[ig]] = f[idf]
        if inp is not None:
            mag_nprot[f[ig]] = f[inp]
        if ind is not None:
            mag_ndark[f[ig]] = f[ind]
print(f"[info] MAG rep dark fractions from {os.path.basename(mag_dark_src)}: {len(mag_dark)}")

with open(f"{P}/results/step3_cultured_dark_rate.tsv", "w") as out:
    out.write("genome\tphylum\tgenus\tbatch\tn_proteins\tn_dark\tdark_fraction\t"
              "n_hypothetical\tn_hypothetical_dark\thypothetical_dark_fraction\t"
              "shares_species_with_MAG\tmag_rep\tmag_ani\tmag_dark_fraction\n")
    for g in sorted(n_prot, key=lambda x: -n_dark[x] / max(n_prot[x], 1)):
        b, ph, ge = meta.get(g, ("NA", "NA", "NA"))
        sh = g in shared_mag
        mrep, mani = shared_mag.get(g, ("NA", float("nan")))
        mdf = mag_dark.get(mrep, "NA") if sh else "NA"
        out.write(
            f"{g}\t{ph or 'NA'}\t{ge or 'NA'}\t{b}\t{n_prot[g]}\t{n_dark[g]}\t"
            f"{n_dark[g] / n_prot[g]:.6f}\t{n_hyp[g]}\t{n_hyp_dark[g]}\t"
            f"{(n_hyp_dark[g] / n_hyp[g]) if n_hyp[g] else 0:.6f}\t"
            f"{'yes' if sh else 'no'}\t{mrep}\t"
            f"{('%.2f' % mani) if sh else 'NA'}\t{mdf}\n")

tot_p = sum(n_prot.values())
tot_d = sum(n_dark.values())
tot_h = sum(n_hyp.values())
tot_hd = sum(n_hyp_dark.values())
sh_p = sum(n_prot[g] for g in n_prot if g in shared_mag)
sh_d = sum(n_dark[g] for g in n_prot if g in shared_mag)
ns_p = tot_p - sh_p
ns_d = tot_d - sh_d
with open(f"{P}/results/step3_cultured_dark_rate_summary.tsv", "w") as out:
    out.write("set\tn_genomes\tn_proteins\tn_dark\tdark_fraction\n")
    out.write(f"all_cultured\t{len(n_prot)}\t{tot_p}\t{tot_d}\t{frac(tot_d, tot_p)}\n")
    ng_sh = sum(1 for g in n_prot if g in shared_mag)
    out.write(f"shares_species_with_MAG\t{ng_sh}\t{sh_p}\t{sh_d}\t{frac(sh_d, sh_p)}\n")
    out.write(f"no_shared_species\t{len(n_prot) - ng_sh}\t{ns_p}\t{ns_d}\t{frac(ns_d, ns_p)}\n")
    out.write(f"hypotheticals_only\t{len(n_prot)}\t{tot_h}\t{tot_hd}\t{frac(tot_hd, tot_h)}\n")
    # CORRECTION 2026-09-17: this was pseudo-replicated. Iterating over cultured genomes counts a
    # MAG representative once per cultured genome that best-matches it, so the shared genomes
    # collapse onto far fewer distinct representatives and the MAG side was weighted by how many
    # strains of a species happened to be cultured. The honest n is the number of distinct
    # representatives. Both are written: _pseudorep is the old figure, kept so the superseded
    # number stays traceable.
    mm_old = [float(mag_dark[shared_mag[g][0]]) for g in n_prot
              if g in shared_mag and shared_mag[g][0] in mag_dark]
    mm_old.sort()
    med_old = mm_old[len(mm_old) // 2] if mm_old else float("nan")
    out.write(f"matched_MAG_reps_pseudorep\t{len(mm_old)}\tNA\tNA\t{med_old:.6f}\n")

    reps = {shared_mag[g][0] for g in n_prot if g in shared_mag and shared_mag[g][0] in mag_dark}
    mm = sorted(float(mag_dark[r]) for r in reps)
    med = mm[len(mm) // 2] if mm else float("nan")
    rp = sum(int(mag_nprot[r]) for r in reps if r in mag_nprot)
    rd = sum(int(mag_ndark[r]) for r in reps if r in mag_ndark)
    out.write(f"matched_MAG_reps_distinct\t{len(reps)}\t{rp}\t{rd}\t{frac(rd, rp)}\n")
    out.write(f"matched_MAG_reps_distinct_median\t{len(reps)}\tNA\tNA\t{med:.6f}\n")

    # Paired comparison, the matched unit. Both sides were pseudo-replicated: several cultured
    # strains map to one MAG representative, so pooling either side weights species by how many
    # strains happened to be cultured. Here each shared species contributes one pair: its cultured
    # genomes aggregated, against its MAG representative.
    cp, cd = defaultdict(int), defaultdict(int)
    for g in n_prot:
        if g in shared_mag and shared_mag[g][0] in mag_dark:
            r = shared_mag[g][0]
            cp[r] += n_prot[g]
            cd[r] += n_dark[g]
    with open(f"{P}/results/step3_paired_by_species.tsv", "w") as pf:
        pf.write("mag_rep\tn_cultured_genomes\tcultured_n_proteins\tcultured_n_dark\t"
                 "cultured_dark_fraction\tmag_n_proteins\tmag_n_dark\tmag_dark_fraction\tdelta\n")
        diffs = []
        by_depth = []
        for r in sorted(cp):
            cf = cd[r] / cp[r] if cp[r] else float("nan")
            mf = float(mag_dark[r])
            ng = sum(1 for g in n_prot if g in shared_mag and shared_mag[g][0] == r)
            diffs.append(cf - mf)
            by_depth.append((ng, cf - mf))
            pf.write(f"{r}\t{ng}\t{cp[r]}\t{cd[r]}\t{cf:.6f}\t{mag_nprot.get(r,'NA')}\t"
                     f"{mag_ndark.get(r,'NA')}\t{mf:.6f}\t{cf - mf:+.6f}\n")
    n_lower = sum(1 for d in diffs if d < 0)
    diffs.sort()
    # median = the upper middle element for an even n, the convention used throughout this file
    med_d = diffs[len(diffs) // 2] if diffs else float("nan")
    out.write(f"paired_by_species_cultured_lower\t{n_lower}\t{len(diffs)}\tNA\t{med_d:+.6f}\n")
    # Significance of the paired difference. Wilcoxon signed rank on the deltas (two sided) and
    # the sign test, which drops the magnitudes and only asks how often cultured is the lower one.
    wil_p = wilcoxon(diffs).pvalue if len(diffs) > 1 else float("nan")
    sign_p = binomtest(n_lower, len(diffs), 0.5).pvalue if diffs else float("nan")
    out.write(f"paired_by_species_wilcoxon_p\t{len(diffs)}\tNA\tNA\t{wil_p:.4g}\n")
    out.write(f"paired_by_species_sign_test_p\t{len(diffs)}\tNA\tNA\t{sign_p:.4g}\n")
    # Cultivation depth: a species cultured as many strains is one the collection sampled deeply,
    # a species cultured once is one it caught by chance.
    for label, sel in (("ge10_strains", lambda n: n >= 10), ("single_strain", lambda n: n == 1)):
        st = sorted(d for n, d in by_depth if sel(n))
        med_s = st[len(st) // 2] if st else float("nan")
        out.write(f"paired_by_species_{label}_median\t{len(st)}\tNA\tNA\t{med_s:+.6f}\n")
        print(f"[info] paired by species, {label}: n={len(st)} median delta {med_s:+.6f}")
    print(f"[info] paired by species: cultured darker in {len(diffs)-n_lower}/{len(diffs)}, "
          f"median delta {med_d:+.6f}, wilcoxon p {wil_p:.4g}, sign test p {sign_p:.4g}")
    # inverse view: the same no-UniRef50 protein seen twice, once in a closed cultured genome
    # and once as a MAG-dark rep. Restricted to the shared-species genomes, where a MAG
    # counterpart exists at all.
    sh_dark_ids = {i for i in dark_ids_set if i.split("@")[0] in shared_mag}
    rec = len(sh_dark_ids & cult_subject_hit)
    out.write(f"shared_genome_dark_proteins_recovered_as_MAG_dark\t{ng_sh}\t"
              f"{len(sh_dark_ids)}\t{rec}\t{frac(rec, len(sh_dark_ids))}\n")
print(f"[info] cultured dark rate overall {tot_d}/{tot_p} = {tot_d / tot_p:.4f}")
print("DONE_STEP3_04")
