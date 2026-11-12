#!/usr/bin/env python3
"""Step 1: representatives, mature sequences, conserved neighbours, complex pairs, controls.

Families: MMseqs2 id30/c80 clusters over the 3,617,011 proteins of 1,957 genomes.
Block = named genus, else ANI95 group (as in analysis2_families/family_sharing.py).
All outputs go to --out (analysis4_structures/work/step1).
"""
import argparse, json, math, re, sys, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sps

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--seed", type=int, default=20261111)
ap.add_argument("--min-controls", type=int, default=120)
ap.add_argument("--max-controls", type=int, default=160)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
R = a.root
rng = np.random.default_rng(a.seed)
EXCLUDED = {"External_Comparators_CrossResource_Duplicate_54", "SPIRE_ProGenomes3_SpecI_Excluded_148"}
MODULES = [  # paper modules: (label, focal gene, partner gene) in harmonised gene identifiers
    ("paper_module_PilJ", "1PRr2__self_SemiBin_81::1PRr2_468829_3", "1PRr2__self_SemiBin_81::1PRr2_468829_2"),
    ("paper_module_envelope", "V47PRr1__self_SemiBin_699_sub::V47PRr1_105096_18",
     "V47PRr1__self_SemiBin_699_sub::V47PRr1_105096_19"),
]
log = lambda *x: print(*x, flush=True)

# ---------------- load ----------------
gm = pd.read_csv(R / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
gm = gm[~gm.cohort.isin(EXCLUDED)].copy()
ani = pd.read_csv(R / "source-metadata/ani_groups/ani_95_groups_1957.tsv", sep="\t", dtype=str)
gm = gm.merge(ani[["genome_id", "ani95_group_id"]], on="genome_id", how="left", validate="1:1")
assert not gm.ani95_group_id.isna().any()
gm["block"] = np.where(gm.tax_genus.fillna("g__").isin(["g__", "", "missing", "none"]),
                       "ani::" + gm.ani95_group_id, gm.tax_genus)
gm["is_eq"] = gm.cohort == "EQ_Thermoplasmatota_237"
gm = gm.reset_index(drop=True)
assert len(gm) == 1957 and gm.is_eq.sum() == 237, (len(gm), gm.is_eq.sum())

g = pd.read_csv(R / "genecall/harmonized_gene_manifest.tsv.gz", sep="\t",
                usecols=["gene_id", "genome_id", "contig_id", "start", "end", "strand", "contig_gene_rank",
                         "contig_n_genes", "partial_flag", "n_internal_stops", "is_complete_orf",
                         "contig_edge_status", "aa_length"], dtype={"partial_flag": str})
g = g[g.genome_id.isin(gm.genome_id)]
assert len(g) == 3617011, len(g)
sp = pd.read_csv(R / "signalp/harmonized_signalp6_predictions.tsv.gz", sep="\t",
                 usecols=["gene_id", "signalp_class", "sp_any", "sp_prob", "cleavage_site_pos"])
uni = pd.read_csv(R / "mmseqs/harmonized_uniref50_per_gene.tsv.gz", sep="\t",
                  usecols=["gene_id", "primary_direct_uniref50_dark", "direct_uniref50_best_target"])
g = g.merge(sp, on="gene_id", validate="1:1").merge(uni, on="gene_id", validate="1:1")
assert len(g) == 3617011
fam = pd.read_csv(R / "analysis2_families/clusters/families_id30_c80.tsv", sep="\t", header=None,
                  names=["fam", "gene_id"])
assert fam.gene_id.is_unique
g = g.merge(fam, on="gene_id", how="left", validate="1:1")
assert not g.fam.isna().any()
g["dark"] = g.primary_direct_uniref50_dark.astype(bool)
g["sp"] = g.sp_any.astype(bool)
g["complete"] = g.is_complete_orf.astype(bool) & (g.n_internal_stops == 0)
g["interior"] = g.complete & (g.contig_edge_status == "interior")
cs = g.cleavage_site_pos.astype(str).str.extract(r"CS pos: (\d+)-(\d+)")
g["cut"] = pd.to_numeric(cs[0], errors="coerce").fillna(0).astype(int)   # residues removed
g.loc[~g.sp, "cut"] = 0
gi = pd.Index(gm.genome_id)
g["gix"] = gi.get_indexer(g.genome_id)
g = g.sort_values(["gix", "contig_id", "contig_gene_rank"], kind="stable").reset_index(drop=True)
g["ctg"] = pd.factorize(g.genome_id + "\t" + g.contig_id)[0]
g["fc"], fam_names = pd.factorize(g.fam)
n = len(g)
ctg = g.ctg.to_numpy(); rank = g.contig_gene_rank.to_numpy(); fc = g.fc.to_numpy(); gix = g.gix.to_numpy()
strand = g.strand.to_numpy(); st = g.start.to_numpy(); en = g.end.to_numpy()
is_eq_gene = gm.is_eq.to_numpy()[gix]
blk_codes, blk_names = pd.factorize(gm.block)
gblk = blk_codes[gix]
same = ctg[1:] == ctg[:-1]
d = np.diff(rank)[same]
log("rank steps within contigs:", dict(zip(*np.unique(d, return_counts=True))))
assert (d == 1).all(), "gene ranks are not consecutive within contigs"
gid = g.gene_id.to_numpy()
gpos = pd.Series(np.arange(n), index=gid)
length = g.aa_length.to_numpy(); cut = g.cut.to_numpy()
mature_len = length - cut
complete = g.complete.to_numpy(); interior = g.interior.to_numpy(); spb = g.sp.to_numpy()
spc = g.signalp_class.to_numpy()

# sequences of EQ genes and of any gene needed later (read lazily: full FASTA index)
def read_fasta_subset(path, wanted):
    out, keep, name, buf = {}, False, None, []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if keep: out[name] = "".join(buf)
                name = line[1:].split()[0]; keep = name in wanted; buf = []
            elif keep:
                buf.append(line.strip())
    if keep: out[name] = "".join(buf)
    return out

# ---------------- focal families ----------------
ft = pd.read_csv(R / "analysis2_families/results/dark_exported_widespread_families_id30_c80.tsv", sep="\t")
assert len(ft) == 981 and ft.family_rep.is_unique
fam_code = pd.Series(np.arange(len(fam_names)), index=fam_names)
focal_codes = fam_code[ft.family_rep].to_numpy()
mod_rows = []
for lab, fg, pg in MODULES:
    ok = fg in gpos.index and pg in gpos.index
    row = {"module": lab, "focal_gene": fg, "partner_gene": pg, "mapped_to_harmonised_gene_calls": ok}
    if ok:
        i, j = gpos[fg], gpos[pg]
        row.update(focal_family=fam_names[fc[i]], partner_family=fam_names[fc[j]],
                   focal_family_in_981=bool(fc[i] in set(focal_codes)),
                   partner_family_in_981=bool(fc[j] in set(focal_codes)),
                   adjacent=bool(ctg[i] == ctg[j] and abs(rank[i] - rank[j]) == 1),
                   focal_signalp=spc[i], partner_signalp=spc[j], focal_dark=bool(g.dark[i]),
                   partner_dark=bool(g.dark[j]), focal_aa=int(length[i]), partner_aa=int(length[j]))
    mod_rows.append(row)
mods = pd.DataFrame(mod_rows)
log(mods.to_string())
extra_focal = [fam_code[r.focal_family] for r in mods.itertuples() if r.mapped_to_harmonised_gene_calls]
analysis_focal = np.unique(np.concatenate([focal_codes, np.array(extra_focal, dtype=int)]))
is_focal981 = np.zeros(len(fam_names), bool); is_focal981[focal_codes] = True

# family level stats (all members, EQ members)
def fam_stats(mask, suffix):
    x = g[mask].groupby("fc").agg(n=("gene_id", "size"), dark=("dark", "mean"), sp=("sp", "mean"),
                                  med_len=("aa_length", "median"), n_genomes=("gix", "nunique"))
    x.columns = [c + suffix for c in x.columns]
    return x
fs = fam_stats(np.ones(n, bool), "_all").join(fam_stats(is_eq_gene, "_eq"), how="left")
fs["n_eq"] = fs.n_eq.fillna(0).astype(int)
def cls(dk, s):
    return (pd.Series(np.where(np.asarray(dk) > 0.5, "dark", "nondark"), dtype=object) + "_"
            + pd.Series(np.where(np.asarray(s) > 0.5, "exported", "not_exported"), dtype=object)).to_numpy()
fs["class_eq_members"] = np.where(fs.n_eq > 0, cls(fs.dark_eq.fillna(0), fs.sp_eq.fillna(0)), "no_eq_member")
fs["class_all_members"] = cls(fs.dark_all, fs.sp_all)

# ---------------- representatives ----------------
eq_focal = g[is_eq_gene & np.isin(fc, analysis_focal)]
rep_rows, queries = [], {}
def add_query(i, role):
    key = (gid[i], int(cut[i]))
    q = queries.setdefault(key, {"gene_id": gid[i], "cut": int(cut[i]), "roles": [], "ix": int(i)})
    q["roles"].append(role)
    return key
rep_of = {}
for f, sub in eq_focal.groupby("fc"):
    med = float(sub.aa_length.median())
    spm = sub[sub.sp]
    modal = spm.signalp_class.value_counts().index[0] if len(spm) else "OTHER"
    sub = sub.assign(dev=(sub.aa_length - med).abs(), is_modal=sub.signalp_class == modal, has_cs=sub.cut > 0)
    tiers = [("T1_complete_interior_modalSP_cleavage", sub.interior & sub.is_modal & sub.has_cs),
             ("T2_complete_edge_modalSP_cleavage", sub.complete & sub.is_modal & sub.has_cs),
             ("T3_complete_anySP_cleavage", sub.complete & sub.has_cs),
             ("T4_complete_no_cleavage_full_length_used", sub.complete),
             ("T5_partial_ORF", np.ones(len(sub), bool))]
    for tname, m in tiers:
        c = sub[m]
        if len(c):
            c = c.sort_values(["dev", "sp_prob", "gene_id"], ascending=[True, False, True])
            r = c.iloc[0]
            break
    i = gpos[r.gene_id]
    rep_of[f] = i
    key = add_query(i, "monomer:" + fam_names[f]) if tname[:2] in ("T1", "T2", "T3", "T4") else None
    rep_rows.append({"family_rep": fam_names[f], "in_981": bool(is_focal981[f]), "rep_gene": r.gene_id,
                     "tier": tname, "usable": key is not None, "n_eq_members": len(sub), "median_eq_aa": med,
                     "modal_sp_class": modal, "rep_sp_class": r.signalp_class, "rep_sp_prob": r.sp_prob,
                     "rep_cleavage": r.cleavage_site_pos, "rep_aa": int(r.aa_length), "signal_peptide_cut": int(r.cut),
                     "mature_aa": int(r.aa_length - r.cut), "rep_complete": bool(r.complete),
                     "rep_interior": bool(r.interior), "rep_partial_flag": r.partial_flag,
                     "n_eq_T1_candidates": int(tiers[0][1].sum())})
reps = pd.DataFrame(rep_rows)
log("representative tiers (981):\n", reps[reps.in_981].tier.value_counts().to_string())

# ---------------- neighbours ----------------
focal_gene_ix = np.flatnonzero(np.isin(fc, analysis_focal))
rows = []
has = {}
for k in (-2, -1, 1, 2):
    j = focal_gene_ix + k
    ok = (j >= 0) & (j < n)
    ok[ok] &= ctg[j[ok]] == ctg[focal_gene_ix[ok]]
    has[k] = ok
    i_, j_ = focal_gene_ix[ok], j[ok]
    if k > 0:
        gap = st[j_] - en[i_] - 1
    else:
        gap = st[i_] - en[j_] - 1
    # intervening gene for |k|==2 must not matter; gap is only reported for |k|==1
    rows.append(pd.DataFrame({"i": i_, "j": j_, "k": k, "F": fc[i_], "P": fc[j_], "gix": gix[i_],
                              "blk": gblk[i_], "is_eq": is_eq_gene[i_], "same_strand": strand[i_] == strand[j_],
                              "gap": gap}))
nb = pd.concat(rows, ignore_index=True)
loc = pd.DataFrame({"i": focal_gene_ix, "F": fc[focal_gene_ix], "is_eq": is_eq_gene[focal_gene_ix],
                    "usable1": has[-1] | has[1], "both1": has[-1] & has[1],
                    "usable2": has[-2] | has[-1] | has[1] | has[2]})

def summarise(nbs, locs, tag):
    den = locs.groupby("F").agg(n_loci=("i", "size"), n_usable=("usable1", "sum"), n_both=("both1", "sum"),
                                n_usable2=("usable2", "sum"))
    a1 = nbs[(nbs.k.abs() == 1) & (nbs.F != nbs.P)]
    s1 = a1.groupby(["F", "P"]).agg(n_adj_loci=("i", "nunique"), n_adj_genomes=("gix", "nunique"),
                                    n_adj_blocks=("blk", "nunique"), frac_same_strand=("same_strand", "mean"),
                                    median_gap_bp=("gap", "median"), min_gap_bp=("gap", "min"),
                                    max_gap_bp=("gap", "max"))
    a2 = nbs[(nbs.F != nbs.P) & nbs.same_strand]
    s2 = a2.groupby(["F", "P"]).agg(n_w2ss_loci=("i", "nunique"), n_w2ss_genomes=("gix", "nunique"),
                                    n_w2ss_blocks=("blk", "nunique"))
    s = s2.join(s1, how="outer").reset_index()
    s = s.merge(den.reset_index(), on="F", how="left")
    for c in ["n_adj_loci", "n_adj_genomes", "n_adj_blocks", "n_w2ss_loci", "n_w2ss_genomes", "n_w2ss_blocks"]:
        s[c] = s[c].fillna(0).astype(int)
    s["frac_adj"] = s.n_adj_loci / s.n_usable.clip(lower=1)
    s["frac_w2ss"] = s.n_w2ss_loci / s.n_usable2.clip(lower=1)
    s["conserved_adj"] = (s.frac_adj >= 0.5) & (s.n_adj_genomes >= 5) & (s.n_adj_blocks >= 3)
    s["conserved_w2ss"] = (s.frac_w2ss >= 0.5) & (s.n_w2ss_genomes >= 5) & (s.n_w2ss_blocks >= 3)
    s = s.set_index(["F", "P"])
    s.columns = [c + "_" + tag for c in s.columns]
    return s
S_all = summarise(nb, loc, "all")
S_eq = summarise(nb[nb.is_eq], loc[loc.is_eq], "eq")
S = S_all.join(S_eq, how="outer").reset_index()
for c in S.columns:
    if c.startswith("conserved"):
        S[c] = S[c].fillna(False).astype(bool)
S["conserved_adj_any"] = S.conserved_adj_all | S.conserved_adj_eq
S["conserved_w2ss_any"] = S.conserved_w2ss_all | S.conserved_w2ss_eq
# tandem same-family neighbours (reported only)
tand = nb[(nb.k.abs() == 1) & (nb.F == nb.P)].groupby("F").agg(n_tandem_loci=("i", "nunique"),
                                                                n_tandem_genomes=("gix", "nunique"))
S["focal_family"] = fam_names[S.F.to_numpy()]
S["partner_family"] = fam_names[S.P.to_numpy()]
S["focal_in_981"] = is_focal981[S.F.to_numpy()]
S["partner_in_981"] = is_focal981[S.P.to_numpy()]
part = fs.loc[S.P.to_numpy()].reset_index(drop=True)
for c in ["n_all", "n_eq", "dark_all", "sp_all", "dark_eq", "sp_eq", "med_len_all", "class_eq_members",
          "class_all_members"]:
    S["partner_" + c] = part[c].to_numpy()
# cheap annotation: most frequent UniRef50 target among partner members
keepP = np.unique(S.P[S.conserved_adj_any | S.conserved_w2ss_any].to_numpy())
ann = {}
subP = g[np.isin(fc, keepP) & (g.direct_uniref50_best_target != "none")]
for p, x in subP.groupby("fc"):
    vc = x.direct_uniref50_best_target.value_counts()
    ann[p] = (vc.index[0], int(vc.iloc[0]), int(len(x)))
S["partner_top_uniref50"] = [ann.get(p, ("none", 0, 0))[0] for p in S.P]
S["partner_top_uniref50_n"] = [ann.get(p, ("none", 0, 0))[1] for p in S.P]
S["partner_n_members_with_uniref50_hit"] = [ann.get(p, ("none", 0, 0))[2] for p in S.P]
keep = (S.conserved_adj_any | S.conserved_w2ss_any | (S.n_adj_genomes_all.fillna(0) >= 5))
S[keep].drop(columns=["F", "P"]).to_csv(a.out / "neighbour_family_pairs.tsv", sep="\t", index=False)

foc = S[S.focal_in_981]
summ = {
    "n_focal": 981,
    "focal_with_conserved_adjacent_partner_all_genomes": int(foc[foc.conserved_adj_all].F.nunique()),
    "focal_with_conserved_adjacent_partner_eq_only": int(foc[foc.conserved_adj_eq].F.nunique()),
    "focal_with_conserved_adjacent_partner_any": int(foc[foc.conserved_adj_any].F.nunique()),
    "pairs_conserved_adjacent_all_genomes": int(foc.conserved_adj_all.sum()),
    "pairs_conserved_adjacent_eq_only": int(foc.conserved_adj_eq.sum()),
    "pairs_conserved_adjacent_any": int(foc.conserved_adj_any.sum()),
    "focal_with_conserved_within2_same_strand_any": int(foc[foc.conserved_w2ss_any].F.nunique()),
    "focal_with_within2_same_strand_but_not_adjacent_conserved": int(
        len(set(foc[foc.conserved_w2ss_any].F) - set(foc[foc.conserved_adj_any].F))),
    "pairs_conserved_within2_same_strand_any": int(foc.conserved_w2ss_any.sum()),
    "focal_with_tandem_same_family_ge5_genomes": int((tand.reindex(focal_codes).n_tandem_genomes.fillna(0) >= 5).sum()),
    "partner_class_eq_members_of_conserved_adjacent_pairs": foc[foc.conserved_adj_any].partner_class_eq_members.value_counts().to_dict(),
}
log(json.dumps(summ, indent=1))

# ---------------- complex pair set ----------------
pair_rows = []
seen = set()
real = S[(S.focal_in_981 & S.conserved_adj_any)].copy()
real["source"] = "conserved_adjacent"
forced = []
for r in mods.itertuples():
    if not r.mapped_to_harmonised_gene_calls: continue
    F, P = fam_code[r.focal_family], fam_code[r.partner_family]
    m = S[(S.F == F) & (S.P == P)]
    forced.append((r.module, F, P, gpos[r.focal_gene], gpos[r.partner_gene], m))
a1 = nb[(nb.k.abs() == 1)]
a1_by = {k: v for k, v in a1.groupby(["F", "P"])}
med_all = fs.med_len_all
def choose_locus(F, P):
    x = a1_by[(F, P)]
    xe = x[x.is_eq] if x.is_eq.any() else x
    i, j = xe.i.to_numpy(), xe.j.to_numpy()
    modal = reps.loc[reps.family_rep == fam_names[F], "modal_sp_class"]
    modal = modal.iloc[0] if len(modal) else None
    dev = np.abs(length[i] - med_all[F]) / med_all[F] + np.abs(length[j] - med_all[P]) / med_all[P]
    key = pd.DataFrame({"i": i, "j": j, "c": -(complete[i] & complete[j]).astype(int),
                        "r": -(i == rep_of.get(F, -1)).astype(int),
                        "m": -((spc[i] == modal) & (cut[i] > 0)).astype(int) if modal else 0,
                        "n": -(interior[i] & interior[j]).astype(int), "dev": dev, "g": gid[i]})
    key = key.sort_values(["c", "r", "m", "n", "dev", "g"])
    return int(key.i.iloc[0]), int(key.j.iloc[0]), bool(xe.is_eq.iloc[0])
def add_pair(F, P, i, j, source, srow, pid):
    ka, kb = add_query(i, f"{pid}:A"), add_query(j, f"{pid}:B")
    d_ = {"pair_id": pid, "type": "real", "source": source, "focal_family": fam_names[F],
          "partner_family": fam_names[P], "focal_gene": gid[i], "partner_gene": gid[j],
          "locus_in_eq": bool(is_eq_gene[i]), "focal_is_monomer_rep": bool(i == rep_of.get(F, -1)),
          "focal_cut": int(cut[i]), "partner_cut": int(cut[j]), "len_A": int(mature_len[i]),
          "len_B": int(mature_len[j]), "focal_complete": bool(complete[i]), "partner_complete": bool(complete[j]),
          "partner_sp_class": spc[j], "locus_same_strand": bool(strand[i] == strand[j])}
    if srow is not None and len(srow):
        s0 = srow.iloc[0]
        for c in ["n_adj_loci_all", "n_adj_genomes_all", "n_adj_blocks_all", "frac_adj_all", "n_usable_all",
                  "n_adj_loci_eq", "n_adj_genomes_eq", "n_adj_blocks_eq", "frac_adj_eq", "n_usable_eq",
                  "conserved_adj_all", "conserved_adj_eq", "frac_same_strand_all", "median_gap_bp_all",
                  "partner_class_eq_members", "partner_class_all_members", "partner_n_all", "partner_in_981",
                  "partner_top_uniref50"]:
            d_[c] = s0[c]
    pair_rows.append(d_)
k_ = 0
for r in real.sort_values(["n_adj_genomes_all"], ascending=False).itertuples():
    und = (min(r.F, r.P), max(r.F, r.P))
    if und in seen:   # both members focal and reciprocal: keep one orientation
        for d_ in pair_rows:
            if {d_["focal_family"], d_["partner_family"]} == {fam_names[r.F], fam_names[r.P]}:
                d_["reciprocal_focal_pair"] = True
        continue
    seen.add(und)
    k_ += 1
    i, j, _ = choose_locus(r.F, r.P)
    add_pair(r.F, r.P, i, j, "conserved_adjacent", S[(S.F == r.F) & (S.P == r.P)], f"R{k_:04d}")
for lab, F, P, i, j, m in forced:
    und = (min(F, P), max(F, P))
    k_ += 1
    note = lab + ("|also_meets_conservation_criteria" if und in seen else "|forced_not_meeting_or_not_focal")
    add_pair(F, P, i, j, note, m, f"R{k_:04d}")
    seen.add(und)
n_real = len(pair_rows)
log("real pairs:", n_real)

# ---------------- controls ----------------
M = sps.csr_matrix((np.ones(n, np.int32), (fc, gix)), shape=(len(fam_names), len(gm)))
M.data[:] = 1; M.sum_duplicates(); M.data[:] = 1
real_partners = {}
for d_ in pair_rows:
    real_partners.setdefault(d_["focal_family"], set()).add(d_["partner_family"])
by_genome = {k: v.to_numpy() for k, v in pd.Series(np.arange(n)).groupby(gix)}
_af = set(analysis_focal.tolist())
fam_idx = {k: np.flatnonzero(fc == k) for k in _af}
n_all = fs.n_all
per = max(1, math.ceil(a.min_controls / max(n_real, 1)))
ctrl_rows, ck = [], 0
LEVELS = [(0.20, 2.0, True), (0.20, 2.0, False), (0.35, 3.0, False), (0.50, 5.0, False)]
ctrl_for = list(pair_rows)
if per == 1 and n_real > a.max_controls:   # random subset of real pairs, forced module pairs always included
    forced_ids = [d_ for d_ in pair_rows if d_["source"] != "conserved_adjacent"]
    rest = [d_ for d_ in pair_rows if d_["source"] == "conserved_adjacent"]
    sel = rng.choice(len(rest), size=a.max_controls - len(forced_ids), replace=False)
    ctrl_for = [rest[k] for k in sorted(sel)] + forced_ids
for d_ in ctrl_for:
    F, P = fam_code[d_["focal_family"]], fam_code[d_["partner_family"]]
    i = gpos[d_["focal_gene"]]
    Fi = fam_idx.get(F)
    if Fi is None: Fi = np.flatnonzero(fc == F)
    near = set()
    for k in range(-10, 11):
        j = Fi + k
        ok = (j >= 0) & (j < n)
        ok[ok] &= ctg[j[ok]] == ctg[Fi[ok]]
        near.update(fc[j[ok]].tolist())
    cand = by_genome[gix[i]]
    cand = cand[complete[cand] & (mature_len[cand] >= 30)]
    cf = fc[cand]
    bad = near | {F} | {fam_code[x] for x in real_partners.get(d_["focal_family"], ())}
    cand = cand[~np.isin(cf, list(bad))]
    cf = fc[cand]
    co = np.asarray((M[cf] @ M[F].T).todense()).ravel()
    target_co = max(float(np.nan_to_num(pd.to_numeric(d_.get("n_adj_genomes_all", 0), errors="coerce"))), 1.0)
    LB, nP = d_["len_B"], float(n_all[P])
    ldev = np.abs(mature_len[cand] - LB) / LB
    fdev = np.abs(length[cand] - med_all.to_numpy()[cf]) / med_all.to_numpy()[cf]
    used = set()
    for rep_ in range(per):
        pick = None
        for lv, (lt, ft_, need_sp) in enumerate(LEVELS):
            m = (ldev <= lt) & (fdev <= 0.2) & (co >= target_co / ft_) & (co <= target_co * ft_)
            m &= (n_all.to_numpy()[cf] >= nP / ft_) & (n_all.to_numpy()[cf] <= nP * ft_)
            if need_sp: m &= (spb[cand] == bool(d_["partner_cut"] > 0))
            m &= ~np.isin(cf, list(used))
            if m.any():
                pick = int(rng.choice(np.flatnonzero(m))); break
        if pick is None: 
            ctrl_rows.append({"matched_real_pair": d_["pair_id"], "status": "no_matched_control_found"})
            continue
        j = int(cand[pick]); used.add(cf[pick]); ck += 1
        pid = f"C{ck:04d}"
        ka, kb = add_query(i, f"{pid}:A"), add_query(j, f"{pid}:B")
        ctrl_rows.append({"pair_id": pid, "type": "control", "status": "ok", "matched_real_pair": d_["pair_id"],
                          "match_level": lv, "focal_family": d_["focal_family"], "partner_family": fam_names[cf[pick]],
                          "focal_gene": gid[i], "partner_gene": gid[j], "focal_cut": int(cut[i]),
                          "partner_cut": int(cut[j]), "len_A": int(mature_len[i]), "len_B": int(mature_len[j]),
                          "real_len_B": LB, "cooccurring_genomes": int(co[pick]), "real_adjacent_genomes": target_co,
                          "control_family_n_all": int(n_all.to_numpy()[cf[pick]]), "real_partner_n_all": nP,
                          "partner_sp_class": spc[j], "same_contig_as_focal": bool(ctg[i] == ctg[j]),
                          "gene_rank_distance_if_same_contig": int(abs(rank[i] - rank[j])) if ctg[i] == ctg[j] else -1,
                          "control_class_eq_members": fs.class_eq_members.to_numpy()[cf[pick]]})
ctrl = pd.DataFrame(ctrl_rows)
log("controls:", ctrl.status.value_counts().to_dict(), "per real pair:", per)

# ---------------- write ----------------
wanted = {q["gene_id"] for q in queries.values()}
seqs = read_fasta_subset(R / "genecall/harmonized_proteins.faa", wanted)
assert wanted <= set(seqs)
qrows = []
with open(a.out / "queries.faa", "w") as f:
    for n_, (key, q) in enumerate(sorted(queries.items()), 1):
        s = seqs[q["gene_id"]].rstrip("*")
        assert "*" not in s and len(s) == length[q["ix"]], (q["gene_id"], len(s), length[q["ix"]])
        m = s[q["cut"]:]
        q["qid"] = f"Q{n_:05d}"
        f.write(f">{q['qid']}\n{m}\n")
        ix = q["ix"]
        qrows.append({"qid": q["qid"], "gene_id": q["gene_id"], "family": fam_names[fc[ix]], "cut": q["cut"],
                      "query_aa": len(m), "full_aa": len(s), "signalp_class": spc[ix], "is_eq": bool(is_eq_gene[ix]),
                      "sha256_query": hashlib.sha256(m.encode()).hexdigest(), "roles": ";".join(q["roles"])})
qt = pd.DataFrame(qrows)
qt.to_csv(a.out / "queries.tsv", sep="\t", index=False)
qid_of = {k: q["qid"] for k, q in queries.items()}
reps["qid"] = [qid_of.get((r.rep_gene, r.signal_peptide_cut), "") if r.usable else "" for r in reps.itertuples()]
reps.to_csv(a.out / "representatives.tsv", sep="\t", index=False)
pairs = pd.concat([pd.DataFrame(pair_rows), ctrl[ctrl.status == "ok"]] if len(ctrl) else [pd.DataFrame(pair_rows)],
                  ignore_index=True)
pairs["qid_A"] = [qid_of[(r.focal_gene, r.focal_cut)] for r in pairs.itertuples()]
pairs["qid_B"] = [qid_of[(r.partner_gene, r.partner_cut)] for r in pairs.itertuples()]
pairs.to_csv(a.out / "pairs.tsv", sep="\t", index=False)
ctrl.to_csv(a.out / "control_selection.tsv", sep="\t", index=False)
mods.to_csv(a.out / "paper_modules_mapping.tsv", sep="\t", index=False)
fs_out = fs.copy(); fs_out.insert(0, "family", fam_names[fs.index.to_numpy()])
need = set(pairs.focal_family) | set(pairs.partner_family) | set(ft.family_rep)
fs_out[fs_out.family.isin(need)].to_csv(a.out / "family_stats.tsv", sep="\t", index=False)
pd.DataFrame({"gene_id": gid, "gix": gix.astype(np.int32), "ctg": ctg.astype(np.int32), "rank": rank.astype(np.int32),
              "fc": fc.astype(np.int32), "is_eq": is_eq_gene}).to_pickle(a.out / "genes_min.pkl")
pd.Series(fam_names).to_csv(a.out / "family_codes.tsv", sep="\t", header=False)
summ.update(n_real_pairs=n_real, n_controls=int((ctrl.status == "ok").sum()) if len(ctrl) else 0,
            n_queries=len(qt), n_usable_representatives_981=int(reps[reps.in_981].usable.sum()),
            representative_tiers_981=reps[reps.in_981].tier.value_counts().to_dict(), seed=a.seed,
            controls_per_real_pair=per)
(a.out / "step1_summary.json").write_text(json.dumps(summ, indent=1, default=str))
log("STEP1_DONE", json.dumps(summ, default=str))
